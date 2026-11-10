from pathlib import Path
import argparse
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Environement import Environment
from frames import (
    body_to_ecef_matrix, quaternion_product, quaternion_conjugate,
    quaternion_from_rotation_vector, quaternion_to_rotation_vector,
)
from kalman import NavigationKF, misalignment_matrix
from sensors import GNSS, GNSSSpec, StarTracker, StarTrackerSpec


def angle_error(estimate, truth):
    return np.linalg.norm(quaternion_to_rotation_vector(quaternion_product(
        quaternion_conjugate(estimate), truth)))


def accel_scale_case(seed=11, duration=180., dt=.05):
    rng = np.random.default_rng(seed)
    scale = np.array([.005, -.003, .002])
    bias = np.array([.01, -.02, .005])
    k = NavigationKF(np.zeros(3), np.zeros(3),
        np.r_[np.full(3, .1), np.full(3, .05), np.full(3, .1)],
        gravity_ecef=np.zeros(3), accel_psd=.001**2*dt,
        estimate_accel_scale=True, accel_sf_sigma0=.02)
    def force(t):
        return np.array([3 + 2*np.sin(.1*t), 2*np.cos(.13*t), 1 + 2*np.sin(.17*t)])
    p, v = np.zeros(3), np.zeros(3)
    last = force(0.)
    k.prime_imu((1 + scale)*last + bias + rng.normal(0, .001, 3), np.zeros(3))
    for i in range(1, round(duration/dt) + 1):
        t = i*dt
        f = force(t)
        a = (last + f)/2
        p += v*dt + .5*a*dt**2
        v += a*dt
        k.predict((1 + scale)*f + bias + rng.normal(0, .001, 3), np.eye(3), dt)
        if i % 10 == 0:
            k.update_gnss(p + rng.normal(0, .01, 3), v + rng.normal(0, .005, 3),
                np.eye(3), np.zeros(3), np.zeros(3), np.full(3, .01), np.full(3, .005))
        last = f
    return {
        'seed':seed, 'duration_s':duration, 'imu_rate_hz':1/dt,
        'true_scale':scale.tolist(), 'estimated_scale':k.accel_sf.tolist(),
        'scale_error_norm':float(np.linalg.norm(k.accel_sf - scale)),
        'true_bias_mps2':bias.tolist(), 'estimated_bias_mps2':k.x[6:9].tolist(),
        'bias_error_norm_mps2':float(np.linalg.norm(k.x[6:9] - bias)),
        'scale_tolerance':4e-4, 'bias_tolerance_mps2':.002,
    }


def attitude_case(seed=12, duration=450., dt=.02):
    rng = np.random.default_rng(seed)
    env = Environment()
    p = np.array([0., 0., env.earth_radius_m + 150000.])
    launch = quaternion_from_rotation_vector([0., -np.pi/2, 0.])
    scale = np.array([.002, -.0015, .001])
    mis = np.array([.001, -.0008, .0006, .0012, -.001, .0005])
    bias = np.array([.0005, -.0003, .0002])
    gain = np.eye(3) + np.diag(scale) + misalignment_matrix(mis)
    spec = StarTrackerSpec(boresight_body=np.array([1., 0., 0.]))
    star = StarTracker(spec, env.earth_rotation_ecef_radps, seed=seed + 7)
    def truth(t):
        angles = np.array([.3, .4, .35])*np.sin(np.array([.031, .047, .063])*t)
        rates = np.array([.3, .4, .35])*np.array([.031, .047, .063])*np.cos(np.array([.031, .047, .063])*t)
        qx, qy, qz = [quaternion_from_rotation_vector(np.eye(3)[i]*angles[i]) for i in range(3)]
        q = quaternion_product(launch, quaternion_product(quaternion_product(qx, qy), qz))
        rz, ry = body_to_ecef_matrix(qz), body_to_ecef_matrix(qy)
        omega = rz.T@ry.T[:, 0]*rates[0] + rz.T[:, 1]*rates[1] + np.array([0., 0., rates[2]])
        return q, omega
    q0, w0 = truth(0.)
    k = NavigationKF(p, np.zeros(3), np.r_[np.ones(3), np.full(3, .1), np.full(3, .001)],
        environment=env, q0=quaternion_product(q0, quaternion_from_rotation_vector([.003, -.002, .004])),
        gyro_psd=1e-4**2*dt, accel_psd=.001**2*dt,
        attitude_sigma0=.01, gyro_bias_sigma0=.001,
        estimate_gyro_calibration=True, gyro_sf_sigma0=.003, gyro_mis_sigma0=.003)
    force_world = -env.gravity_ecef(p) - env.frame_acceleration_ecef(p, np.zeros(3))
    f0 = body_to_ecef_matrix(q0).T@force_world
    k.prime_imu(f0 + rng.normal(0, .001, 3), gain@w0 + bias + rng.normal(0, 1e-4, 3))
    trace = []
    raw_errors, corrected_errors = [], []
    for i in range(round(duration/dt) + 1):
        t = i*dt
        q, omega = truth(t)
        rot = body_to_ecef_matrix(q)
        raw = gain@omega + bias + rng.normal(0, 1e-4, 3)
        if i:
            k.predict(rot.T@force_world + rng.normal(0, .001, 3), None, dt, gyro_meas=raw)
        state = type('Truth', (), {
            'quaternion_body_to_ecef_wxyz':q,
            'omega_body_radps':omega,
            'position_body_m':rot.T@p,
        })()
        sample = star.update(t, state)
        if sample is not None and sample.valid:
            k.update_startracker(sample.q_meas_sensor_to_eci_wxyz, sample.t_meas_s,
                spec.q_body_from_sensor, sample.sigma_sensor_rad,
                env.earth_rotation_ecef_radps, sample.latency_s, gate_nis=16.27)
        if i % 5 == 0:
            k.update_gnss(p + rng.normal(0, .01, 3), rng.normal(0, .005, 3),
                None, raw, np.zeros(3), np.full(3, .01), np.full(3, .005))
        if t >= duration - 100:
            raw_errors.append(np.sum((raw - omega)**2))
            corrected_errors.append(np.sum((k.calibrated_omega(raw) - omega)**2))
        if i % 25 == 0:
            trace.append([t, np.rad2deg(angle_error(k.q, q)),
                float(np.linalg.norm(k.bg - bias)), float(np.linalg.norm(k.gyro_sf - scale)),
                float(np.linalg.norm(k.gyro_mis - mis))])
    return {
        'seed':seed, 'duration_s':duration, 'imu_rate_hz':1/dt,
        'motion':'Bounded analytic three-axis sinusoidal attitude; fixed 150 km position',
        'mounting':'Fixture star boresight is Body +x; other star defaults unchanged',
        'final_attitude_error_deg':float(np.rad2deg(angle_error(k.q, q))),
        'true_gyro_bias_radps':bias.tolist(), 'estimated_gyro_bias_radps':k.bg.tolist(),
        'true_gyro_scale':scale.tolist(), 'estimated_gyro_scale':k.gyro_sf.tolist(),
        'true_gyro_misalignment':mis.tolist(), 'estimated_gyro_misalignment':k.gyro_mis.tolist(),
        'bias_error_norm_radps':float(np.linalg.norm(k.bg - bias)),
        'scale_error_norm':float(np.linalg.norm(k.gyro_sf - scale)),
        'misalignment_error_norm_rad':float(np.linalg.norm(k.gyro_mis - mis)),
        'raw_rate_rms_degps':float(np.rad2deg(np.sqrt(np.mean(raw_errors)))),
        'calibrated_rate_rms_degps':float(np.rad2deg(np.sqrt(np.mean(corrected_errors)))),
        'accepted_star_updates':k.startracker_accepted, 'rejected_star_updates':k.startracker_rejected,
        'attitude_tolerance_deg':.03, 'bias_tolerance_radps':2e-5,
        'scale_tolerance':4e-4, 'misalignment_tolerance_rad':8e-4,
        'trace_columns':['time_s','attitude_error_deg','gyro_bias_error_radps','scale_error','misalignment_error_rad'],
        'trace':trace,
    }


def covariance_case(seeds=20, duration=120., dt=.05):
    nees, nis, per_seed = [], [], []
    for seed in range(1, seeds + 1):
        rng = np.random.default_rng(1000 + seed)
        spec = GNSSSpec(np.zeros(3), rate_hz=2., latency_s=0., gm_tau_s=30.)
        sensor = GNSS(spec, seed=2000 + seed)
        bias = rng.normal(0, .03, 3)
        v = np.array([1., -.2, .1])
        k = NavigationKF(rng.normal(0, 2., 3), v + rng.normal(0, .1, 3),
            np.r_[np.full(3, 2.), np.full(3, .1), np.full(3, .03)],
            gravity_ecef=np.zeros(3), accel_psd=.01**2*dt,
            gnss_bias_sigma_m=spec.pos_gm_sigma_m, gnss_bias_tau_s=spec.gm_tau_s)
        k.prime_imu(bias + rng.normal(0, .01, 3), np.zeros(3))
        case_nis = []
        for i in range(round(duration/dt) + 1):
            t = i*dt
            if i:
                k.predict(bias + rng.normal(0, .01, 3), np.eye(3), dt)
            if i % 10 == 0:
                sample = sensor.measure(t, v*t, v)
                k.update_gnss(sample.pos_meas_m, sample.vel_meas_mps, np.eye(3),
                    np.zeros(3), np.zeros(3), spec.pos_sigma_m, spec.vel_sigma_mps)
                case_nis.append(k.last_nis)
        error = k.x[:3] - v*duration
        value = float(error@np.linalg.solve(k.P[:3, :3], error))
        nees.append(value)
        nis.extend(case_nis)
        per_seed.append({'seed':seed,'position_nees':value,'gnss_nis_mean':float(np.mean(case_nis))})
    bounds = np.array([40.481748, 83.297675])/20
    if seeds != 20:
        raise ValueError('The tabulated confidence interval is for 20 trials.')
    return {
        'seeds':seeds,'duration_s':duration,'imu_rate_hz':1/dt,'gnss_rate_hz':2.,
        'gnss_correlation_time_s':30.,'position_anees':float(np.mean(nees)),
        'position_anees_95_bounds':bounds.tolist(),'gnss_nis_mean':float(np.mean(nis)),
        'position_dof':3,'gnss_dof':6,'independent_samples':'20 final-epoch errors, one per independent noise/prior realization',
        'confidence_source':'NIST chi-square CDF 0.025/0.975 for 60 degrees of freedom, divided by 20',
        'per_seed':per_seed,
    }


def integrator_case():
    from dynamics import State, Loads, rk4_step
    from mass_properties import combine_components, fixed_cylinder
    class Model:
        mass_model_name = 'fixed_rotation_fixture'
        def mass_properties(self, configuration, propellant, consumption):
            return combine_components([fixed_cylinder('fixture', 1000., 0., 0., 10., 1.)])
    model = Model()
    omega = np.array([0., 0., .4])
    exact = quaternion_from_rotation_vector(omega*10)
    errors = []
    for step in [.4, .2, .1]:
        state = State([1000., 0., 0.], [0., 0., 0.], [1., 0., 0., 0.], omega, [0., 0.])
        for i in range(round(10/step)):
            state = rk4_step(i*step, state, step, model, None,
                            lambda t,s: Loads(np.zeros(3), np.zeros(3), np.zeros(3), [0., 0.]))
        errors.append(angle_error(state.quaternion_body_to_ecef_wxyz, exact))
    ratios = np.array(errors[:-1])/np.array(errors[1:])
    return {'steps_s':[.4, .2, .1], 'duration_s':10., 'omega_radps':omega.tolist(),
            'attitude_errors_rad':errors, 'error_ratios':ratios.tolist(),
            'observed_orders':np.log2(ratios).tolist()}


def run_all():
    results = {'accelerometer_scale':accel_scale_case(), 'gyro_calibration':attitude_case(),
               'covariance_consistency':covariance_case()}
    a, g, c = results.values()
    results['integrator'] = integrator_case()
    results['checks'] = {
        'accelerometer_scale':a['scale_error_norm'] < a['scale_tolerance'],
        'accelerometer_bias':a['bias_error_norm_mps2'] < a['bias_tolerance_mps2'],
        'attitude':g['final_attitude_error_deg'] < g['attitude_tolerance_deg'],
        'gyro_bias':g['bias_error_norm_radps'] < g['bias_tolerance_radps'],
        'gyro_scale':g['scale_error_norm'] < g['scale_tolerance'],
        'gyro_misalignment':g['misalignment_error_norm_rad'] < g['misalignment_tolerance_rad'],
        'angular_rate':g['calibrated_rate_rms_degps'] < g['raw_rate_rms_degps'],
        'position_anees':c['position_anees_95_bounds'][0] <= c['position_anees'] <= c['position_anees_95_bounds'][1],
        'gnss_nis':5.5 <= c['gnss_nis_mean'] <= 6.5,
        'rk4_order':all(3.8 <= x <= 4.2 for x in results['integrator']['observed_orders']),
    }
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'report/navigation_validation.json')
    args = parser.parse_args()
    result = run_all()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'gyro_calibration'},indent=2))
    print(json.dumps({k:v for k,v in result['gyro_calibration'].items() if k!='trace'},indent=2))
    if not all(result['checks'].values()):
        raise SystemExit(1)
