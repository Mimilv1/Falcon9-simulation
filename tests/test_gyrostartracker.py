import sys
sys.path.insert(0, "stubs")
import numpy as np
from types import SimpleNamespace
from Environement import Environment
from frames import (body_to_ecef_matrix, quaternion_from_rotation_vector, quaternion_product,
                    quaternion_conjugate, quaternion_to_rotation_vector, normalize_quaternion)
from kalman import NavigationKF, misalignment_matrix, _misalignment_jacobian, SF_SL, MIS_SL, BG_SL, TH_SL
from sensors import StarTracker, StarTrackerSpec, attitude_error_rotation_vector

# ---------- 1. Jacobian check: dtheta response to dbg / dsf / dmis -------------
rng = np.random.default_rng(0)
env = Environment()
q0 = normalize_quaternion(rng.normal(size=4))
sf_t, mis_t = np.array([2e-3, -1.5e-3, 1e-3]), rng.normal(0, 1e-3, 6)
bg_t = np.array([5e-4, -3e-4, 2e-4])
omega = np.array([0.03, -0.02, 0.04])
G_t = np.eye(3) + np.diag(sf_t) + misalignment_matrix(mis_t)
z = G_t @ omega + bg_t
dt = 1.0
kf = NavigationKF(np.array([0, 0, 6.5e6]), np.zeros(3), np.r_[np.ones(9), 1e-2*np.ones(3), 1e-3*np.ones(3), 1e-3*np.ones(3), 1e-3*np.ones(6)],
                  environment=env, q0=q0, estimate_gyro_calibration=True)
def att_err_after(dbg, dsf, dmis):
    k = NavigationKF(np.array([0,0,6.5e6]), np.zeros(3), np.ones(24), environment=env, q0=q0,
                     estimate_gyro_calibration=True, bg0=bg_t - dbg, gyro_sf0=sf_t - dsf, gyro_mis0=mis_t - dmis)
    # truth rate over the step = omega (constant); est rate from calibrated gyro
    w_est = k.calibrated_omega(z)
    q_est = quaternion_product(q0, quaternion_from_rotation_vector(w_est * dt))
    q_true = quaternion_product(q0, quaternion_from_rotation_vector(omega * dt))
    return attitude_error_rotation_vector(q_true, q_est)  # dtheta: q_true = q_est (x) exp(dtheta)
eps = 1e-5
F = kf.__class__  # noqa
kf.gyro_sf, kf.gyro_mis, kf.bg = sf_t.copy(), mis_t.copy(), bg_t.copy()
kf._omega_prev = omega
kf.predict(np.zeros(3), None, dt, gyro_meas=z)  # builds F via P propagation only; rebuild explicitly below
g_inv = np.linalg.inv(kf.gyro_matrix)
for name, blk, n in (("bg", g_inv, 3), ("sf", g_inv @ np.diag(omega), 3), ("mis", g_inv @ _misalignment_jacobian(omega), 6)):
    num = np.zeros((3, n))
    for i in range(n):
        d = np.zeros(n); d[i] = eps
        args = [np.zeros(3), np.zeros(3), np.zeros(6)]
        args[{"bg": 0, "sf": 1, "mis": 2}[name]] = d
        num[:, i] = att_err_after(*args) / eps
    print(f"d(theta)/d{name}: max |analytic(-blk*dt) - numeric| = {np.abs(-blk*dt - num).max():.2e}  (scale {np.abs(blk).max():.2e})")

# ---------- 2. Closed-loop: MEKF + star tracker, hovering vehicle, rotating -----
R = env.earth_radius_m
pos = np.array([0.0, 0.0, R + 150e3])
T, imu_dt = 900.0, 0.05
sf_true = np.array([2e-3, -1.5e-3, 1e-3])
mis_true = np.array([1.0e-3, -0.8e-3, 0.6e-3, 1.2e-3, -1.0e-3, 0.5e-3])
G_true = np.eye(3) + np.diag(sf_true) + misalignment_matrix(mis_true)
bg_true = np.array([5e-4, -3e-4, 2e-4])
gyro_sigma, rate_hz = 1e-4, 1 / imu_dt

def run(use_cal, use_st, seed=1):
    rng = np.random.default_rng(seed)
    q_true = normalize_quaternion(np.array([0.9, 0.1, -0.2, 0.3]))
    q_est0 = quaternion_product(q_true, quaternion_from_rotation_vector([3e-3, -2e-3, 4e-3]))
    kf = NavigationKF(pos, np.zeros(3), np.r_[np.full(3, 1.0), np.full(3, 0.1), np.zeros(3) + 1e-6],
                      environment=env, q0=q_est0, gyro_psd=gyro_sigma**2 / rate_hz,
                      accel_psd=1e-8, estimate_gyro_calibration=use_cal,
                      attitude_sigma0=1e-2, gyro_bias_sigma0=1e-3, gyro_sf_sigma0=3e-3, gyro_mis_sigma0=3e-3)
    spec = StarTrackerSpec(boresight_body=[0, 1, 0], max_rate_radps=np.deg2rad(10), min_altitude_m=60e3,
                           earth_exclusion_deg=10.0, alignment_bias_rotvec_rad=np.zeros(3))
    st = StarTracker(spec, env.earth_rotation_ecef_radps, seed=seed + 7)
    t, n = 0.0, int(T / imu_dt)
    errs, nis, reasons = [], [], {}
    for _ in range(n):
        wf = lambda tt: 0.05 * np.array([np.sin(0.05*tt), np.sin(0.037*tt + 1), np.sin(0.023*tt + 2)])
        w = wf(t)
        w_end, w_mid = wf(t + imu_dt), wf(t + imu_dt / 2)
        state = SimpleNamespace(quaternion_body_to_ecef_wxyz=q_true, omega_body_radps=w,
                                position_body_m=pos.copy())
        rot_t = body_to_ecef_matrix(q_true)
        a_model = env.gravity_ecef(pos) + env.frame_acceleration_ecef(pos, np.zeros(3))
        f_meas = rot_t.T @ (-a_model)
        gyro = G_true @ w_end + bg_true + rng.normal(0, gyro_sigma, 3)
        if use_st:
            s = st.update(t, state)
            if s is not None:
                reasons[s.reason] = reasons.get(s.reason, 0) + 1
                if s.valid:
                    kf.update_startracker(s.q_meas_sensor_to_eci_wxyz, s.t_meas_s, spec.q_body_from_sensor,
                                          s.sigma_sensor_rad, env.earth_rotation_ecef_radps, s.latency_s,
                                          gate_nis=16.27)
                    nis.append(kf.last_st_nis)
        kf.predict(f_meas, None, imu_dt, gyro_meas=gyro)
        q_true = normalize_quaternion(quaternion_product(q_true, quaternion_from_rotation_vector(w_mid * imu_dt)))
        t += imu_dt
        errs.append(np.rad2deg(np.linalg.norm(attitude_error_rotation_vector(kf.q, q_true))))
    return kf, np.array(errs), nis, reasons, st

for use_cal, use_st in ((False, False), (False, True), (True, True)):
    kf, errs, nis, reasons, st = run(use_cal, use_st)
    msg = f"cal={use_cal!s:5} st={use_st!s:5}: attitude err final={errs[-1]*3600:8.1f}\"  rms(last 300s)={np.sqrt(np.mean(errs[-6000:]**2))*3600:8.1f}\""
    if use_st:
        msg += f"  ST NIS mean={np.mean(nis):.2f} (3 dof) rejected={kf.startracker_rejected} n={len(nis)}"
    print(msg)
    if use_cal:
        print("  sf   true", sf_true, "\n  sf   est ", np.round(kf.gyro_sf, 5), "\n  sf   1sig", np.round(kf.gyro_sf_sigma, 5))
        print("  mis  true", mis_true, "\n  mis  est ", np.round(kf.gyro_mis, 5), "\n  mis  1sig", np.round(kf.gyro_mis_sigma_rad, 5))
        print("  bg   true", bg_true, " est", np.round(kf.bg, 6))
print("ST reasons:", reasons)