"""Navigation Kalman filter V1: IMU propagation + GNSS update, initial frame.

State x = [p(3), v(3), b_a(3)]
  p, v : IMU station position/velocity, initial frame
  b_a  : accelerometer bias, Body
Predict : bias-corrected specific force rotated to initial frame + gravity
Update  : GNSS antenna position/velocity, lever arm + latency compensation
V1 limits : attitude supplied externally (truth), no gyro bias, no scale factor / misalignment states,
            constant gravity vector, GNSS correlated error not modeled
"""

import numpy as np

from sensors import quat_to_matrix

N = 9  # state size
P_SL, V_SL, B_SL = slice(0, 3), slice(3, 6), slice(6, 9)


class NavigationKF:
    def __init__(self, p0, v0, sigma0, gravity_init=(-9.81, 0.0, 0.0), accel_psd=1e-6,
                 accel_bias_rw=0.0, ba0=None):
        """sigma0: 9 initial 1-sigma [p, v, b_a]; accel_psd: (m/s2)^2*s white noise PSD."""
        self.x = np.r_[np.asarray(p0, float), np.asarray(v0, float),
                       np.zeros(3) if ba0 is None else np.asarray(ba0, float)]
        self.P = np.diag(np.asarray(sigma0, dtype=float) ** 2)
        self.g = np.asarray(gravity_init, dtype=float)
        self.accel_psd = float(accel_psd)
        self.accel_bias_rw = float(accel_bias_rw)  # m/s2/sqrt(s)
        self.last_accel = np.zeros(3)  # initial-frame acceleration estimate, for latency compensation
        self.last_nis = np.nan  # normalized innovation squared, last update
        self.log = []

    @classmethod
    def from_imu_spec(cls, imu_spec, p0, v0, sigma0, **kwargs):
        """Process noise from IMUSpec: per-sample sigma -> PSD = sigma^2 / rate."""
        return cls(p0, v0, sigma0,
                   accel_psd=imu_spec.accel_sigma_mps2 ** 2 / imu_spec.rate_hz,
                   accel_bias_rw=imu_spec.accel_bias_rw, **kwargs)

    def predict(self, f_meas, rot, dt):
        """Propagate over dt with IMU specific force (Body) and attitude rot (Body -> initial)."""
        p, v, ba = self.x[P_SL].copy(), self.x[V_SL].copy(), self.x[B_SL].copy()
        a = rot @ (f_meas - ba) + self.g
        self.x[P_SL] = p + v * dt + 0.5 * a * dt ** 2
        self.x[V_SL] = v + a * dt
        self.last_accel = a

        eye = np.eye(3)
        F = np.eye(N)
        F[P_SL, V_SL] = eye * dt
        F[P_SL, B_SL] = -0.5 * rot * dt ** 2
        F[V_SL, B_SL] = -rot * dt
        S = self.accel_psd
        Q = np.zeros((N, N))
        Q[P_SL, P_SL] = S * dt ** 3 / 3 * eye
        Q[P_SL, V_SL] = Q[V_SL, P_SL] = S * dt ** 2 / 2 * eye
        Q[V_SL, V_SL] = S * dt * eye
        Q[B_SL, B_SL] = self.accel_bias_rw ** 2 * dt * eye
        self.P = F @ self.P @ F.T + Q

    def update_gnss(self, pos_meas, vel_meas, rot, omega, arm_body, sigma_pos, sigma_vel, latency_s=0.0):
        """GNSS fix: antenna = IMU station + lever arm (fixed in Body), measurement delayed by latency_s."""
        a, tau = self.last_accel, latency_s
        z_p = pos_meas + vel_meas * tau + 0.5 * a * tau ** 2  # extrapolate to now
        z_v = vel_meas + a * tau
        z = np.r_[z_p - rot @ arm_body, z_v - rot @ np.cross(omega, arm_body)]  # remove lever arm

        H = np.zeros((6, N))
        H[0:3, P_SL] = np.eye(3)
        H[3:6, V_SL] = np.eye(3)
        R = np.diag(np.r_[np.asarray(sigma_pos, float) ** 2, np.asarray(sigma_vel, float) ** 2])

        y = z - H @ self.x  # innovation
        S = H @ self.P @ H.T + R
        K = np.linalg.solve(S, H @ self.P).T
        self.x = self.x + K @ y
        I_KH = np.eye(N) - K @ H
        self.P = I_KH @ self.P @ I_KH.T + K @ R @ K.T  # Joseph form
        self.P = 0.5 * (self.P + self.P.T)
        self.last_nis = float(y @ np.linalg.solve(S, y))
        return y

    def record(self, t):
        self.log.append((t, self.x.copy(), np.sqrt(np.diag(self.P))))

    def to_arrays(self):
        """Log -> (t, x (N,9), sigma (N,9))."""
        return (np.array([e[0] for e in self.log]),
                np.array([e[1] for e in self.log]),
                np.array([e[2] for e in self.log]))


if __name__ == "__main__":
    # demo: fake model, constant thrust + pitch maneuver, IMU 100 Hz + GNSS 1 Hz
    from types import SimpleNamespace
    from dynamics import State, Loads, rk4_step
    from sensors import IMU, IMUSpec, GNSS, GNSSSpec, truth_inputs

    THRUST_N, MASS_KG, G, TORQUE = 20000.0, 1000.0, 9.81, 2.0e4
    G_INIT = np.array([-G, 0.0, 0.0])  # x = vertical at launch
    CG = np.array([20.0, 0.0, 0.0])
    IMU_STATION, ANTENNA_STATION = [10.0, 0.0, 0.0], [30.0, 0.0, 0.0]

    model = SimpleNamespace(
        mass_properties=lambda configuration, propellant, consumption, **kw: SimpleNamespace(
            mass_kg=MASS_KG, cg_station_body_m=CG.copy(), cg_rate_body_mps=np.zeros(3),
            cg_acceleration_body_mps2=np.zeros(3), mass_rate_kgps=0.0, mass_acceleration_kgps2=0.0,
            inertia_body_kgm2=np.diag([1e4, 1e6, 1e6]), inertia_rate_body_kgm2ps=np.zeros((3, 3)),
            arm_from_station=lambda station: np.asarray(station) - CG,
        )
    )

    def provider(t, s):
        rot = quat_to_matrix(s.quaternion_body_to_initial_body_wxyz)
        torque = TORQUE if t < 3.0 else (-TORQUE if t < 6.0 else 0.0)  # pitch up then stop rotation
        return Loads(force_body_n=[THRUST_N, 0, 0], moment_body_nm=[0, 0, torque],
                     gravity_body_mps2=rot.T @ G_INIT, consumption_kgps=[0.0, 0.0])

    def structural_truth(state, props, station):
        """Truth position/velocity of a fixed structural point, initial frame."""
        rot = quat_to_matrix(state.quaternion_body_to_initial_body_wxyz)
        arm = np.asarray(station) - props.cg_station_body_m
        return (rot @ (state.position_body_m + arm),
                rot @ (state.velocity_body_mps + np.cross(state.omega_body_radps, arm) - props.cg_rate_body_mps))

    def run_demo(use_gnss, duration_s=60.0, dt=None):
        imu_spec = IMUSpec(IMU_STATION, accel_bias_mps2=[0.05, -0.03, 0.02], accel_scale_factor=1e-3,
                           accel_bias_rw=1e-4, gyro_bias_radps=[0.0, 0.0, 1e-4])  # PLACEHOLDER errors
        gnss_spec = GNSSSpec(ANTENNA_STATION, outages=((20.0, 30.0),))  # no fix 20-30 s
        dt = 1.0 / imu_spec.rate_hz if dt is None else dt  # one step per IMU sample
        imu, gnss = IMU(imu_spec, seed=1), GNSS(gnss_spec, seed=2)
        arm_body = np.asarray(ANTENNA_STATION) - np.asarray(IMU_STATION)

        state = State([0, 0, 0], [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0.0, 0.0])
        t = 0.0
        loads, props = truth_inputs(t, state, model, None, provider)
        p_true0, v_true0 = structural_truth(state, props, IMU_STATION)
        kf = NavigationKF.from_imu_spec(  # initial error: pad survey offset, at rest
            imu_spec, p_true0 + [1.0, -1.0, 0.5], v_true0, [2, 2, 2, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1],
            gravity_init=G_INIT)
        truth_log, nis_log, gnss_log, last = [], [], [], {"t": None}

        def process(t, state, loads, props):
            rot = quat_to_matrix(state.quaternion_body_to_initial_body_wxyz)  # attitude = truth (V1)
            s = imu.update(t, state, loads, props)
            if s is not None:
                if last["t"] is not None:
                    kf.predict(s.accel_meas_mps2, rot, s.t - last["t"])
                last["t"], last["w"] = s.t, s.gyro_meas_radps
            g = gnss.update(t, state, props)  # every step, latency history
            if g is not None and g.valid:  # GNSS-only solution: fix - lever arm, latency extrapolated with own velocity
                z_p = g.pos_meas_m + g.vel_meas_mps * gnss_spec.latency_s - rot @ arm_body
                z_v = g.vel_meas_mps - rot @ np.cross(last["w"], arm_body)
                gnss_log.append(np.r_[z_p, z_v, *structural_truth(state, props, IMU_STATION)])
            if use_gnss and g is not None and g.valid:
                kf.update_gnss(g.pos_meas_m, g.vel_meas_mps, rot, last["w"], arm_body,
                               g.pos_sigma_m, g.vel_sigma_mps, gnss_spec.latency_s)
                nis_log.append(kf.last_nis)
            if s is not None:
                kf.record(t)
                truth_log.append(np.r_[structural_truth(state, props, IMU_STATION)])

        process(t, state, loads, props)
        while t < duration_s - 1e-9:
            state = rk4_step(t, state, dt, model, None, provider)  # accepted step
            t += dt
            loads, props = truth_inputs(t, state, model, None, provider)
            process(t, state, loads, props)

        ts, X, SIG = kf.to_arrays()
        truth = np.array(truth_log)
        return ts, X, SIG, truth, np.array(nis_log), state, np.array(gnss_log)

    results = {flag: run_demo(flag) for flag in (False, True)}  # IMU only / IMU + GNSS, same seeds

    def report(label, ts, X, SIG, truth):
        err_p, err_v = X[:, 0:3] - truth[:, 0:3], X[:, 3:6] - truth[:, 3:6]
        outage = (ts >= 20.0) & (ts < 30.0)
        print(f"[{label:<10}] final pos err {np.linalg.norm(err_p[-1]):7.2f} m, vel err {np.linalg.norm(err_v[-1]):6.3f} m/s, "
              f"pos rms {np.sqrt((err_p ** 2).sum(axis=1).mean()):6.2f} m, "
              f"rms in outage {np.sqrt((err_p[outage] ** 2).sum(axis=1).mean()):6.2f} m, "
              f"in 3-sigma {np.mean(np.abs(err_p) < 3 * SIG[:, 0:3]):.2f}")
        return np.linalg.norm(err_p[-1])

    ts, X, SIG, truth, nis, final, gnss_arr = results[False]
    final_ins = report("IMU only", ts, X, SIG, truth)

    # GNSS only: fixes at GNSS rate, none during outage, stats at fix epochs only
    gnss_arr = results[True][6]
    ge_p, ge_v = gnss_arr[:, 0:3] - gnss_arr[:, 6:9], gnss_arr[:, 3:6] - gnss_arr[:, 9:12]
    print(f"[{'GNSS only':<10}] final pos err {np.linalg.norm(ge_p[-1]):7.2f} m, vel err {np.linalg.norm(ge_v[-1]):6.3f} m/s, "
          f"pos rms {np.sqrt((ge_p ** 2).sum(axis=1).mean()):6.2f} m, rms in outage    n/a, "
          f"{len(gnss_arr)} fixes, no solution in outage")

    ts, X, SIG, truth, nis, final, gnss_arr = results[True]
    final_kf = report("IMU + GNSS", ts, X, SIG, truth)
    print(f"  mean NIS {nis.mean():.2f} (6 dof expected ~6), est. accel bias {np.round(X[-1, 6:9], 3)}")
    assert final_kf < 10.0, "filter diverged"
    assert final_kf < final_ins, "GNSS did not help"
    print("KALMAN OK")