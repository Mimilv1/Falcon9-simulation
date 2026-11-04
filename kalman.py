"""Navigation Kalman filter V1: IMU propagation + GNSS update + LIDAR update, ECEF.

State x = [p(3), v(3), b_a(3)]
  p, v : IMU station position/velocity, ECEF (origin at the Earth centre)
  b_a  : accelerometer bias, Body
Predict : bias-corrected specific force rotated to ECEF + gravity
GNSS    : antenna position/velocity, lever arm + latency compensation (linear)
LIDAR   : range along boresight to spherical planet, used only when a return is valid
          (planet within reach), lever arm, EKF (Jacobian on position)
V1 limits : attitude supplied externally (truth), no gyro bias, no scale factor / misalignment states,
            no LIDAR bias / scale factor states, constant gravity vector, spherical planet,
            GNSS correlated error not modeled, Earth rotation not modeled (ECEF treated as inertial)
"""

import numpy as np

from sensors import quat_to_matrix

N = 9  # state size
P_SL, V_SL, B_SL = slice(0, 3), slice(3, 6), slice(6, 9)


class NavigationKF:
    def __init__(self, p0, v0, sigma0, gravity_ecef=(0.0, 0.0, -9.81), accel_psd=1e-6,
                 accel_bias_rw=0.0, ba0=None):
        """sigma0: 9 initial 1-sigma [p, v, b_a]; accel_psd: (m/s2)^2*s white noise PSD."""
        self.x = np.r_[np.asarray(p0, float), np.asarray(v0, float),
                       np.zeros(3) if ba0 is None else np.asarray(ba0, float)]
        self.P = np.diag(np.asarray(sigma0, dtype=float) ** 2)
        self.g = np.asarray(gravity_ecef, dtype=float)
        self.accel_psd = float(accel_psd)
        self.accel_bias_rw = float(accel_bias_rw)  # m/s2/sqrt(s)
        self.last_accel = np.zeros(3)  # initial-frame acceleration estimate, for latency compensation
        self.last_nis = np.nan  # normalized innovation squared, last GNSS update
        self.last_lidar_nis = np.nan  # normalized innovation squared, last LIDAR update
        self.log = []

    @classmethod
    def from_imu_spec(cls, imu_spec, p0, v0, sigma0, **kwargs):
        """Process noise from IMUSpec: per-sample sigma -> PSD = sigma^2 / rate."""
        return cls(p0, v0, sigma0,
                   accel_psd=imu_spec.accel_sigma_mps2 ** 2 / imu_spec.rate_hz,
                   accel_bias_rw=imu_spec.accel_bias_rw, **kwargs)

    def predict(self, f_meas, rot, dt):
        """Propagate over dt with IMU specific force (Body) and attitude rot (Body -> ECEF)."""
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

    def _correct(self, H, y, R):
        """Kalman correction, Joseph form. H (m,N), innovation y (m,), noise R (m,m); returns NIS."""
        S = H @ self.P @ H.T + R
        K = np.linalg.solve(S, H @ self.P).T
        self.x = self.x + K @ y
        I_KH = np.eye(N) - K @ H
        self.P = I_KH @ self.P @ I_KH.T + K @ R @ K.T
        self.P = 0.5 * (self.P + self.P.T)
        return float(y @ np.linalg.solve(S, y))

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
        self.last_nis = self._correct(H, y, R)
        return y

    def update_lidar(self, range_meas, rot, arm_body, boresight_body, planet_center, planet_radius, sigma_range):
        """LIDAR range to spherical planet along boresight. Call only with a valid return.

        h(p) = ray-sphere hit distance from sensor = IMU station + lever arm; gradient wrt p.
        Returns innovation, or None if the estimate predicts no surface hit (update skipped).
        """
        rel = self.x[P_SL] + rot @ arm_body - np.asarray(planet_center, float)  # sensor from planet center
        d = rot @ np.asarray(boresight_body, float)
        b = float(rel @ d)
        disc = b * b - (float(rel @ rel) - planet_radius ** 2)
        if b >= 0.0 or disc <= 0.0:
            return None
        root = np.sqrt(disc)
        h = -b - root  # predicted range
        H = np.zeros((1, N))
        H[0, P_SL] = -d - (b * d - rel) / root  # dh/dp
        y = np.array([range_meas - h])
        self.last_lidar_nis = self._correct(H, y, np.array([[float(sigma_range) ** 2]]))
        return y

    def record(self, t):
        self.log.append((t, self.x.copy(), np.sqrt(np.diag(self.P))))

    def to_arrays(self):
        """Log -> (t, x (N,9), sigma (N,9))."""
        return (np.array([e[0] for e in self.log]),
                np.array([e[1] for e in self.log]),
                np.array([e[2] for e in self.log]))


if __name__ == "__main__":
    # demo: fake model, constant thrust + pitch maneuver, IMU + GNSS + LIDAR (reach-limited)
    from types import SimpleNamespace
    from dynamics import State, Loads, rk4_step
    from sensors import IMU, IMUSpec, GNSS, GNSSSpec, Lidar, LidarSpec, truth_inputs
    from frames import (LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ as Q_LAUNCH,
                        launch_position_body_m)

    EARTH_RADIUS_M = 6_371_000.0
    THRUST_N, MASS_KG, G, TORQUE = 20000.0, 1000.0, 9.81, 2.0e4
    G_ECEF = np.array([0.0, 0.0, -G])  # ECEF: z = vertical at the launch site (0,0,R)
    CG = np.array([20.0, 0.0, 0.0])
    IMU_STATION, ANTENNA_STATION, LIDAR_STATION = [10.0, 0.0, 0.0], [30.0, 0.0, 0.0], [3.0, 0.0, 0.0]

    model = SimpleNamespace(
        mass_properties=lambda configuration, propellant, consumption, **kw: SimpleNamespace(
            mass_kg=MASS_KG, cg_station_body_m=CG.copy(), cg_rate_body_mps=np.zeros(3),
            cg_acceleration_body_mps2=np.zeros(3), mass_rate_kgps=0.0, mass_acceleration_kgps2=0.0,
            inertia_body_kgm2=np.diag([1e4, 1e6, 1e6]), inertia_rate_body_kgm2ps=np.zeros((3, 3)),
            arm_from_station=lambda station: np.asarray(station) - CG,
        )
    )

    def provider(t, s):
        rot = quat_to_matrix(s.quaternion_body_to_ecef_wxyz)
        torque = TORQUE if t < 3.0 else (-TORQUE if t < 6.0 else 0.0)  # pitch up then stop rotation
        return Loads(force_body_n=[THRUST_N, 0, 0], moment_body_nm=[0, 0, torque],
                     gravity_body_mps2=rot.T @ G_ECEF, consumption_kgps=[0.0, 0.0])

    def structural_truth(state, props, station):
        """Truth position/velocity of a fixed structural point, ECEF."""
        rot = quat_to_matrix(state.quaternion_body_to_ecef_wxyz)
        arm = np.asarray(station) - props.cg_station_body_m
        return (rot @ (state.position_body_m + arm),
                rot @ (state.velocity_body_mps + np.cross(state.omega_body_radps, arm) - props.cg_rate_body_mps))

    def run_demo(use_gnss, use_lidar, duration_s=60.0, dt=None):
        imu_spec = IMUSpec(IMU_STATION, accel_bias_mps2=[0.05, -0.03, 0.02], accel_scale_factor=1e-3,
                           accel_bias_rw=1e-4, gyro_bias_radps=[0.0, 0.0, 1e-4])  # PLACEHOLDER errors
        gnss_spec = GNSSSpec(ANTENNA_STATION, outages=((20.0, 30.0),))  # no fix 20-30 s
        lidar_spec = LidarSpec(LIDAR_STATION)  # default reach / noise, planet center at the ECEF origin
        dt = 1.0 / imu_spec.rate_hz if dt is None else dt  # one step per IMU sample
        imu, gnss, lidar = IMU(imu_spec, seed=1), GNSS(gnss_spec, seed=2), Lidar(lidar_spec, seed=3)
        gnss_arm = np.asarray(ANTENNA_STATION) - np.asarray(IMU_STATION)
        lidar_arm = np.asarray(LIDAR_STATION) - np.asarray(IMU_STATION)

        state = State(launch_position_body_m(EARTH_RADIUS_M, 20.0), [0, 0, 0], Q_LAUNCH, [0, 0, 0], [0.0, 0.0])  # CG 20 m up: base at the launch site
        t = 0.0
        loads, props = truth_inputs(t, state, model, None, provider)
        p_true0, v_true0 = structural_truth(state, props, IMU_STATION)
        kf = NavigationKF.from_imu_spec(  # initial error: pad survey offset, at rest
            imu_spec, p_true0 + [1.0, -1.0, 0.5], v_true0, [2, 2, 2, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1],
            gravity_ecef=G_ECEF)
        truth_log, nis_log, lidar_nis_log, gnss_log, lidar_valid_t, last = [], [], [], [], [], {"t": None}

        def process(t, state, loads, props):
            rot = quat_to_matrix(state.quaternion_body_to_ecef_wxyz)  # attitude = truth (V1)
            s = imu.update(t, state, loads, props)
            if s is not None:
                if last["t"] is not None:
                    kf.predict(s.accel_meas_mps2, rot, s.t - last["t"])
                last["t"], last["w"] = s.t, s.gyro_meas_radps
            g = gnss.update(t, state, props)  # every step, latency history
            if g is not None and g.valid:  # GNSS-only solution: fix - lever arm, latency extrapolated with own velocity
                z_p = g.pos_meas_m + g.vel_meas_mps * gnss_spec.latency_s - rot @ gnss_arm
                z_v = g.vel_meas_mps - rot @ np.cross(last["w"], gnss_arm)
                gnss_log.append(np.r_[z_p, z_v, *structural_truth(state, props, IMU_STATION), t])
            if use_gnss and g is not None and g.valid:
                kf.update_gnss(g.pos_meas_m, g.vel_meas_mps, rot, last["w"], gnss_arm,
                               g.pos_sigma_m, g.vel_sigma_mps, gnss_spec.latency_s)
                nis_log.append(kf.last_nis)
            lid = lidar.update(t, state, props)
            if lid is not None and lid.valid:  # valid return = planet within reach
                lidar_valid_t.append(t)
                if use_lidar and kf.update_lidar(lid.range_meas_m, rot, lidar_arm, lidar_spec.boresight_body,
                                                 lidar_spec.planet_center_ecef_m, lidar_spec.planet_radius_m,
                                                 lid.sigma_m) is not None:
                    lidar_nis_log.append(kf.last_lidar_nis)
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
        return {"ts": ts, "X": X, "SIG": SIG, "truth": np.array(truth_log), "nis": np.array(nis_log),
                "lidar_nis": np.array(lidar_nis_log), "gnss": np.array(gnss_log),
                "t_low": max(lidar_valid_t)}  # last valid LIDAR return = end of lidar reach

    runs = {"ins": run_demo(False, False), "gnss": run_demo(True, False), "lidar": run_demo(True, True)}
    t_low = runs["lidar"]["t_low"]

    def report(label, run):
        ts, X, SIG, truth = run["ts"], run["X"], run["SIG"], run["truth"]
        err_p, err_v = X[:, 0:3] - truth[:, 0:3], X[:, 3:6] - truth[:, 3:6]
        outage, low = (ts >= 20.0) & (ts < 30.0), ts <= t_low
        rms = lambda e: np.sqrt((e ** 2).sum(axis=1).mean())
        print(f"[{label:<18}] final pos err {np.linalg.norm(err_p[-1]):7.2f} m, vel err {np.linalg.norm(err_v[-1]):6.3f} m/s, "
              f"pos rms {rms(err_p):6.2f} m, rms in outage {rms(err_p[outage]):6.2f} m, "
              f"rms t<={t_low:.1f}s {rms(err_p[low]):5.2f} m (vertical {np.sqrt((err_p[low, 2] ** 2).mean()):5.2f} m), "
              f"in 3-sigma {np.mean(np.abs(err_p) < 3 * SIG[:, 0:3]):.2f}")
        return np.linalg.norm(err_p[-1]), np.sqrt((err_p[low, 2] ** 2).mean())

    final_ins, _ = report("IMU only", runs["ins"])

    # GNSS only: fixes at GNSS rate, none during outage, stats at fix epochs only
    ga = runs["gnss"]["gnss"]
    ge_p, ge_v, low = ga[:, 0:3] - ga[:, 6:9], ga[:, 3:6] - ga[:, 9:12], ga[:, 12] <= t_low
    print(f"[{'GNSS only':<18}] final pos err {np.linalg.norm(ge_p[-1]):7.2f} m, vel err {np.linalg.norm(ge_v[-1]):6.3f} m/s, "
          f"pos rms {np.sqrt((ge_p ** 2).sum(axis=1).mean()):6.2f} m, rms in outage    n/a, "
          f"rms t<={t_low:.1f}s {np.sqrt((ge_p[low] ** 2).sum(axis=1).mean()):5.2f} m "
          f"(vertical {np.sqrt((ge_p[low, 2] ** 2).mean()):5.2f} m), {len(ga)} fixes, no solution in outage")

    final_kf, vert_kf = report("IMU + GNSS", runs["gnss"])
    print(f"  mean NIS {runs['gnss']['nis'].mean():.2f} (6 dof expected ~6), "
          f"est. accel bias {np.round(runs['gnss']['X'][-1, 6:9], 3)}")
    final_lid, vert_lid = report("IMU + GNSS + LIDAR", runs["lidar"])
    ln = runs["lidar"]["lidar_nis"]
    print(f"  lidar updates {len(ln)} (t <= {t_low:.1f} s, reach 1800 m), mean NIS {ln.mean():.2f} (1 dof expected ~1)")

    assert final_kf < 10.0, "filter diverged"
    assert final_kf < final_ins, "GNSS did not help"
    assert vert_lid < vert_kf, "LIDAR did not improve vertical error below its reach"
    print("KALMAN OK")