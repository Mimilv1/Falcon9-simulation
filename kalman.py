"""Error-state navigation with attitude and sensor-error estimation."""

from bisect import bisect_left, bisect_right
import numpy as np

from frames import (
    body_to_ecef_matrix,
    normalize_quaternion,
    quaternion_from_rotation_vector,
    quaternion_product,
    skew,
    quaternion_conjugate,
    quaternion_to_rotation_vector,
    earth_rotation_quaternion,
)

N_NAV = 9
N_ATT = 15
N_CAL = 24
N = N_NAV  # backward-compatible alias
P_SL, V_SL, B_SL = slice(0, 3), slice(3, 6), slice(6, 9)
TH_SL, BG_SL = slice(9, 12), slice(12, 15)
SF_SL, MIS_SL = slice(15, 18), slice(18, 24)
# Off-diagonal gyro misalignment terms (row i, column j of G), in state order
MIS_PAIRS = ((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1))


def misalignment_matrix(terms):
    """3x3 zero-diagonal matrix from the 6 off-diagonal terms (MIS_PAIRS order)."""
    m = np.zeros((3, 3))
    for value, (i, j) in zip(np.asarray(terms, dtype=float), MIS_PAIRS):
        m[i, j] = value
    return m


def misalignment_terms(matrix):
    """Inverse of misalignment_matrix (diagonal ignored)."""
    m = np.asarray(matrix, dtype=float)
    return np.array([m[i, j] for i, j in MIS_PAIRS])


def _misalignment_jacobian(omega):
    """d(M(mis) @ omega)/d(mis), shape (3, 6)."""
    jac = np.zeros((3, 6))
    for k, (i, j) in enumerate(MIS_PAIRS):
        jac[i, k] = omega[j]
    return jac


def _right_jacobian(vector):
    angle = np.linalg.norm(vector)
    cross = skew(vector)
    if angle < 1e-5:
        return np.eye(3) - 0.5 * cross + cross @ cross / 6
    return (
        np.eye(3)
        - (1 - np.cos(angle)) / angle**2 * cross
        + (angle - np.sin(angle)) / angle**3 * cross @ cross
    )


class NavigationKF:
    def __init__(
        self,
        p0,
        v0,
        sigma0,
        gravity_ecef=(0.0, 0.0, -9.81),
        accel_psd=1e-6,
        accel_bias_rw=0.0,
        ba0=None,
        environment=None,
        q0=None,
        bg0=None,
        gyro_psd=1e-8,
        gyro_bias_rw=0.0,
        attitude_sigma0=1e-2,
        gyro_bias_sigma0=1e-3,
        estimate_gyro_calibration=False,
        gyro_sf0=None,
        gyro_mis0=None,
        gyro_sf_sigma0=1e-3,
        gyro_mis_sigma0=1e-3,
        gyro_sf_rw=0.0,
        gyro_mis_rw=0.0,
        estimate_accel_scale=False,
        accel_sf_sigma0=1e-3,
        accel_sf_rw=0.0,
        gnss_bias_sigma_m=None,
        gnss_bias_tau_s=300.0,
    ):
        self.x = np.r_[
            np.asarray(p0, float),
            np.asarray(v0, float),
            np.zeros(3) if ba0 is None else np.asarray(ba0, float),
        ]
        self.attitude_enabled = q0 is not None
        self.calibration_enabled = bool(estimate_gyro_calibration)
        if self.calibration_enabled and not self.attitude_enabled:
            raise ValueError("Gyro calibration needs the attitude block (give q0).")
        sigma0 = np.asarray(sigma0, dtype=float)
        # Gyro calibration estimate (always defined; stays 0 if not estimated)
        self.gyro_sf = (
            np.zeros(3)
            if gyro_sf0 is None
            else np.broadcast_to(np.asarray(gyro_sf0, float), (3,)).copy()
        )
        if gyro_mis0 is None:
            self.gyro_mis = np.zeros(6)
        else:
            mis0 = np.asarray(gyro_mis0, dtype=float)
            self.gyro_mis = (
                misalignment_terms(mis0) if mis0.shape == (3, 3) else mis0.copy()
            )
            if self.gyro_mis.shape != (6,):
                raise ValueError("gyro_mis0 must be a 3x3 matrix or 6 terms.")
        if self.attitude_enabled:
            self.n = N_CAL if self.calibration_enabled else N_ATT
            self.q = normalize_quaternion(q0)
            self.bg = np.zeros(3) if bg0 is None else np.asarray(bg0, dtype=float).copy()
            if sigma0.shape == (N_NAV,):
                sigma0 = np.r_[
                    sigma0,
                    np.broadcast_to(np.asarray(attitude_sigma0, float), (3,)),
                    np.broadcast_to(np.asarray(gyro_bias_sigma0, float), (3,)),
                ]
            if self.calibration_enabled and sigma0.shape == (N_ATT,):
                sigma0 = np.r_[
                    sigma0,
                    np.broadcast_to(np.asarray(gyro_sf_sigma0, float), (3,)),
                    np.broadcast_to(np.asarray(gyro_mis_sigma0, float), (6,)),
                ]
        else:
            self.n = N_NAV
            self.q = None
            self.bg = np.zeros(3)
        extend_sigma = sigma0.shape == (self.n,)
        self.accel_sf_enabled = bool(estimate_accel_scale)
        self.gnss_bias_enabled = gnss_bias_sigma_m is not None
        self.accel_sf = np.zeros(3)
        self.gnss_bias = np.zeros(3)
        self.accel_sf_slice = slice(self.n, self.n + 3)
        if self.accel_sf_enabled:
            self.n += 3
            if extend_sigma:
                sigma0 = np.r_[sigma0, np.full(3, accel_sf_sigma0)]
        self.gnss_bias_slice = slice(self.n, self.n + 3)
        self.gnss_bias_sigma = (
            np.zeros(3)
            if gnss_bias_sigma_m is None
            else np.broadcast_to(np.asarray(gnss_bias_sigma_m, float), (3,)).copy()
        )
        self.gnss_bias_tau_s = float(gnss_bias_tau_s)
        self.accel_sf_rw = float(accel_sf_rw)
        if self.gnss_bias_tau_s <= 0 or np.any(self.gnss_bias_sigma < 0):
            raise ValueError("Require positive correlation time and nonnegative sigma.")
        if self.gnss_bias_enabled:
            self.n += 3
            if extend_sigma:
                sigma0 = np.r_[sigma0, self.gnss_bias_sigma]
        if sigma0.shape != (self.n,):
            raise ValueError(f"sigma0 must have {self.n} entries.")
        self.P = np.diag(sigma0**2)
        self.g = np.asarray(gravity_ecef, dtype=float)
        self.environment = environment
        self.accel_psd = float(accel_psd)
        self.accel_bias_rw = float(accel_bias_rw)
        self.gyro_psd = float(gyro_psd)
        self.gyro_bias_rw = float(gyro_bias_rw)
        self.gyro_sf_rw = float(gyro_sf_rw)
        self.gyro_mis_rw = float(gyro_mis_rw)
        self._omega_prev = None
        self._gyro_prev = None
        self._force_prev = None
        self._time_s = 0.0
        self._gyro_times = []
        self._gyro_values = []
        self.last_st_nis = np.nan
        self.startracker_rejected = 0
        self.startracker_accepted = 0
        self.startracker_diagnostics = []
        self.last_accel = np.zeros(3)
        self.last_nis = np.nan
        self.last_lidar_nis = np.nan
        self.log = []

    def prime_imu(self, force, gyro):
        self._force_prev = np.asarray(force, float).copy()
        self._gyro_prev = np.asarray(gyro, float).copy()
        self._omega_prev = self.calibrated_omega(gyro)
        self._gyro_times = [self._time_s]
        self._gyro_values = [self._gyro_prev.copy()]

    def _gyro_increment(self, start, end):
        if end <= start:
            return np.array([1., 0, 0, 0]), {}, np.zeros((3, 3))
        if not self._gyro_times:
            return (
                quaternion_from_rotation_vector(self._omega_prev * (end - start)),
                {}, self.gyro_psd * (end - start) * np.eye(3),
            )
        first = max(0, bisect_right(self._gyro_times, start) - 1)
        last = min(len(self._gyro_times), bisect_left(self._gyro_times, end) + 1)
        times = np.asarray(self._gyro_times[first:last])
        values = np.asarray(self._gyro_values[first:last])
        knots = np.r_[start, times[(times > start) & (times < end)], end]
        raw = np.column_stack([np.interp(knots, times, values[:, i]) for i in range(3)])
        delta = np.array([1., 0, 0, 0])
        maps = {BG_SL.start: np.zeros((3, 3))}
        if self.calibration_enabled:
            maps[SF_SL.start] = np.zeros((3, 3))
            maps[MIS_SL.start] = np.zeros((3, 6))
        noise = np.zeros((3, 3))
        inverse = np.linalg.inv(self.gyro_matrix)
        for i, dt in enumerate(np.diff(knots)):
            omega = self.calibrated_omega((raw[i] + raw[i + 1]) / 2)
            phi = omega * dt
            dq = quaternion_from_rotation_vector(phi)
            rotation = body_to_ecef_matrix(dq)
            jac = _right_jacobian(phi) @ inverse
            derivatives = {BG_SL.start: np.eye(3)}
            if self.calibration_enabled:
                derivatives[SF_SL.start] = np.diag(omega)
                derivatives[MIS_SL.start] = _misalignment_jacobian(omega)
            for index, derivative in derivatives.items():
                maps[index] = rotation.T @ maps[index] + jac @ derivative * dt
            noise = rotation.T @ noise @ rotation + self.gyro_psd * jac @ jac.T * dt
            delta = normalize_quaternion(quaternion_product(delta, dq))
        return delta, maps, noise

    @classmethod
    def from_imu_spec(cls, imu_spec, p0, v0, sigma0, **kwargs):
        kwargs.setdefault("gyro_psd", imu_spec.gyro_sigma_radps**2 / imu_spec.rate_hz)
        kwargs.setdefault("gyro_bias_rw", imu_spec.gyro_bias_rw)
        return cls(
            p0,
            v0,
            sigma0,
            accel_psd=imu_spec.accel_sigma_mps2**2 / imu_spec.rate_hz,
            accel_bias_rw=imu_spec.accel_bias_rw,
            **kwargs
        )

    # ------------------------------------------------------------------ helpers
    @property
    def rotation(self):
        """Estimated body -> ECEF rotation matrix (None if attitude is not estimated)."""
        return body_to_ecef_matrix(self.q) if self.attitude_enabled else None

    @property
    def attitude_sigma_rad(self):
        return np.sqrt(np.diag(self.P)[TH_SL]) if self.attitude_enabled else None

    @property
    def gyro_bias_sigma_radps(self):
        return np.sqrt(np.diag(self.P)[BG_SL]) if self.attitude_enabled else None

    @property
    def gyro_matrix(self):
        """Estimated gyro matrix G = I + diag(sf) + M (identity if not estimated)."""
        if not self.calibration_enabled:
            return np.eye(3)
        return np.eye(3) + np.diag(self.gyro_sf) + misalignment_matrix(self.gyro_mis)

    @property
    def gyro_sf_sigma(self):
        return np.sqrt(np.diag(self.P)[SF_SL]) if self.calibration_enabled else None

    @property
    def gyro_mis_sigma_rad(self):
        return np.sqrt(np.diag(self.P)[MIS_SL]) if self.calibration_enabled else None

    def calibrated_omega(self, gyro_meas):
        """Body rate estimate G^-1 (w_meas - bg)."""
        raw = np.asarray(gyro_meas, dtype=float) - self.bg
        if not self.calibration_enabled:
            return raw
        return np.linalg.solve(self.gyro_matrix, raw)

    def _model_acceleration(self, p, v):
        env = self.environment
        if env is None:
            return self.g, np.zeros((3, 3)), np.zeros((3, 3))
        d_cent, d_cor = env.frame_acceleration_gradients_ecef()
        return (
            env.gravity_ecef(p) + env.frame_acceleration_ecef(p, v),
            env.gravity_gradient_ecef(p) + d_cent,
            d_cor,
        )

    # ------------------------------------------------------------------ predict
    def predict(self, f_meas, rot, dt, gyro_meas=None):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Prediction interval must be positive and finite.")
        p, v, ba = self.x[P_SL].copy(), self.x[V_SL].copy(), self.x[B_SL].copy()
        model_a, d_dp, d_dv = self._model_acceleration(p, v)
        force = np.asarray(f_meas, float)
        force_mean = force if self._force_prev is None else (force + self._force_prev) / 2
        half_rotation = np.eye(3)
        if self.attitude_enabled:
            if gyro_meas is None:
                raise ValueError("gyro_meas is required when attitude is estimated.")
            gyro = np.asarray(gyro_meas, float)
            previous = gyro if self._gyro_prev is None else self._gyro_prev
            omega_avg = self.calibrated_omega((gyro + previous) / 2)
            phi = omega_avg * dt
            half_rotation = body_to_ecef_matrix(quaternion_from_rotation_vector(phi / 2))
            rot = self.rotation @ half_rotation
        gain = 1 + self.accel_sf
        if np.any(np.abs(gain) < 1e-8):
            raise ValueError("Singular accelerometer calibration.")
        f_corrected = (force_mean - ba) / gain
        a = rot @ f_corrected + model_a
        self.x[P_SL] = p + v * dt + 0.5 * a * dt**2
        self.x[V_SL] = v + a * dt
        self.last_accel = a
        eye = np.eye(3)
        F = np.eye(self.n)
        F[P_SL, P_SL] = eye + 0.5 * d_dp * dt**2
        F[P_SL, V_SL] = eye * dt + 0.5 * d_dv * dt**2
        F[V_SL, P_SL] = d_dp * dt
        F[V_SL, V_SL] = eye + d_dv * dt
        accel_bias_map = rot @ np.diag(1 / gain)
        F[P_SL, B_SL] = -0.5 * accel_bias_map * dt**2
        F[V_SL, B_SL] = -accel_bias_map * dt
        Q = np.zeros((self.n, self.n))
        noise = self.accel_psd * accel_bias_map @ accel_bias_map.T
        Q[P_SL, P_SL] = noise * dt**3 / 3
        Q[P_SL, V_SL] = Q[V_SL, P_SL] = noise * dt**2 / 2
        Q[V_SL, V_SL] = noise * dt
        Q[B_SL, B_SL] = self.accel_bias_rw**2 * dt * eye
        if self.accel_sf_enabled:
            sf_map = rot @ np.diag(f_corrected / gain)
            F[V_SL, self.accel_sf_slice] = -sf_map * dt
            F[P_SL, self.accel_sf_slice] = -0.5 * sf_map * dt**2
            Q[self.accel_sf_slice, self.accel_sf_slice] = self.accel_sf_rw**2 * dt * eye
        if self.gnss_bias_enabled:
            decay = np.exp(-dt / self.gnss_bias_tau_s)
            self.gnss_bias *= decay
            F[self.gnss_bias_slice, self.gnss_bias_slice] = decay * eye
            Q[self.gnss_bias_slice, self.gnss_bias_slice] = np.diag(
                self.gnss_bias_sigma**2 * (1 - decay**2)
            )
        if self.attitude_enabled:
            r_fx = rot @ skew(f_corrected)
            F[V_SL, TH_SL] = -r_fx @ half_rotation.T * dt
            F[P_SL, TH_SL] = -0.5 * r_fx @ half_rotation.T * dt**2
            dq = quaternion_from_rotation_vector(phi)
            F[TH_SL, TH_SL] = body_to_ecef_matrix(dq).T
            g_inv = np.linalg.inv(self.gyro_matrix)
            full_map = _right_jacobian(phi) @ g_inv
            half_map = _right_jacobian(phi / 2) @ g_inv
            blocks = [(BG_SL, eye)]
            if self.calibration_enabled:
                blocks += [(SF_SL, np.diag(omega_avg)),
                           (MIS_SL, _misalignment_jacobian(omega_avg))]
            for block, rate_map in blocks:
                F[TH_SL, block] = -full_map @ rate_map * dt
                F[V_SL, block] = 0.5 * r_fx @ half_map @ rate_map * dt**2
                F[P_SL, block] = 0.25 * r_fx @ half_map @ rate_map * dt**3
            Q[TH_SL, TH_SL] = self.gyro_psd * full_map @ full_map.T * dt
            Q[BG_SL, BG_SL] = self.gyro_bias_rw**2 * dt * eye
            if self.calibration_enabled:
                Q[SF_SL, SF_SL] = self.gyro_sf_rw**2 * dt * eye
                Q[MIS_SL, MIS_SL] = self.gyro_mis_rw**2 * dt * np.eye(6)
            self.q = normalize_quaternion(quaternion_product(self.q, dq))
            self._gyro_prev = gyro.copy()
            self._omega_prev = self.calibrated_omega(gyro)
            if not self._gyro_times:
                self._gyro_times.append(self._time_s)
                self._gyro_values.append(previous.copy())
            self._gyro_times.append(self._time_s + dt)
            self._gyro_values.append(gyro.copy())
        self._force_prev = force.copy()
        self._time_s += dt
        self.last_transition = F
        self.last_process_covariance = Q
        self.P = F @ self.P @ F.T + Q
        self.P = (self.P + self.P.T) / 2

    # ------------------------------------------------------------------ correct
    def _correct(self, H, y, R):
        S = H @ self.P @ H.T + R
        self.last_measurement_jacobian = H.copy()
        self.last_innovation = y.copy()
        self.last_innovation_covariance = S.copy()
        K = np.linalg.solve(S, H @ self.P).T
        dx = K @ y
        self.x = self.x + dx[:N_NAV]
        if self.accel_sf_enabled:
            self.accel_sf += dx[self.accel_sf_slice]
        if self.gnss_bias_enabled:
            self.gnss_bias += dx[self.gnss_bias_slice]
        I_KH = np.eye(self.n) - K @ H
        P = I_KH @ self.P @ I_KH.T + K @ R @ K.T
        if self.attitude_enabled:
            dtheta = dx[TH_SL]
            self.q = normalize_quaternion(
                quaternion_product(self.q, quaternion_from_rotation_vector(dtheta))
            )
            self.bg = self.bg + dx[BG_SL]
            if self.calibration_enabled:
                self.gyro_sf = self.gyro_sf + dx[SF_SL]
                self.gyro_mis = self.gyro_mis + dx[MIS_SL]
            # Covariance reset: express P in the updated attitude error frame
            G = np.eye(self.n)
            G[TH_SL, TH_SL] = np.eye(3) - skew(0.5 * dtheta)
            P = G @ P @ G.T
            if self._gyro_prev is not None:
                self._omega_prev = self.calibrated_omega(self._gyro_prev)
        self.P = 0.5 * (P + P.T)
        return float(y @ np.linalg.solve(S, y))

    def update_gnss(
        self,
        pos_meas,
        vel_meas,
        rot,
        omega,
        arm_body,
        sigma_pos,
        sigma_vel,
        latency_s=0.0,
    ):
        if self.attitude_enabled:
            rot = body_to_ecef_matrix(self.q)
            omega = self.calibrated_omega(omega)
        a, tau = self.last_accel, latency_s
        z_p = pos_meas + vel_meas * tau + 0.5 * a * tau**2
        z_v = vel_meas + a * tau
        lever_velocity = np.cross(omega, arm_body)
        z = np.r_[z_p - rot @ arm_body, z_v - rot @ lever_velocity]

        H = np.zeros((6, self.n))
        H[0:3, P_SL] = np.eye(3)
        H[3:6, V_SL] = np.eye(3)
        if self.attitude_enabled:
            # Lever-arm terms depend on attitude and gyro bias
            H[0:3, TH_SL] = -rot @ skew(arm_body)
            H[3:6, TH_SL] = -rot @ skew(lever_velocity)
            lever = rot @ skew(arm_body) @ np.linalg.inv(self.gyro_matrix)
            H[3:6, BG_SL] = lever
            if self.calibration_enabled:
                H[3:6, SF_SL] = lever @ np.diag(omega)
                H[3:6, MIS_SL] = lever @ _misalignment_jacobian(omega)
        velocity_cov = np.diag(np.asarray(sigma_vel, float) ** 2)
        R = np.zeros((6, 6))
        R[:3, :3] = np.diag(np.asarray(sigma_pos, float) ** 2) + tau**2 * velocity_cov
        R[:3, 3:] = R[3:, :3] = tau * velocity_cov
        R[3:, 3:] = velocity_cov
        prediction = np.r_[self.x[P_SL], self.x[V_SL]]
        if self.gnss_bias_enabled:
            prediction[:3] += self.gnss_bias
            H[:3, self.gnss_bias_slice] = np.eye(3)
        y = z - prediction
        self.last_nis = self._correct(H, y, R)
        return y

    def update_lidar(
        self,
        range_meas,
        rot,
        arm_body,
        boresight_body,
        planet_center,
        planet_radius,
        sigma_range,
    ):
        if self.attitude_enabled:
            rot = body_to_ecef_matrix(self.q)
        boresight = np.asarray(boresight_body, float)
        rel = self.x[P_SL] + rot @ arm_body - np.asarray(planet_center, float)
        d = rot @ boresight
        b = float(rel @ d)
        disc = b * b - (float(rel @ rel) - planet_radius**2)
        if b >= 0.0 or disc <= 0.0:
            return None
        root = np.sqrt(disc)
        h = -b - root
        H = np.zeros((1, self.n))
        d_h_d_rel = -d - (b * d - rel) / root
        H[0, P_SL] = d_h_d_rel
        if self.attitude_enabled:
            d_h_d_dir = -rel * (1.0 + b / root)
            H[0, TH_SL] = d_h_d_rel @ (-rot @ skew(arm_body)) + d_h_d_dir @ (
                -rot @ skew(boresight)
            )
        y = np.array([range_meas - h])
        self.last_lidar_nis = self._correct(H, y, np.array([[float(sigma_range) ** 2]]))
        return y

    def update_startracker(
        self,
        q_meas_sensor_to_eci,
        t_meas,
        q_body_from_sensor,
        sigma_sensor_rad,
        earth_rotation_ecef_radps=None,
        latency_s=0.0,
        gate_nis=None,
    ):
        """Absolute attitude update from a star tracker quaternion (sensor -> ECI).

        t_meas              epoch of the attitude (time of exposure), for ECI -> ECEF
        q_body_from_sensor  known mounting quaternion (sensor -> body)
        sigma_sensor_rad    1-sigma noise about the sensor x, y, z axes (z = boresight)
        latency_s           age of the measurement: it is propagated to the current
                            time with the last estimated body rate
        gate_nis            reject if the 3-dof NIS exceeds this value (None = off)

        Measurement model: conj(q) (x) q_meas_body = exp(dtheta + n_body), so
        H = [0 ... I(dtheta) ... 0] and R = R_bs diag(sigma^2) R_bs^T (+ latency term).
        """
        if not self.attitude_enabled:
            return None
        if earth_rotation_ecef_radps is None:
            earth_rotation_ecef_radps = (
                np.zeros(3)
                if self.environment is None
                else self.environment.earth_rotation_ecef_radps
            )
        q_eci_from_ecef = earth_rotation_quaternion(t_meas, earth_rotation_ecef_radps)
        q_body_meas = quaternion_product(
            quaternion_product(
                quaternion_conjugate(q_eci_from_ecef), q_meas_sensor_to_eci
            ),
            quaternion_conjugate(normalize_quaternion(q_body_from_sensor)),
        )
        tau = float(latency_s)
        delta = np.array([1., 0, 0, 0])
        maps = {}
        noise = np.zeros((3, 3))
        if tau > 0.0 and self._omega_prev is not None:
            delta, maps, noise = self._gyro_increment(t_meas, t_meas + tau)
            q_body_meas = quaternion_product(q_body_meas, delta)
        y = quaternion_to_rotation_vector(
            quaternion_product(quaternion_conjugate(self.q), q_body_meas)
        )
        r_bs = body_to_ecef_matrix(delta).T @ body_to_ecef_matrix(q_body_from_sensor)
        R = r_bs @ np.diag(np.asarray(sigma_sensor_rad, float) ** 2) @ r_bs.T
        R = R + noise
        H = np.zeros((3, self.n))
        H[:, TH_SL] = np.eye(3)
        for index, mapping in maps.items():
            H[:, index:index + mapping.shape[1]] = mapping
        S = H @ self.P @ H.T + R
        nis = float(y @ np.linalg.solve(S, y))
        accepted = gate_nis is None or nis <= gate_nis
        self.last_st_nis = nis
        self.startracker_diagnostics.append(
            (float(t_meas + tau), nis, bool(accepted), float(np.linalg.norm(y)))
        )
        if not accepted:
            self.startracker_rejected += 1
            return None
        self.startracker_accepted += 1
        self.last_st_nis = self._correct(H, y, R)
        return y

    def record(self, t):
        self.log.append((t, self.x.copy(), np.sqrt(np.diag(self.P))))

    def to_arrays(self):
        return (
            np.array([e[0] for e in self.log]),
            np.array([e[1] for e in self.log]),
            np.array([e[2] for e in self.log]),
        )
