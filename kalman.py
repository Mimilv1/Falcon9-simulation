"""Navigation filter: error-state EKF with optional multiplicative attitude (MEKF).

Nominal state
    x[0:3]  position, ECEF [m]            (IMU station)
    x[3:6]  velocity, ECEF [m/s]
    x[6:9]  accelerometer bias, body [m/s^2]
    q       attitude quaternion body -> ECEF (Hamilton, [w, x, y, z])  [MEKF only]
    bg      gyro bias, body [rad/s]                                      [MEKF only]

Error state (true = nominal + error, attitude error is multiplicative):
    dx = [dp(3), dv(3), dba(3), dtheta(3), dbg(3)]
    q_true = q (x) exp(dtheta)       (body-frame error rotation)
    bg_true = bg + dbg

The attitude block is active only when an initial quaternion `q0` is given.
Without `q0` the filter is the original 9-state filter and the attitude must
be supplied by the caller (`rot` arguments), exactly as before.
"""

import numpy as np

from frames import (
    body_to_ecef_matrix,
    normalize_quaternion,
    quaternion_from_rotation_vector,
    quaternion_product,
    skew,
)

N_NAV = 9
N_ATT = 15
N = N_NAV  # backward-compatible alias
P_SL, V_SL, B_SL = slice(0, 3), slice(3, 6), slice(6, 9)
TH_SL, BG_SL = slice(9, 12), slice(12, 15)


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
    ):
        self.x = np.r_[
            np.asarray(p0, float),
            np.asarray(v0, float),
            np.zeros(3) if ba0 is None else np.asarray(ba0, float),
        ]
        self.attitude_enabled = q0 is not None
        sigma0 = np.asarray(sigma0, dtype=float)
        if self.attitude_enabled:
            self.n = N_ATT
            self.q = normalize_quaternion(q0)
            self.bg = np.zeros(3) if bg0 is None else np.asarray(bg0, dtype=float).copy()
            if sigma0.shape == (N_NAV,):
                sigma0 = np.r_[
                    sigma0,
                    np.broadcast_to(np.asarray(attitude_sigma0, float), (3,)),
                    np.broadcast_to(np.asarray(gyro_bias_sigma0, float), (3,)),
                ]
        else:
            self.n = N_NAV
            self.q = None
            self.bg = np.zeros(3)
        if sigma0.shape != (self.n,):
            raise ValueError(f"sigma0 must have {self.n} entries.")
        self.P = np.diag(sigma0**2)
        self.g = np.asarray(gravity_ecef, dtype=float)
        self.environment = environment
        self.accel_psd = float(accel_psd)
        self.accel_bias_rw = float(accel_bias_rw)
        self.gyro_psd = float(gyro_psd)
        self.gyro_bias_rw = float(gyro_bias_rw)
        self._omega_prev = None
        self.last_accel = np.zeros(3)
        self.last_nis = np.nan
        self.last_lidar_nis = np.nan
        self.log = []

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
        """Propagate with one IMU sample.

        f_meas     accelerometer specific force, body [m/s^2]
        rot        body -> ECEF rotation supplied by the caller; ignored when the
                   attitude is estimated (the filter uses its own quaternion)
        gyro_meas  gyro rate, body [rad/s]; required when the attitude is estimated
        """
        p, v, ba = self.x[P_SL].copy(), self.x[V_SL].copy(), self.x[B_SL].copy()
        model_a, d_dp, d_dv = self._model_acceleration(p, v)

        if self.attitude_enabled:
            if gyro_meas is None:
                raise ValueError("gyro_meas is required when the attitude is estimated.")
            omega = np.asarray(gyro_meas, dtype=float) - self.bg
            omega_avg = omega if self._omega_prev is None else 0.5 * (self._omega_prev + omega)
            # Mid-interval attitude used for the specific-force rotation
            q_mid = normalize_quaternion(
                quaternion_product(
                    self.q, quaternion_from_rotation_vector(0.5 * omega_avg * dt)
                )
            )
            rot = body_to_ecef_matrix(q_mid)

        f_corrected = f_meas - ba
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
        F[P_SL, B_SL] = -0.5 * rot * dt**2
        F[V_SL, B_SL] = -rot * dt
        S = self.accel_psd
        Q = np.zeros((self.n, self.n))
        Q[P_SL, P_SL] = S * dt**3 / 3 * eye
        Q[P_SL, V_SL] = Q[V_SL, P_SL] = S * dt**2 / 2 * eye
        Q[V_SL, V_SL] = S * dt * eye
        Q[B_SL, B_SL] = self.accel_bias_rw**2 * dt * eye

        if self.attitude_enabled:
            # Attitude error enters the acceleration as -R [f]x dtheta
            r_fx = rot @ skew(f_corrected)
            F[V_SL, TH_SL] = -r_fx * dt
            F[P_SL, TH_SL] = -0.5 * r_fx * dt**2
            # dtheta_dot = -[omega]x dtheta - dbg - gyro noise
            dq = quaternion_from_rotation_vector(omega_avg * dt)
            F[TH_SL, TH_SL] = body_to_ecef_matrix(dq).T
            F[TH_SL, BG_SL] = -eye * dt
            Q[TH_SL, TH_SL] = self.gyro_psd * dt * eye
            Q[BG_SL, BG_SL] = self.gyro_bias_rw**2 * dt * eye
            # Nominal attitude propagation over the full interval
            self.q = normalize_quaternion(quaternion_product(self.q, dq))
            self._omega_prev = omega

        self.P = F @ self.P @ F.T + Q

    # ------------------------------------------------------------------ correct
    def _correct(self, H, y, R):
        S = H @ self.P @ H.T + R
        K = np.linalg.solve(S, H @ self.P).T
        dx = K @ y
        self.x = self.x + dx[:N_NAV]
        I_KH = np.eye(self.n) - K @ H
        P = I_KH @ self.P @ I_KH.T + K @ R @ K.T
        if self.attitude_enabled:
            dtheta = dx[TH_SL]
            self.q = normalize_quaternion(
                quaternion_product(self.q, quaternion_from_rotation_vector(dtheta))
            )
            self.bg = self.bg + dx[BG_SL]
            # Covariance reset: express P in the updated attitude error frame
            G = np.eye(self.n)
            G[TH_SL, TH_SL] = np.eye(3) - skew(0.5 * dtheta)
            P = G @ P @ G.T
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
            omega = np.asarray(omega, dtype=float) - self.bg
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
            H[3:6, BG_SL] = rot @ skew(arm_body)
        R = np.diag(
            np.r_[np.asarray(sigma_pos, float) ** 2, np.asarray(sigma_vel, float) ** 2]
        )
        y = z - np.r_[self.x[P_SL], self.x[V_SL]]
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

    def record(self, t):
        self.log.append((t, self.x.copy(), np.sqrt(np.diag(self.P))))

    def to_arrays(self):
        return (
            np.array([e[0] for e in self.log]),
            np.array([e[1] for e in self.log]),
            np.array([e[2] for e in self.log]),
        )