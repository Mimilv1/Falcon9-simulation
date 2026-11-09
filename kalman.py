import numpy as np

N = 9
P_SL, V_SL, B_SL = slice(0, 3), slice(3, 6), slice(6, 9)


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
    ):
        self.x = np.r_[
            np.asarray(p0, float),
            np.asarray(v0, float),
            np.zeros(3) if ba0 is None else np.asarray(ba0, float),
        ]
        self.P = np.diag(np.asarray(sigma0, dtype=float) ** 2)
        self.g = np.asarray(gravity_ecef, dtype=float)
        self.environment = environment
        self.accel_psd = float(accel_psd)
        self.accel_bias_rw = float(accel_bias_rw)
        self.last_accel = np.zeros(3)
        self.last_nis = np.nan
        self.last_lidar_nis = np.nan
        self.log = []

    @classmethod
    def from_imu_spec(cls, imu_spec, p0, v0, sigma0, **kwargs):
        return cls(
            p0,
            v0,
            sigma0,
            accel_psd=imu_spec.accel_sigma_mps2**2 / imu_spec.rate_hz,
            accel_bias_rw=imu_spec.accel_bias_rw,
            **kwargs
        )

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

    def predict(self, f_meas, rot, dt):
        p, v, ba = self.x[P_SL].copy(), self.x[V_SL].copy(), self.x[B_SL].copy()
        model_a, d_dp, d_dv = self._model_acceleration(p, v)
        a = rot @ (f_meas - ba) + model_a
        self.x[P_SL] = p + v * dt + 0.5 * a * dt**2
        self.x[V_SL] = v + a * dt
        self.last_accel = a

        eye = np.eye(3)
        F = np.eye(N)
        F[P_SL, P_SL] = eye + 0.5 * d_dp * dt**2
        F[P_SL, V_SL] = eye * dt + 0.5 * d_dv * dt**2
        F[V_SL, P_SL] = d_dp * dt
        F[V_SL, V_SL] = eye + d_dv * dt
        F[P_SL, B_SL] = -0.5 * rot * dt**2
        F[V_SL, B_SL] = -rot * dt
        S = self.accel_psd
        Q = np.zeros((N, N))
        Q[P_SL, P_SL] = S * dt**3 / 3 * eye
        Q[P_SL, V_SL] = Q[V_SL, P_SL] = S * dt**2 / 2 * eye
        Q[V_SL, V_SL] = S * dt * eye
        Q[B_SL, B_SL] = self.accel_bias_rw**2 * dt * eye
        self.P = F @ self.P @ F.T + Q

    def _correct(self, H, y, R):
        S = H @ self.P @ H.T + R
        K = np.linalg.solve(S, H @ self.P).T
        self.x = self.x + K @ y
        I_KH = np.eye(N) - K @ H
        self.P = I_KH @ self.P @ I_KH.T + K @ R @ K.T
        self.P = 0.5 * (self.P + self.P.T)
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
        a, tau = self.last_accel, latency_s
        z_p = pos_meas + vel_meas * tau + 0.5 * a * tau**2
        z_v = vel_meas + a * tau
        z = np.r_[z_p - rot @ arm_body, z_v - rot @ np.cross(omega, arm_body)]

        H = np.zeros((6, N))
        H[0:3, P_SL] = np.eye(3)
        H[3:6, V_SL] = np.eye(3)
        R = np.diag(
            np.r_[np.asarray(sigma_pos, float) ** 2, np.asarray(sigma_vel, float) ** 2]
        )
        y = z - H @ self.x
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
        rel = self.x[P_SL] + rot @ arm_body - np.asarray(planet_center, float)
        d = rot @ np.asarray(boresight_body, float)
        b = float(rel @ d)
        disc = b * b - (float(rel @ rel) - planet_radius**2)
        if b >= 0.0 or disc <= 0.0:
            return None
        root = np.sqrt(disc)
        h = -b - root
        H = np.zeros((1, N))
        H[0, P_SL] = -d - (b * d - rel) / root
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
