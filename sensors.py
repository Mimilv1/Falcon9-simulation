from dataclasses import dataclass, field
from typing import Callable, Optional
import bisect
import numpy as np

from dynamics import acceleration_truth
from frames import (
    vector3,
    normalize_quaternion,
    quaternion_conjugate,
    quaternion_from_rotation_vector,
    quaternion_product,
    body_to_ecef_matrix,
)


def _axes3(value):
    return np.broadcast_to(np.asarray(value, dtype=float), (3,)).copy()


def _cross_axis(value):
    m = np.zeros((3, 3)) if value is None else np.asarray(value, dtype=float).copy()
    if m.shape != (3, 3):
        raise ValueError("misalignment must be a 3x3 matrix")
    np.fill_diagonal(m, 0.0)
    return m


@dataclass
class IMUSpec:
    station_body_m: np.ndarray
    rate_hz: float = 200.0

    accel_sigma_mps2: float = 0.01
    gyro_sigma_radps: float = 1e-4

    accel_bias_mps2: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gyro_bias_radps: np.ndarray = field(default_factory=lambda: np.zeros(3))

    accel_scale_factor: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gyro_scale_factor: np.ndarray = field(default_factory=lambda: np.zeros(3))

    accel_bias_rw: float = 0.0
    gyro_bias_rw: float = 0.0

    accel_misalignment: Optional[np.ndarray] = None
    gyro_misalignment: Optional[np.ndarray] = None

    temperature_c: Optional[Callable[[float], float]] = None
    temp_ref_c: float = 20.0

    accel_bias_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gyro_bias_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))
    accel_sf_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gyro_sf_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self):
        self.station_body_m = vector3(self.station_body_m).copy()
        for name in (
            "accel_bias_mps2",
            "gyro_bias_radps",
            "accel_scale_factor",
            "gyro_scale_factor",
            "accel_bias_tempco",
            "gyro_bias_tempco",
            "accel_sf_tempco",
            "gyro_sf_tempco",
        ):
            setattr(self, name, _axes3(getattr(self, name)))
        self.accel_misalignment = _cross_axis(self.accel_misalignment)
        self.gyro_misalignment = _cross_axis(self.gyro_misalignment)
        if self.temperature_c is None:
            ref = self.temp_ref_c
            self.temperature_c = lambda t: ref
        if self.rate_hz <= 0:
            raise ValueError("rate_hz must be positive")

    @classmethod
    def from_density(
        cls, station_body_m, rate_hz, accel_density, gyro_density, **kwargs
    ):
        bandwidth = rate_hz / 2.0
        return cls(
            station_body_m,
            rate_hz,
            accel_density * np.sqrt(bandwidth),
            gyro_density * np.sqrt(bandwidth),
            **kwargs,
        )


@dataclass(frozen=True)
class IMUSample:
    t: float
    accel_meas_mps2: np.ndarray
    gyro_meas_radps: np.ndarray
    accel_true_mps2: np.ndarray
    gyro_true_radps: np.ndarray
    temperature_c: float = 20.0


class IMU:
    def __init__(self, spec, seed=None):
        self.spec = spec
        self.rng = np.random.default_rng(seed)
        self.period_s = 1.0 / spec.rate_hz
        self.next_t = 0.0
        self.accel_rw = np.zeros(3)
        self.gyro_rw = np.zeros(3)
        self.log = []

    def measure(self, t, state, loads, properties, dt=None):
        s = self.spec
        dt = self.period_s if dt is None else dt
        temp = float(s.temperature_c(t))
        d_temp = temp - s.temp_ref_c
        truth = acceleration_truth(state, loads, properties, s.station_body_m)
        f_true = truth["specific_force_body_mps2"]
        w_true = state.omega_body_radps

        mat_a = (
            np.eye(3)
            + np.diag(s.accel_scale_factor + s.accel_sf_tempco * d_temp)
            + s.accel_misalignment
        )
        mat_g = (
            np.eye(3)
            + np.diag(s.gyro_scale_factor + s.gyro_sf_tempco * d_temp)
            + s.gyro_misalignment
        )
        f_meas = (
            mat_a @ f_true
            + s.accel_bias_mps2
            + s.accel_bias_tempco * d_temp
            + self.accel_rw
            + self.rng.normal(0.0, s.accel_sigma_mps2, 3)
        )
        w_meas = (
            mat_g @ w_true
            + s.gyro_bias_radps
            + s.gyro_bias_tempco * d_temp
            + self.gyro_rw
            + self.rng.normal(0.0, s.gyro_sigma_radps, 3)
        )

        self.accel_rw += self.rng.normal(0.0, s.accel_bias_rw * np.sqrt(dt), 3)
        self.gyro_rw += self.rng.normal(0.0, s.gyro_bias_rw * np.sqrt(dt), 3)
        return IMUSample(t, f_meas, w_meas, f_true.copy(), w_true.copy(), temp)

    def update(self, t, state, loads, properties):
        if t + 1e-9 < self.next_t:
            return None
        sample = self.measure(t, state, loads, properties)
        self.log.append(sample)
        self.next_t += self.period_s
        return sample

    def to_arrays(self):
        return (
            np.array([s.t for s in self.log]),
            np.array([s.accel_meas_mps2 for s in self.log]),
            np.array([s.gyro_meas_radps for s in self.log]),
            np.array([s.accel_true_mps2 for s in self.log]),
            np.array([s.gyro_true_radps for s in self.log]),
        )


def truth_inputs(t, state, model, configuration, provider):
    loads = provider(t, state)
    if np.any(loads.consumption_rate_kgps2):
        properties = model.mass_properties(
            configuration,
            state.propellant_kg,
            loads.consumption_kgps,
            consumption_rate_kgps2=loads.consumption_rate_kgps2,
        )
    else:
        properties = model.mass_properties(
            configuration, state.propellant_kg, loads.consumption_kgps
        )
    return loads, properties


def quat_to_matrix(q):
    w, x, y, z = np.asarray(q, dtype=float) / np.linalg.norm(q)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


@dataclass
class GNSSSpec:
    antenna_station_body_m: np.ndarray
    rate_hz: float = 10.0
    latency_s: float = 0.1

    pos_sigma_m: np.ndarray = field(default_factory=lambda: np.array([1.5, 1.5, 2.5]))
    pos_gm_sigma_m: np.ndarray = field(
        default_factory=lambda: np.array([0.6, 0.6, 1.0])
    )
    gm_tau_s: float = 300.0
    vel_sigma_mps: np.ndarray = field(default_factory=lambda: np.full(3, 0.1))
    first_fix_time_s: float = 0.0
    dropout_prob: float = 0.0
    outages: tuple = ()

    def __post_init__(self):
        self.antenna_station_body_m = vector3(self.antenna_station_body_m).copy()
        for name in ("pos_sigma_m", "pos_gm_sigma_m", "vel_sigma_mps"):
            setattr(self, name, _axes3(getattr(self, name)))
        self.outages = tuple((float(a), float(b)) for a, b in self.outages)
        if self.rate_hz <= 0 or self.gm_tau_s <= 0 or self.latency_s < 0:
            raise ValueError(
                "rate_hz, gm_tau_s must be positive, latency_s nonnegative"
            )
        if not 0.0 <= self.dropout_prob <= 1.0:
            raise ValueError("dropout_prob must be in [0, 1]")


@dataclass(frozen=True)
class GNSSSample:
    t: float
    valid: bool
    pos_meas_m: np.ndarray
    vel_meas_mps: np.ndarray
    pos_true_m: np.ndarray
    vel_true_mps: np.ndarray
    pos_sigma_m: np.ndarray
    vel_sigma_mps: np.ndarray


class GNSS:

    def __init__(self, spec, seed=None):
        self.spec = spec
        self.rng = np.random.default_rng(seed)
        self.period_s = 1.0 / spec.rate_hz
        self.next_t = 0.0
        self.gm = self.rng.normal(0.0, spec.pos_gm_sigma_m)
        self._hist_t, self._hist_p, self._hist_v = [], [], []
        self.log = []

    def _truth(self, state, properties):
        rot = quat_to_matrix(state.quaternion_body_to_ecef_wxyz)
        arm = self.spec.antenna_station_body_m - properties.cg_station_body_m
        cg_rate = (
            np.zeros(3)
            if properties.cg_rate_body_mps is None
            else properties.cg_rate_body_mps
        )
        p_body = state.position_body_m + arm
        v_body = (
            state.velocity_body_mps + np.cross(state.omega_body_radps, arm) - cg_rate
        )
        return rot @ p_body, rot @ v_body

    def _push(self, t, p, v):
        self._hist_t.append(t)
        self._hist_p.append(p)
        self._hist_v.append(v)
        keep_from = t - self.spec.latency_s - self.period_s
        while len(self._hist_t) > 2 and self._hist_t[1] <= keep_from:
            self._hist_t.pop(0)
            self._hist_p.pop(0)
            self._hist_v.pop(0)

    def _delayed(self, t_q):
        ts = self._hist_t
        if t_q <= ts[0]:
            return self._hist_p[0], self._hist_v[0]
        i = bisect.bisect_left(ts, t_q)
        if i >= len(ts):
            return self._hist_p[-1], self._hist_v[-1]
        w = (t_q - ts[i - 1]) / (ts[i] - ts[i - 1])
        return (
            (1 - w) * self._hist_p[i - 1] + w * self._hist_p[i],
            (1 - w) * self._hist_v[i - 1] + w * self._hist_v[i],
        )

    def _available(self, t):
        s = self.spec
        if t < s.first_fix_time_s - 1e-9:
            return False
        if any(t0 - 1e-9 <= t < t1 - 1e-9 for t0, t1 in s.outages):
            return False
        return not (s.dropout_prob > 0.0 and self.rng.random() < s.dropout_prob)

    def measure(self, t, p_true, v_true, dt=None):
        s = self.spec
        dt = self.period_s if dt is None else dt
        pos_meas = p_true + self.gm + self.rng.normal(0.0, s.pos_sigma_m, 3)
        vel_meas = v_true + self.rng.normal(0.0, s.vel_sigma_mps, 3)
        valid = self._available(t)
        if not valid:
            pos_meas = np.full(3, np.nan)
            vel_meas = np.full(3, np.nan)

        phi = np.exp(-dt / s.gm_tau_s)
        self.gm = phi * self.gm + s.pos_gm_sigma_m * np.sqrt(
            1.0 - phi**2
        ) * self.rng.normal(size=3)
        pos_sigma = np.sqrt(s.pos_sigma_m**2 + s.pos_gm_sigma_m**2)
        return GNSSSample(
            t,
            valid,
            pos_meas,
            vel_meas,
            p_true.copy(),
            v_true.copy(),
            pos_sigma,
            s.vel_sigma_mps.copy(),
        )

    def update(self, t, state, properties):
        p, v = self._truth(state, properties)
        self._push(t, p, v)
        if t + 1e-9 < self.next_t:
            return None
        p_d, v_d = self._delayed(t - self.spec.latency_s)
        sample = self.measure(t, p_d, v_d)
        self.log.append(sample)
        self.next_t += self.period_s
        return sample

    def to_arrays(self):
        return (
            np.array([s.t for s in self.log]),
            np.array([s.valid for s in self.log]),
            np.array([s.pos_meas_m for s in self.log]),
            np.array([s.vel_meas_mps for s in self.log]),
            np.array([s.pos_true_m for s in self.log]),
            np.array([s.vel_true_mps for s in self.log]),
        )


@dataclass
class LidarSpec:
    station_body_m: np.ndarray
    boresight_body: np.ndarray = field(
        default_factory=lambda: np.array([-1.0, 0.0, 0.0])
    )
    rate_hz: float = 10.0
    max_range_m: float = 1800.0
    min_range_m: float = 1.0
    planet_radius_m: float = 6_371_000.0

    planet_center_ecef_m: Optional[np.ndarray] = None
    range_sigma_m: float = 0.05
    range_sigma_per_m: float = 1e-4
    range_bias_m: float = 0.0
    range_scale_factor: float = 0.0
    dropout_prob: float = 0.0

    def __post_init__(self):
        self.station_body_m = vector3(self.station_body_m).copy()
        axis = vector3(self.boresight_body).astype(float)
        if np.linalg.norm(axis) == 0.0:
            raise ValueError("boresight_body must be nonzero")
        self.boresight_body = axis / np.linalg.norm(axis)
        if self.planet_center_ecef_m is None:
            self.planet_center_ecef_m = np.zeros(3)
        self.planet_center_ecef_m = vector3(self.planet_center_ecef_m).copy()
        if self.rate_hz <= 0 or self.planet_radius_m <= 0:
            raise ValueError("rate_hz and planet_radius_m must be positive")
        if not 0.0 <= self.min_range_m < self.max_range_m:
            raise ValueError("need 0 <= min_range_m < max_range_m")
        if not 0.0 <= self.dropout_prob <= 1.0:
            raise ValueError("dropout_prob must be in [0, 1]")


@dataclass(frozen=True)
class LidarSample:
    t: float
    in_reach: bool
    valid: bool
    range_meas_m: float
    range_true_m: float
    sigma_m: float


class Lidar:

    def __init__(self, spec, seed=None):
        self.spec = spec
        self.rng = np.random.default_rng(seed)
        self.period_s = 1.0 / spec.rate_hz
        self.next_t = 0.0
        self.log = []

    def _geometry(self, state, properties):
        s = self.spec
        rot = quat_to_matrix(state.quaternion_body_to_ecef_wxyz)
        arm = s.station_body_m - properties.cg_station_body_m
        rel = rot @ (state.position_body_m + arm) - s.planet_center_ecef_m
        return rel, rot @ s.boresight_body

    def measure(self, t, state, properties):
        s = self.spec
        rel, d = self._geometry(state, properties)
        radius = s.planet_radius_m
        dist = float(np.sqrt(rel[0] ** 2 + rel[1] ** 2 + rel[2] ** 2))
        in_reach = dist < radius + s.max_range_m
        range_true = np.nan
        if in_reach and dist > radius:
            b = float(rel @ d)
            disc = b * b - (dist * dist - radius * radius)
            if b < 0.0 and disc >= 0.0:
                hit = -b - np.sqrt(disc)
                if s.min_range_m <= hit <= s.max_range_m:
                    range_true = float(hit)
        geometric = not np.isnan(range_true)
        valid = geometric and not (
            s.dropout_prob > 0.0 and self.rng.random() < s.dropout_prob
        )
        range_meas, sigma = np.nan, np.nan
        if valid:
            sigma = s.range_sigma_m + s.range_sigma_per_m * range_true
            range_meas = (
                (1.0 + s.range_scale_factor) * range_true
                + s.range_bias_m
                + self.rng.normal(0.0, sigma)
            )
        return LidarSample(
            t, bool(in_reach), bool(valid), range_meas, range_true, sigma
        )

    def update(self, t, state, properties):
        if t + 1e-9 < self.next_t:
            return None
        sample = self.measure(t, state, properties)
        self.log.append(sample)
        self.next_t += self.period_s
        return sample

    def to_arrays(self):
        return (
            np.array([s.t for s in self.log]),
            np.array([s.in_reach for s in self.log]),
            np.array([s.valid for s in self.log]),
            np.array([s.range_meas_m for s in self.log]),
            np.array([s.range_true_m for s in self.log]),
        )


# ---------------------------------------------------------------------------
# Gyro-based attitude measurement (quaternion propagation)
#
# Conventions (identical to frames.py / dynamics.py):
#   - Hamilton quaternions, scalar-first [w, x, y, z]
#   - q maps body -> ECEF, omega is expressed in the body frame, so
#     q_dot = 0.5 * q (x) [0, omega_body]
#   - IMU.measure() returns the body rate relative to ECEF, so no Earth-rate
#     term is needed.
# ---------------------------------------------------------------------------


def attitude_error_rotation_vector(q_est, q_true):
    """Body-frame rotation vector [rad] that takes q_true to q_est."""
    q_err = quaternion_product(quaternion_conjugate(q_true), q_est)
    if q_err[0] < 0.0:
        q_err = -q_err  # shortest rotation
    vector_norm = float(np.linalg.norm(q_err[1:]))
    if vector_norm < 1e-15:
        return np.zeros(3)
    return 2.0 * np.arctan2(vector_norm, q_err[0]) * q_err[1:] / vector_norm


@dataclass
class GyroAttitudeSpec:
    # Initial misalignment of the estimate, body frame [rad]
    initial_error_rotvec_rad: np.ndarray = field(default_factory=lambda: np.zeros(3))
    # Constant gyro bias the estimator subtracts (0 = uncompensated gyro)
    gyro_bias_estimate_radps: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self):
        for name in ("initial_error_rotvec_rad", "gyro_bias_estimate_radps"):
            setattr(self, name, _axes3(getattr(self, name)))


@dataclass(frozen=True)
class GyroAttitudeSample:
    t: float
    quaternion_est_wxyz: np.ndarray
    quaternion_true_wxyz: np.ndarray
    error_rotvec_body_rad: np.ndarray
    error_angle_rad: float


class GyroAttitude:
    """Attitude measurement by propagating the quaternion from noisy gyros.

    Each IMU sample advances the estimate with the exact exponential map,
    using the trapezoidal mean of two consecutive gyro samples over the
    interval. There is no external aiding: gyro noise, bias and scale-factor
    errors make the estimate drift (angle random walk + bias * t).
    """

    def __init__(self, spec, q0):
        self.spec = spec
        q0 = normalize_quaternion(q0)
        error = quaternion_from_rotation_vector(spec.initial_error_rotvec_rad)
        self.q = normalize_quaternion(quaternion_product(q0, error))
        self.last_t = None
        self.last_omega = np.zeros(3)
        self.log = []

    @property
    def rotation(self):
        """Estimated body -> ECEF rotation matrix."""
        return body_to_ecef_matrix(self.q)

    def update(self, t, imu_sample, q_true):
        """Advance to the IMU sample time; q_true is used only for error logging."""
        omega = imu_sample.gyro_meas_radps - self.spec.gyro_bias_estimate_radps
        if not np.all(np.isfinite(omega)):
            raise ValueError("Gyro sample must be finite.")
        if self.last_t is not None:
            dt = t - self.last_t
            if dt < 0.0:
                raise ValueError("Gyro samples must be time-ordered.")
            if dt > 0.0:
                dq = quaternion_from_rotation_vector(0.5 * (self.last_omega + omega) * dt)
                self.q = normalize_quaternion(quaternion_product(self.q, dq))
        self.last_t, self.last_omega = float(t), omega.copy()
        q_true = normalize_quaternion(q_true)
        error = attitude_error_rotation_vector(self.q, q_true)
        sample = GyroAttitudeSample(
            t, self.q.copy(), q_true, error, float(np.linalg.norm(error))
        )
        self.log.append(sample)
        return sample

    def to_arrays(self):
        return (
            np.array([s.t for s in self.log]),
            np.array([s.quaternion_est_wxyz for s in self.log]),
            np.array([s.quaternion_true_wxyz for s in self.log]),
            np.array([s.error_rotvec_body_rad for s in self.log]),
        )