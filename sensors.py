"""Sensors (IMU, GNSS): truth from dynamics + error models.

IMU, truth from dynamics.acceleration_truth:

dT   = T(t) - T_ref
meas = (I + diag(scale_factor + sf_tempco*dT) + misalignment) @ truth
       + bias_const + bias_tempco*dT + bias_random_walk + white_noise
"""

from dataclasses import dataclass, field
from typing import Callable, Optional
import bisect
import numpy as np

from dynamics import acceleration_truth
from frames import vector3


def _axes3(value):
    """Scalar or 3-vector -> per-axis array (3,)."""
    return np.broadcast_to(np.asarray(value, dtype=float), (3,)).copy()


def _cross_axis(value):
    """3x3 cross-axis coupling (dimensionless, rad), diagonal forced to 0."""
    m = np.zeros((3, 3)) if value is None else np.asarray(value, dtype=float).copy()
    if m.shape != (3, 3):
        raise ValueError("misalignment must be a 3x3 matrix")
    np.fill_diagonal(m, 0.0)
    return m


@dataclass
class IMUSpec:
    station_body_m: np.ndarray  # fixed structural station, Body, same datum as cg_station_body_m
    rate_hz: float = 100.0
    # white noise, 1-sigma per sample, PLACEHOLDER (datasheet)
    accel_sigma_mps2: float = 0.01
    gyro_sigma_radps: float = 1e-4
    # constant bias (turn-on), per axis
    accel_bias_mps2: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gyro_bias_radps: np.ndarray = field(default_factory=lambda: np.zeros(3))
    # scale factor error, dimensionless (1e-3 = 1000 ppm), scalar or per axis
    accel_scale_factor: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gyro_scale_factor: np.ndarray = field(default_factory=lambda: np.zeros(3))
    # bias random walk, unit/sqrt(s): accel m/s2/sqrt(s), gyro rad/s/sqrt(s)
    accel_bias_rw: float = 0.0
    gyro_bias_rw: float = 0.0
    # misalignment, 3x3 cross-axis coupling, row = output axis, col = input axis, rad
    accel_misalignment: Optional[np.ndarray] = None
    gyro_misalignment: Optional[np.ndarray] = None
    # temperature, T(t) in deg C, None -> constant at temp_ref_c
    temperature_c: Optional[Callable[[float], float]] = None
    temp_ref_c: float = 20.0
    # temperature coefficients, per deg C, scalar or per axis
    accel_bias_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))  # m/s2/degC
    gyro_bias_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))  # rad/s/degC
    accel_sf_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))  # 1/degC
    gyro_sf_tempco: np.ndarray = field(default_factory=lambda: np.zeros(3))  # 1/degC

    def __post_init__(self):
        self.station_body_m = vector3(self.station_body_m).copy()
        for name in ("accel_bias_mps2", "gyro_bias_radps", "accel_scale_factor", "gyro_scale_factor",
                     "accel_bias_tempco", "gyro_bias_tempco", "accel_sf_tempco", "gyro_sf_tempco"):
            setattr(self, name, _axes3(getattr(self, name)))
        self.accel_misalignment = _cross_axis(self.accel_misalignment)
        self.gyro_misalignment = _cross_axis(self.gyro_misalignment)
        if self.temperature_c is None:
            ref = self.temp_ref_c
            self.temperature_c = lambda t: ref
        if self.rate_hz <= 0:
            raise ValueError("rate_hz must be positive")

    @classmethod
    def from_density(cls, station_body_m, rate_hz, accel_density, gyro_density, **kwargs):
        """Noise density (unit/sqrt(Hz)) -> per-sample sigma, bandwidth = rate/2."""
        bandwidth = rate_hz / 2.0
        return cls(station_body_m, rate_hz,
                   accel_density * np.sqrt(bandwidth),
                   gyro_density * np.sqrt(bandwidth), **kwargs)


@dataclass(frozen=True)
class IMUSample:
    t: float
    accel_meas_mps2: np.ndarray  # specific force, Body
    gyro_meas_radps: np.ndarray  # angular rate, Body
    accel_true_mps2: np.ndarray
    gyro_true_radps: np.ndarray
    temperature_c: float = 20.0


class IMU:
    def __init__(self, spec, seed=None):
        self.spec = spec
        self.rng = np.random.default_rng(seed)
        self.period_s = 1.0 / spec.rate_hz
        self.next_t = 0.0  # next sample time
        self.accel_rw = np.zeros(3)  # bias random walk state
        self.gyro_rw = np.zeros(3)
        self.log = []

    def measure(self, t, state, loads, properties, dt=None):
        """One noisy sample, no sampling-rate logic; advances random walk by dt."""
        s = self.spec
        dt = self.period_s if dt is None else dt
        temp = float(s.temperature_c(t))
        d_temp = temp - s.temp_ref_c
        truth = acceleration_truth(state, loads, properties, s.station_body_m)
        f_true = truth["specific_force_body_mps2"]  # accelerometer truth
        w_true = state.omega_body_radps  # gyro truth
        # scale factor + cross-axis matrix, temperature-dependent scale factor
        mat_a = np.eye(3) + np.diag(s.accel_scale_factor + s.accel_sf_tempco * d_temp) + s.accel_misalignment
        mat_g = np.eye(3) + np.diag(s.gyro_scale_factor + s.gyro_sf_tempco * d_temp) + s.gyro_misalignment
        f_meas = (mat_a @ f_true + s.accel_bias_mps2 + s.accel_bias_tempco * d_temp + self.accel_rw
                  + self.rng.normal(0.0, s.accel_sigma_mps2, 3))
        w_meas = (mat_g @ w_true + s.gyro_bias_radps + s.gyro_bias_tempco * d_temp + self.gyro_rw
                  + self.rng.normal(0.0, s.gyro_sigma_radps, 3))
        # random walk step after sample, bias starts at 0
        self.accel_rw += self.rng.normal(0.0, s.accel_bias_rw * np.sqrt(dt), 3)
        self.gyro_rw += self.rng.normal(0.0, s.gyro_bias_rw * np.sqrt(dt), 3)
        return IMUSample(t, f_meas, w_meas, f_true.copy(), w_true.copy(), temp)

    def update(self, t, state, loads, properties):
        """Call on each ACCEPTED step; returns sample when due, else None. Needs step <= period."""
        if t + 1e-9 < self.next_t:
            return None
        sample = self.measure(t, state, loads, properties)
        self.log.append(sample)
        self.next_t += self.period_s
        return sample

    def to_arrays(self):
        """Log -> (t, accel_meas, gyro_meas, accel_true, gyro_true) numpy arrays."""
        return (np.array([s.t for s in self.log]),
                np.array([s.accel_meas_mps2 for s in self.log]),
                np.array([s.gyro_meas_radps for s in self.log]),
                np.array([s.accel_true_mps2 for s in self.log]),
                np.array([s.gyro_true_radps for s in self.log]))


def truth_inputs(t, state, model, configuration, provider):
    """Rebuild loads + mass properties like dynamics.state_derivative."""
    loads = provider(t, state)
    if np.any(loads.consumption_rate_kgps2):
        properties = model.mass_properties(configuration, state.propellant_kg, loads.consumption_kgps,
                                           consumption_rate_kgps2=loads.consumption_rate_kgps2)
    else:
        properties = model.mass_properties(configuration, state.propellant_kg, loads.consumption_kgps)
    return loads, properties


def _quat_to_matrix(q):
    """Quaternion wxyz (Body -> initial Body) -> rotation matrix."""
    w, x, y, z = np.asarray(q, dtype=float) / np.linalg.norm(q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


@dataclass
class GNSSSpec:
    antenna_station_body_m: np.ndarray  # fixed structural station, Body, same datum as cg_station_body_m
    rate_hz: float = 1.0  # typical 1-10 Hz
    latency_s: float = 0.1  # sample at t reports truth at t - latency
    # position error, initial-frame axes (x = vertical at launch), 1-sigma, PLACEHOLDER (datasheet)
    pos_sigma_m: np.ndarray = field(default_factory=lambda: np.array([2.5, 1.5, 1.5]))  # white
    pos_gm_sigma_m: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.6, 0.6]))  # correlated
    gm_tau_s: float = 300.0  # Gauss-Markov correlation time (atmosphere, multipath, ephemeris)
    vel_sigma_mps: np.ndarray = field(default_factory=lambda: np.full(3, 0.1))  # Doppler, white
    first_fix_time_s: float = 0.0  # time to first fix, no valid output before
    dropout_prob: float = 0.0  # random loss of fix, per sample
    outages: tuple = ()  # ((t0, t1), ...) forced no-fix windows

    def __post_init__(self):
        self.antenna_station_body_m = vector3(self.antenna_station_body_m).copy()
        for name in ("pos_sigma_m", "pos_gm_sigma_m", "vel_sigma_mps"):
            setattr(self, name, _axes3(getattr(self, name)))
        self.outages = tuple((float(a), float(b)) for a, b in self.outages)
        if self.rate_hz <= 0 or self.gm_tau_s <= 0 or self.latency_s < 0:
            raise ValueError("rate_hz, gm_tau_s must be positive, latency_s nonnegative")
        if not 0.0 <= self.dropout_prob <= 1.0:
            raise ValueError("dropout_prob must be in [0, 1]")


@dataclass(frozen=True)
class GNSSSample:
    t: float
    valid: bool  # fix available
    pos_meas_m: np.ndarray  # initial frame, NaN if no fix
    vel_meas_mps: np.ndarray  # initial frame, NaN if no fix
    pos_true_m: np.ndarray  # truth at t - latency
    vel_true_mps: np.ndarray
    pos_sigma_m: np.ndarray  # reported 1-sigma (white + correlated), for filter R
    vel_sigma_mps: np.ndarray


class GNSS:
    """GNSS receiver: antenna position/velocity in initial frame.

    Lever arm + CG migration, latency, white + Gauss-Markov position error,
    Doppler velocity noise, TTFF / outages / random dropouts.
    """

    def __init__(self, spec, seed=None):
        self.spec = spec
        self.rng = np.random.default_rng(seed)
        self.period_s = 1.0 / spec.rate_hz
        self.next_t = 0.0  # next sample time
        self.gm = self.rng.normal(0.0, spec.pos_gm_sigma_m)  # correlated error, stationary start
        self._hist_t, self._hist_p, self._hist_v = [], [], []  # truth history for latency
        self.log = []

    def _truth(self, state, properties):
        """Antenna position/velocity, initial frame."""
        rot = _quat_to_matrix(state.quaternion_body_to_initial_body_wxyz)
        arm = self.spec.antenna_station_body_m - properties.cg_station_body_m
        cg_rate = np.zeros(3) if properties.cg_rate_body_mps is None else properties.cg_rate_body_mps
        p_body = state.position_body_m + arm
        v_body = state.velocity_body_mps + np.cross(state.omega_body_radps, arm) - cg_rate
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
        """Truth at t_q, linear interpolation, clamped to history."""
        ts = self._hist_t
        if t_q <= ts[0]:
            return self._hist_p[0], self._hist_v[0]
        i = bisect.bisect_left(ts, t_q)
        if i >= len(ts):
            return self._hist_p[-1], self._hist_v[-1]
        w = (t_q - ts[i - 1]) / (ts[i] - ts[i - 1])
        return ((1 - w) * self._hist_p[i - 1] + w * self._hist_p[i],
                (1 - w) * self._hist_v[i - 1] + w * self._hist_v[i])

    def _available(self, t):
        s = self.spec
        if t < s.first_fix_time_s - 1e-9:
            return False
        if any(t0 - 1e-9 <= t < t1 - 1e-9 for t0, t1 in s.outages):
            return False
        return not (s.dropout_prob > 0.0 and self.rng.random() < s.dropout_prob)

    def measure(self, t, p_true, v_true, dt=None):
        """One sample from given truth, no sampling/latency logic; advances Gauss-Markov by dt."""
        s = self.spec
        dt = self.period_s if dt is None else dt
        pos_meas = p_true + self.gm + self.rng.normal(0.0, s.pos_sigma_m, 3)
        vel_meas = v_true + self.rng.normal(0.0, s.vel_sigma_mps, 3)
        valid = self._available(t)
        if not valid:
            pos_meas = np.full(3, np.nan)
            vel_meas = np.full(3, np.nan)
        # Gauss-Markov step after sample, error evolves even without fix
        phi = np.exp(-dt / s.gm_tau_s)
        self.gm = phi * self.gm + s.pos_gm_sigma_m * np.sqrt(1.0 - phi ** 2) * self.rng.normal(size=3)
        pos_sigma = np.sqrt(s.pos_sigma_m ** 2 + s.pos_gm_sigma_m ** 2)
        return GNSSSample(t, valid, pos_meas, vel_meas, p_true.copy(), v_true.copy(),
                          pos_sigma, s.vel_sigma_mps.copy())

    def update(self, t, state, properties):
        """Call on EVERY accepted step (truth history); returns sample when due, else None."""
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
        """Log -> (t, valid, pos_meas, vel_meas, pos_true, vel_true) numpy arrays."""
        return (np.array([s.t for s in self.log]),
                np.array([s.valid for s in self.log]),
                np.array([s.pos_meas_m for s in self.log]),
                np.array([s.vel_meas_mps for s in self.log]),
                np.array([s.pos_true_m for s in self.log]),
                np.array([s.vel_true_mps for s in self.log]))


if __name__ == "__main__":
    # self-test: fake model, constant thrust along +x, no burn, no ports
    from types import SimpleNamespace
    from dynamics import State, Loads, rk4_step

    THRUST_N, MASS_KG, G = 20000.0, 1000.0, 9.81  # expected specific force = 20 m/s2

    model = SimpleNamespace(
        mass_properties=lambda configuration, propellant, consumption, **kw: SimpleNamespace(
            mass_kg=MASS_KG,
            cg_station_body_m=np.array([20.0, 0.0, 0.0]),
            cg_rate_body_mps=np.zeros(3),
            cg_acceleration_body_mps2=np.zeros(3),
            mass_rate_kgps=0.0,
            mass_acceleration_kgps2=0.0,
            inertia_body_kgm2=np.diag([1e4, 1e6, 1e6]),
            inertia_rate_body_kgm2ps=np.zeros((3, 3)),
            arm_from_station=lambda station: np.asarray(station) - np.array([20.0, 0.0, 0.0]),
        )
    )
    configuration = None
    provider = lambda t, s: Loads(
        force_body_n=[THRUST_N, 0, 0],
        moment_body_nm=[0, 0, 0],
        gravity_body_mps2=[-G, 0, 0],
        consumption_kgps=[0.0, 0.0],
    )

    def run(spec, duration_s=5.0, dt=0.01, seed=0):
        state = State([0, 0, 0], [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0.0, 0.0])
        imu = IMU(spec, seed=seed)
        t = 0.0
        loads, props = truth_inputs(t, state, model, configuration, provider)
        imu.update(t, state, loads, props)  # sample at t=0
        while t < duration_s - 1e-9:
            state = rk4_step(t, state, dt, model, configuration, provider)  # accepted step
            t += dt
            loads, props = truth_inputs(t, state, model, configuration, provider)
            imu.update(t, state, loads, props)
        return imu.to_arrays()

    f_expected = THRUST_N / MASS_KG
    station = [10.0, 0.0, 0.0]
    misalign = np.zeros((3, 3))
    misalign[1, 0] = 1e-3  # x input leaks into y output
    temp_profile = lambda t: 20.0 + 2.0 * t  # heating ramp, PLACEHOLDER

    # 1) deterministic: no noise, no random walk -> exact model
    spec = IMUSpec(station, accel_sigma_mps2=0.0, gyro_sigma_radps=0.0,
                   accel_scale_factor=1e-3, accel_bias_mps2=[0.05, 0.0, 0.0],
                   accel_misalignment=misalign,
                   temperature_c=temp_profile, temp_ref_c=20.0,
                   accel_bias_tempco=[0.002, 0.0, 0.0], accel_sf_tempco=1e-5)
    ts, f_meas, w_meas, f_true, w_true = run(spec)
    d_temp = 2.0 * ts
    assert np.allclose(f_true[:, 0], f_expected), "specific force mismatch"
    expected_x = f_expected * (1 + 1e-3 + 1e-5 * d_temp) + 0.05 + 0.002 * d_temp
    assert np.allclose(f_meas[:, 0], expected_x), "scale factor / bias / tempco mismatch"
    assert np.allclose(f_meas[:, 1], 1e-3 * f_expected), "misalignment mismatch"
    print(f"deterministic OK: f_meas x {f_meas[0, 0]:.4f} -> {f_meas[-1, 0]:.4f}, y = {f_meas[0, 1]:.4f}")

    # 2) full error model, PLACEHOLDER values
    spec = IMUSpec(station, accel_scale_factor=1e-3, accel_bias_mps2=[0.05, 0.0, 0.0],
                   gyro_bias_radps=[1e-4, 0.0, 0.0], accel_bias_rw=1e-3, gyro_bias_rw=1e-5,
                   accel_misalignment=misalign, temperature_c=temp_profile,
                   accel_bias_tempco=[0.002, 0.0, 0.0], accel_sf_tempco=1e-5)
    ts, f_meas, w_meas, f_true, w_true = run(spec, duration_s=60.0)
    err_f, err_w = f_meas - f_true, w_meas - w_true
    print(f"samples: {len(ts)}")
    print(f"accel err x: start {err_f[:100, 0].mean():.4f}, end {err_f[-100:, 0].mean():.4f}")
    print(f"gyro err x: start {err_w[:100, 0].mean():.2e}, end {err_w[-100:, 0].mean():.2e}")
    print("OK")

    # ---- GNSS ----
    def run_gnss(spec, duration_s=5.0, dt=0.01, seed=0):
        state = State([0, 0, 0], [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0.0, 0.0])
        gnss = GNSS(spec, seed=seed)
        t = 0.0
        loads, props = truth_inputs(t, state, model, configuration, provider)
        gnss.update(t, state, props)  # sample at t=0
        while t < duration_s - 1e-9:
            state = rk4_step(t, state, dt, model, configuration, provider)  # accepted step
            t += dt
            loads, props = truth_inputs(t, state, model, configuration, provider)
            gnss.update(t, state, props)
        return gnss.to_arrays()

    # quaternion convention: 90 deg about z maps x -> y
    q90 = [np.cos(np.pi / 4), 0.0, 0.0, np.sin(np.pi / 4)]
    assert np.allclose(_quat_to_matrix(q90) @ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), "quaternion convention"

    a_x = THRUST_N / MASS_KG - G  # inertial acceleration, gravity included
    antenna = [30.0, 0.0, 0.0]  # arm = +10 m from CG (cg station 20 m)

    # 3) deterministic: no noise, latency 0.2 s, lever arm +10 m
    spec = GNSSSpec(antenna, rate_hz=1.0, latency_s=0.2, pos_sigma_m=0.0, pos_gm_sigma_m=0.0, vel_sigma_mps=0.0)
    ts, valid, pos_m, vel_m, pos_t, vel_t = run_gnss(spec)
    t_q = np.maximum(ts - 0.2, 0.0)
    assert valid.all()
    assert np.allclose(pos_m[:, 0], 0.5 * a_x * t_q ** 2 + 10.0, atol=1e-3), "gnss position mismatch"
    assert np.allclose(vel_m[:, 0], a_x * t_q, atol=1e-3), "gnss velocity mismatch"
    print(f"gnss deterministic OK: pos x {pos_m[0, 0]:.3f} -> {pos_m[-1, 0]:.3f}")

    # 4) fix logic + full error model, PLACEHOLDER values
    spec = GNSSSpec(antenna, rate_hz=1.0, first_fix_time_s=2.0, outages=((10.0, 12.0),))
    ts, valid, pos_m, vel_m, pos_t, vel_t = run_gnss(spec, duration_s=60.0)
    expected_valid = ~((ts < 2.0 - 1e-9) | ((ts >= 10.0 - 1e-9) & (ts < 12.0 - 1e-9)))
    assert np.array_equal(valid, expected_valid), "fix flags mismatch"
    assert np.isnan(pos_m[~valid]).all() and np.isfinite(pos_m[valid]).all()
    print(f"gnss samples: {len(ts)}, valid: {valid.sum()}")
    print(f"pos err std (m): {np.round((pos_m - pos_t)[valid].std(axis=0), 2)}")
    print(f"vel err std (m/s): {np.round((vel_m - vel_t)[valid].std(axis=0), 3)}")
    print("GNSS OK")