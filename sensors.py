"""IMU model: truth from dynamics.acceleration_truth + error model.

meas = (1 + scale_factor) * truth + bias_const + bias_random_walk + white_noise
"""

from dataclasses import dataclass, field
import numpy as np

from dynamics import acceleration_truth
from frames import vector3


def _axes3(value):
    """Scalar or 3-vector -> per-axis array (3,)."""
    return np.broadcast_to(np.asarray(value, dtype=float), (3,)).copy()


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

    def __post_init__(self):
        self.station_body_m = vector3(self.station_body_m).copy()
        for name in ("accel_bias_mps2", "gyro_bias_radps", "accel_scale_factor", "gyro_scale_factor"):
            setattr(self, name, _axes3(getattr(self, name)))
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
        truth = acceleration_truth(state, loads, properties, s.station_body_m)
        f_true = truth["specific_force_body_mps2"]  # accelerometer truth
        w_true = state.omega_body_radps  # gyro truth
        f_meas = ((1.0 + s.accel_scale_factor) * f_true + s.accel_bias_mps2 + self.accel_rw
                  + self.rng.normal(0.0, s.accel_sigma_mps2, 3))
        w_meas = ((1.0 + s.gyro_scale_factor) * w_true + s.gyro_bias_radps + self.gyro_rw
                  + self.rng.normal(0.0, s.gyro_sigma_radps, 3))
        # random walk step after sample, bias starts at 0
        self.accel_rw += self.rng.normal(0.0, s.accel_bias_rw * np.sqrt(dt), 3)
        self.gyro_rw += self.rng.normal(0.0, s.gyro_bias_rw * np.sqrt(dt), 3)
        return IMUSample(t, f_meas, w_meas, f_true.copy(), w_true.copy())

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

    # 1) deterministic: no noise, no random walk -> exact scale factor + bias
    spec = IMUSpec(station, accel_sigma_mps2=0.0, gyro_sigma_radps=0.0,
                   accel_scale_factor=1e-3, accel_bias_mps2=[0.05, 0.0, 0.0])
    ts, f_meas, w_meas, f_true, w_true = run(spec)
    assert np.allclose(f_true[:, 0], f_expected), "specific force mismatch"
    assert np.allclose(f_meas[:, 0], f_expected * 1.001 + 0.05), "scale factor / bias mismatch"
    print(f"deterministic OK: f_meas x = {f_meas[0, 0]:.4f} (expected {f_expected * 1.001 + 0.05:.4f})")

    # 2) full error model, PLACEHOLDER values
    spec = IMUSpec(station, accel_scale_factor=1e-3, accel_bias_mps2=[0.05, 0.0, 0.0],
                   gyro_bias_radps=[1e-4, 0.0, 0.0], accel_bias_rw=1e-3, gyro_bias_rw=1e-5)
    ts, f_meas, w_meas, f_true, w_true = run(spec, duration_s=60.0)
    err_f, err_w = f_meas - f_true, w_meas - w_true
    print(f"samples: {len(ts)}")
    print(f"accel err x: start {err_f[:100, 0].mean():.4f}, end {err_f[-100:, 0].mean():.4f} (drift = random walk)")
    print(f"gyro err x: start {err_w[:100, 0].mean():.2e}, end {err_w[-100:, 0].mean():.2e}")
    print("OK")