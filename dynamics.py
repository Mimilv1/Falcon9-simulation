"""A-owned six-DOF core: all shared vectors in Body, +x toward nose.

Force already includes propulsive momentum/pressure reaction. Explicit exhaust
ports supply convective and rotational outflow terms, not a second thrust.
The open-system equations assume axial CG migration and balanced axial outflow.
"""

from dataclasses import dataclass
from typing import Callable
import numpy as np

from frames import (
    vector3, normalize_quaternion, quaternion_rate_body_to_initial_body,
)
from mass_flow import ExhaustPort, point_gyration, validate_ports, PropellantLimitError


@dataclass(frozen=True)
class State:
    position_body_m: np.ndarray  # initial fixed origin -> CG, expressed in CURRENT Body
    velocity_body_mps: np.ndarray
    quaternion_body_to_initial_body_wxyz: np.ndarray  # current Body -> initial Body basis
    omega_body_radps: np.ndarray
    propellant_kg: np.ndarray  # [Stage 1, attached Stage 2]

    def __post_init__(self):
        for field in ("position_body_m", "velocity_body_mps", "omega_body_radps"):
            object.__setattr__(self, field, vector3(getattr(self, field)).copy())
        object.__setattr__(self, "quaternion_body_to_initial_body_wxyz",
                           normalize_quaternion(self.quaternion_body_to_initial_body_wxyz))
        prop = np.asarray(self.propellant_kg, dtype=float)
        if prop.shape != (2,) or not np.all(np.isfinite(prop)) or np.any(prop < -1e-7):
            raise ValueError("State needs nonnegative [Stage 1, Stage 2] propellant masses.")
        object.__setattr__(self, "propellant_kg", np.maximum(prop, 0.0))

    def to_vector(self):
        return np.r_[self.position_body_m, self.velocity_body_mps,
                     self.quaternion_body_to_initial_body_wxyz, self.omega_body_radps,
                     self.propellant_kg]

    @classmethod
    def from_vector(cls, vector):
        y = np.asarray(vector, dtype=float)
        if y.shape != (15,) or not np.all(np.isfinite(y)):
            raise ValueError("State vector must contain 15 finite values.")
        return cls(y[:3], y[3:6], y[6:10], y[10:13], y[13:15])


@dataclass(frozen=True)
class Loads:
    force_body_n: np.ndarray       # Thrust + aerodynamics; gravity excluded.
    moment_body_nm: np.ndarray     # About current CG, Body components.
    gravity_body_mps2: np.ndarray   # Current Body gravity, added exactly once by A.
    consumption_kgps: np.ndarray   # Nonnegative [Stage 1, Stage 2] consumption.
    exhaust_ports: tuple[ExhaustPort, ...] = ()

    def __post_init__(self):
        for field in ("force_body_n", "moment_body_nm", "gravity_body_mps2"):
            object.__setattr__(self, field, vector3(getattr(self, field)).copy())
        rate = np.asarray(self.consumption_kgps, dtype=float)
        if rate.shape != (2,) or not np.all(np.isfinite(rate)) or np.any(rate < 0):
            raise ValueError("Consumption is a nonnegative two-component vector in kg/s.")
        object.__setattr__(self, "consumption_kgps", rate.copy())
        ports = tuple(self.exhaust_ports)
        validate_ports(ports,rate)
        object.__setattr__(self,'exhaust_ports',ports)

    @property
    def consumption_rate_kgps2(self):
        return validate_ports(self.exhaust_ports,self.consumption_kgps)


LoadProvider = Callable[[float, State], Loads]


def state_derivative(time_s, state_vector, model, configuration, provider):
    if hasattr(model,'propellant_bounds'):
        lower,upper = model.propellant_bounds(configuration)
        inventory = np.asarray(state_vector)[13:15]
        if np.any(inventory < np.asarray(lower)-1e-7) or np.any(inventory > np.asarray(upper)+1e-7):
            raise PropellantLimitError('RK trial state crossed a propellant limit.')
    state = State.from_vector(state_vector)
    loads = provider(time_s, state)
    if np.any(loads.consumption_rate_kgps2):
        properties = model.mass_properties(configuration,state.propellant_kg,loads.consumption_kgps,
                                            consumption_rate_kgps2=loads.consumption_rate_kgps2)
    else:
        properties = model.mass_properties(configuration,state.propellant_kg,loads.consumption_kgps)
    omega = state.omega_body_radps
    position_rate = state.velocity_body_mps-np.cross(omega,state.position_body_m)
    acceleration_body,angular_rate = accelerations(state,loads,properties)
    body_velocity_rate = acceleration_body-np.cross(omega,state.velocity_body_mps)
    attitude_rate = quaternion_rate_body_to_initial_body(state.quaternion_body_to_initial_body_wxyz, omega)
    return np.r_[position_rate, body_velocity_rate, attitude_rate, angular_rate,
                 -loads.consumption_kgps]


def accelerations(state,loads,properties):
    """CG inertial acceleration in Body components, plus Body angular acceleration.

    Signed mdot=-q. See docs/powered_dynamics.md for the control-volume balance,
    fixed-port/axial-balance assumptions, and distinction from IMU acceleration.
    """
    omega = state.omega_body_radps
    mass = properties.mass_kg
    migration = (np.zeros(3) if properties.cg_rate_body_mps is None else properties.cg_rate_body_mps)
    if not np.allclose(properties.cg_station_body_m[1:],0,atol=1e-10) or not np.allclose(migration[1:],0,atol=1e-10):
        raise ValueError('This mass-flow model requires axial CG migration.')
    if not np.isclose(properties.mass_rate_kgps,-np.sum(loads.consumption_kgps),rtol=1e-10,atol=1e-9):
        raise ValueError('Mass-property derivative and load consumption disagree.')
    if not np.isclose(properties.mass_acceleration_kgps2,-np.sum(loads.consumption_rate_kgps2),rtol=1e-10,atol=1e-9):
        raise ValueError('Mass-property second derivative and exhaust data disagree.')
    correction = np.zeros(3)
    angular_flux = np.zeros((3,3))
    for port in loads.exhaust_ports:
        arm = properties.arm_from_station(port.station_body_m)
        mdot = -port.mass_flow_kgps
        mddot = -port.mass_flow_rate_kgps2
        correction += 2*mdot*(np.cross(omega,arm)-migration)+mddot*arm
        angular_flux += mdot*(point_gyration(arm)+port.gyration_body_m2)
    acceleration = loads.force_body_n/mass+loads.gravity_body_mps2+correction/mass
    angular = np.linalg.solve(properties.inertia_body_kgm2,
        loads.moment_body_nm-np.cross(omega,properties.inertia_body_kgm2@omega)
        - properties.inertia_rate_body_kgm2ps@omega+angular_flux@omega)
    return acceleration,angular


def rk4_step(time_s, state, step_s, model, configuration, provider):
    """One smooth segment; caller splits steps at cutoff, separation, or fuel limits."""
    if not np.isfinite(step_s) or step_s <= 0:
        raise ValueError("Integration step must be finite and positive.")
    y = state.to_vector()
    f = lambda t, v: state_derivative(t, v, model, configuration, provider)
    k1 = f(time_s, y)
    k2 = f(time_s+step_s/2, y+step_s*k1/2)
    k3 = f(time_s+step_s/2, y+step_s*k2/2)
    k4 = f(time_s+step_s, y+step_s*k3)
    return State.from_vector(y+step_s*(k1+2*k2+2*k3+k4)/6)


def navigation_truth(state, loads):
    """Truth at the instantaneous CG for C; a real displaced IMU needs lever-arm terms."""
    return {
        "position_body_m": state.position_body_m.copy(),
        "velocity_body_mps": state.velocity_body_mps.copy(),
        "quaternion_body_to_initial_body_wxyz": state.quaternion_body_to_initial_body_wxyz.copy(),
        # Angular velocity relative to the fixed initial basis; Earth rotation omitted.
        "omega_body_radps": state.omega_body_radps.copy(),
        "force_nongrav_body_n": loads.force_body_n.copy(),
        "gravity_body_mps2": loads.gravity_body_mps2.copy(),
    }


def acceleration_truth(state, loads, properties, sensor_station_body_m=None):
    """Geometric CG acceleration and specific force at a fixed structural station.

    Supply the same sensor station at every time for a physically fixed IMU.
    Default: fixed-body point coincident with CG at this instant.
    """
    acceleration,angular = accelerations(state,loads,properties)
    station = properties.cg_station_body_m if sensor_station_body_m is None else vector3(sensor_station_body_m)
    arm = station-properties.cg_station_body_m
    migration = np.zeros(3) if properties.cg_rate_body_mps is None else properties.cg_rate_body_mps
    migration_acceleration = np.zeros(3) if properties.cg_acceleration_body_mps2 is None else properties.cg_acceleration_body_mps2
    sensor_acceleration = (acceleration+np.cross(angular,arm)
                           +np.cross(state.omega_body_radps,np.cross(state.omega_body_radps,arm))
                           -2*np.cross(state.omega_body_radps,migration)-migration_acceleration)
    specific_force = sensor_acceleration-loads.gravity_body_mps2
    return {
        "specific_force_body_mps2": specific_force,
        "acceleration_body_mps2": acceleration,
        "sensor_acceleration_body_mps2": sensor_acceleration,
        "angular_acceleration_body_radps2": angular,
    }
