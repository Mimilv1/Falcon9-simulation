"""A-owned six-DOF core: Body +x nose, fixed ENU position, no Earth rotation.

Force already includes propulsive momentum reaction. Do not also add an
exhaust-velocity times mass-flow term. The preliminary angular balance retains
Jdot but assumes no separate angular-momentum outflow or jet-damping term.
"""

from dataclasses import dataclass
from typing import Callable
import numpy as np

from frames import (
    vector3, normalize_quaternion, body_to_enu_matrix, quaternion_rate_body_to_enu,
)


@dataclass(frozen=True)
class State:
    position_enu_m: np.ndarray
    velocity_body_mps: np.ndarray
    quaternion_body_to_enu_wxyz: np.ndarray
    omega_body_radps: np.ndarray
    propellant_kg: np.ndarray  # [Stage 1, attached Stage 2]

    def __post_init__(self):
        for field in ("position_enu_m", "velocity_body_mps", "omega_body_radps"):
            object.__setattr__(self, field, vector3(getattr(self, field)).copy())
        object.__setattr__(self, "quaternion_body_to_enu_wxyz",
                           normalize_quaternion(self.quaternion_body_to_enu_wxyz))
        prop = np.asarray(self.propellant_kg, dtype=float)
        if prop.shape != (2,) or not np.all(np.isfinite(prop)) or np.any(prop < -1e-7):
            raise ValueError("State needs nonnegative [Stage 1, Stage 2] propellant masses.")
        object.__setattr__(self, "propellant_kg", np.maximum(prop, 0.0))

    def to_vector(self):
        return np.r_[self.position_enu_m, self.velocity_body_mps,
                     self.quaternion_body_to_enu_wxyz, self.omega_body_radps,
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
    gravity_enu_mps2: np.ndarray   # Gravity added exactly once by A.
    consumption_kgps: np.ndarray   # Nonnegative [Stage 1, Stage 2] consumption.

    def __post_init__(self):
        for field in ("force_body_n", "moment_body_nm", "gravity_enu_mps2"):
            object.__setattr__(self, field, vector3(getattr(self, field)).copy())
        rate = np.asarray(self.consumption_kgps, dtype=float)
        if rate.shape != (2,) or not np.all(np.isfinite(rate)) or np.any(rate < 0):
            raise ValueError("Consumption is a nonnegative two-component vector in kg/s.")
        object.__setattr__(self, "consumption_kgps", rate.copy())


LoadProvider = Callable[[float, State], Loads]


def state_derivative(time_s, state_vector, model, configuration, provider):
    state = State.from_vector(state_vector)
    loads = provider(time_s, state)
    if model.mass_model_name == "dual_tank" and np.any(loads.consumption_kgps > 0):
        raise ValueError("Powered dual-tank propagation requires a consistent open-system model. "
                         "Use fixed-inventory diagnostics; do not silently use Jdot alone.")
    properties = model.mass_properties(configuration, state.propellant_kg, loads.consumption_kgps)
    rotation = body_to_enu_matrix(state.quaternion_body_to_enu_wxyz)
    omega = state.omega_body_radps
    position_rate = rotation @ state.velocity_body_mps
    # A rotating Body velocity needs this transport term even when Earth does not rotate.
    body_velocity_rate = (
        loads.force_body_n/properties.mass_kg
        + rotation.T @ loads.gravity_enu_mps2
        - np.cross(omega, state.velocity_body_mps)
    )
    attitude_rate = quaternion_rate_body_to_enu(state.quaternion_body_to_enu_wxyz, omega)
    angular_rate = np.linalg.solve(properties.inertia_body_kgm2,
        loads.moment_body_nm
        - np.cross(omega, properties.inertia_body_kgm2 @ omega)
        - properties.inertia_rate_body_kgm2ps @ omega
    )
    return np.r_[position_rate, body_velocity_rate, attitude_rate, angular_rate,
                 -loads.consumption_kgps]


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
    rotation = body_to_enu_matrix(state.quaternion_body_to_enu_wxyz)
    return {
        "position_enu_m": state.position_enu_m.copy(),
        "velocity_enu_mps": rotation @ state.velocity_body_mps,
        "quaternion_body_to_enu_wxyz": state.quaternion_body_to_enu_wxyz.copy(),
        # Body/ENU equals Body/inertial ONLY under the documented zero-Earth-rotation approximation.
        "omega_body_radps": state.omega_body_radps.copy(),
        "force_nongrav_body_n": loads.force_body_n.copy(),
        "gravity_enu_mps2": loads.gravity_enu_mps2.copy(),
    }


def acceleration_truth(state, loads, mass_kg):
    rotation = body_to_enu_matrix(state.quaternion_body_to_enu_wxyz)
    specific_force = loads.force_body_n/mass_kg
    return {
        "specific_force_body_mps2": specific_force,
        "acceleration_enu_mps2": rotation @ specific_force + loads.gravity_enu_mps2,
    }
