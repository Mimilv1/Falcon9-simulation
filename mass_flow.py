"""Explicit exhaust boundary data; no propulsion law or vehicle constants.

Positive q denotes mass leaving a stage. Thrust, including relative exhaust
momentum and pressure reaction, is supplied separately in Loads.force_body_n.
"""
from dataclasses import dataclass, field
import numpy as np
from frames import vector3


class PropellantLimitError(ValueError):
    """An RK trial state crossed an inventory boundary."""


def point_gyration(position):
    """Per-unit-mass parallel-axis tensor: |r|² I - r rᵀ, in m²."""
    r = vector3(position)
    return np.dot(r,r)*np.eye(3)-np.outer(r,r)


def uniform_disk_gyration(radius_m, normal_body=(1,0,0)):
    """Uniform exit mass-flux disk approximation; not a measured flow profile."""
    normal = vector3(normal_body)
    norm = np.linalg.norm(normal)
    if not np.isfinite(radius_m) or radius_m < 0 or norm < 1e-12:
        raise ValueError('Invalid disk radius or normal.')
    normal = normal/norm
    return radius_m**2/4*(np.eye(3)+np.outer(normal,normal))


@dataclass(frozen=True)
class ExhaustPort:
    stage_index: int                 # 0=stage 1, 1=stage 2
    station_body_m: np.ndarray       # fixed vehicle-base datum, not a CG arm
    mass_flow_kgps: float            # q >= 0 outflow
    mass_flow_rate_kgps2: float = 0.0 # dq/dt; total derivative along the trajectory
    gyration_body_m2: np.ndarray = field(default_factory=lambda: np.zeros((3,3)))

    def __post_init__(self):
        if self.stage_index not in (0,1):
            raise ValueError('stage_index must be 0 or 1.')
        if not np.all(np.isfinite([self.mass_flow_kgps,self.mass_flow_rate_kgps2])) or self.mass_flow_kgps < 0:
            raise ValueError('Exhaust flow must be nonnegative and its derivative finite.')
        station = vector3(self.station_body_m).copy()
        tensor = np.asarray(self.gyration_body_m2,dtype=float)
        if (tensor.shape != (3,3) or not np.all(np.isfinite(tensor))
                or not np.allclose(tensor,tensor.T,atol=1e-12)
                or np.min(np.linalg.eigvalsh(tensor)) < -1e-12):
            raise ValueError('Exit gyration must be a finite symmetric positive-semidefinite tensor.')
        station.setflags(write=False)
        tensor = tensor.copy(); tensor.setflags(write=False)
        object.__setattr__(self,'station_body_m',station)
        object.__setattr__(self,'gyration_body_m2',tensor)


def validate_ports(ports, consumption_kgps):
    """Check mass balance and the scope of the balanced axial-flow equations."""
    rates = np.zeros(2)
    derivatives = np.zeros(2)
    first = np.zeros(3)
    first_dot = np.zeros(3)
    for port in ports:
        if not isinstance(port,ExhaustPort):
            raise TypeError('Every exhaust entry must be an ExhaustPort.')
        rates[port.stage_index] += port.mass_flow_kgps
        derivatives[port.stage_index] += port.mass_flow_rate_kgps2
        first += port.mass_flow_kgps*port.station_body_m
        first_dot += port.mass_flow_rate_kgps2*port.station_body_m
    if not np.allclose(rates,consumption_kgps,rtol=1e-10,atol=1e-9):
        raise ValueError('Exhaust ports must account for every stage consumption rate.')
    scale = 1+sum(p.mass_flow_kgps*np.linalg.norm(p.station_body_m) for p in ports)
    scale_dot = 1+sum(abs(p.mass_flow_rate_kgps2)*np.linalg.norm(p.station_body_m) for p in ports)
    if np.linalg.norm(first[1:]) > 1e-10*scale or np.linalg.norm(first_dot[1:]) > 1e-10*scale_dot:
        raise ValueError('This model requires an axially balanced exhaust centroid and derivative.')
    return derivatives
