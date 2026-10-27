"""A/B boundary: structural force stations -> resultant force and moment at CG.

Body +x is toward the nose. All input points use the same fixed structural
datum as mass properties. This is a load-assembly utility, not B's engine law.
"""
from dataclasses import dataclass
import numpy as np
from frames import vector3
from dynamics import Loads


@dataclass(frozen=True)
class AppliedForce:
    force_body_n: np.ndarray
    station_body_m: np.ndarray

    def __post_init__(self):
        for field in ('force_body_n', 'station_body_m'):
            object.__setattr__(self, field, vector3(getattr(self, field)).copy())


def assemble_loads(properties, applied_forces, *, gravity_body_mps2,
                   free_moment_body_nm=(0, 0, 0), consumption_kgps=(0, 0),exhaust_ports=()):
    """Sum non-gravitational forces; add each arm cross force exactly once.

    Pure couples belong in free_moment_body_nm. Gravity stays separate and is
    applied by A. B must NOT supply an additional r cross F for a force already
    represented here. Current properties must be recalculated at RK4 trial states.
    """
    force = np.zeros(3)
    moment = vector3(free_moment_body_nm).copy()
    for applied in applied_forces:
        force += applied.force_body_n
        moment += np.cross(properties.arm_from_station(applied.station_body_m),
                           applied.force_body_n)
    return Loads(force, moment, gravity_body_mps2, consumption_kgps,tuple(exhaust_ports))
