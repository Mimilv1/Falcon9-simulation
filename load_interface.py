from dataclasses import dataclass
import numpy as np
from frames import vector3
from dynamics import Loads


@dataclass(frozen=True)
class AppliedForce:
    force_body_n: np.ndarray
    station_body_m: np.ndarray

    def __post_init__(self):
        for field in ("force_body_n", "station_body_m"):
            object.__setattr__(self, field, vector3(getattr(self, field)).copy())


def assemble_loads(
    properties,
    applied_forces,
    *,
    gravity_body_mps2,
    free_moment_body_nm=(0, 0, 0),
    consumption_kgps=(0, 0),
    exhaust_ports=(),
    frame_acceleration_body_mps2=(0, 0, 0)
):
    force = np.zeros(3)
    moment = vector3(free_moment_body_nm).copy()
    for applied in applied_forces:
        force += applied.force_body_n
        moment += np.cross(
            properties.arm_from_station(applied.station_body_m), applied.force_body_n
        )
    return Loads(
        force,
        moment,
        gravity_body_mps2,
        consumption_kgps,
        tuple(exhaust_ports),
        frame_acceleration_body_mps2,
    )
