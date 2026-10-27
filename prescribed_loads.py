"""Explicit verification loads, independent of engine or aerodynamic models."""
from dataclasses import dataclass
import numpy as np
from frames import vector3
from load_interface import AppliedForce


@dataclass(frozen=True)
class SmoothForceDoublet:
    """Fixed Body station/direction, positive then negative sin-cubed pulse.

    This is an A-class verification input, not a measured actuator/wind law.
    Force and its first two derivatives are continuous at the endpoints.
    """
    start_s: float
    duration_s: float
    peak_force_body_n: np.ndarray
    station_body_m: np.ndarray

    def __post_init__(self):
        if not np.all(np.isfinite([self.start_s,self.duration_s])) or self.duration_s<=0:
            raise ValueError('A force doublet needs finite start and positive duration.')
        for name in ('peak_force_body_n','station_body_m'):
            object.__setattr__(self,name,vector3(getattr(self,name)).copy())

    @property
    def breakpoints_s(self):
        return (self.start_s,self.start_s+self.duration_s/2,self.start_s+self.duration_s)

    def __call__(self,time_s):
        if not np.isfinite(time_s):
            raise ValueError('Load evaluation time must be finite.')
        fraction=(time_s-self.start_s)/self.duration_s
        scale=np.sin(2*np.pi*fraction)**3 if 0<fraction<1 else 0.0
        return AppliedForce(scale*self.peak_force_body_n,self.station_body_m)


def force_doublet_from_config(record):
    if record['kind']!='smooth_force_doublet':
        raise ValueError('Unknown prescribed-load kind.')
    return SmoothForceDoublet(record['start_s'],record['duration_s'],
                              record['peak_force_body_n'],record['station_body_m'])
