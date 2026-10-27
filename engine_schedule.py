"""Smooth A-class shutdown verification envelope; not measured engine response."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SmoothShutdown:
    start_s: float
    duration_s: float

    def __post_init__(self):
        if not np.all(np.isfinite([self.start_s,self.duration_s])) or self.start_s<0 or self.duration_s<=0:
            raise ValueError('Require nonnegative shutdown start and positive duration.')

    @property
    def end_s(self):return self.start_s+self.duration_s

    @property
    def breakpoints_s(self):return (self.start_s,self.end_s)

    def __call__(self,time_s):
        """Return amplitude and its analytic derivative in 1/s."""
        if not np.isfinite(time_s):raise ValueError('Time must be finite.')
        if time_s<=self.start_s:return 1.0,0.0
        if time_s>=self.end_s:return 0.0,0.0
        u=(time_s-self.start_s)/self.duration_s
        # Quintic smoothstep: amplitude, first/second derivative join at both ends.
        amplitude=(1-u)**3*(1+3*u+6*u*u)
        derivative=-30*u*u*(1-u)**2/self.duration_s
        return float(amplitude),float(derivative)
