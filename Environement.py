"""First B environment model for low-altitude powered ascent.

Dry standard troposphere (0..11 km), with local geometric height used as
geopotential height (A). Constant gravity magnitude and initial-basis wind are
scenario assumptions. Shared vectors remain current Body components.
"""
from dataclasses import dataclass
import numpy as np
from frames import body_to_initial_body_matrix, vector3


@dataclass(frozen=True)
class Atmosphere:
    temperature_k: float
    pressure_pa: float
    density_kgm3: float
    sound_speed_mps: float


class TroposphereEnvironment:
    def __init__(self, *, gravity_mps2=9.80665, ground_altitude_m=0.0,
                 up_initial_body=(1,0,0),wind_initial_body_mps=(0,0,0)):
        self.gravity_mps2=float(gravity_mps2)
        self.ground_altitude_m=float(ground_altitude_m)
        up=vector3(up_initial_body)
        if (not np.all(np.isfinite([gravity_mps2,ground_altitude_m]))
                or gravity_mps2 < 0 or np.linalg.norm(up) < 1e-12):
            raise ValueError('Invalid gravity, ground altitude or initial vertical.')
        self.up_initial_body=up/np.linalg.norm(up)
        self.wind_initial_body_mps=vector3(wind_initial_body_mps).copy()

    @staticmethod
    def atmosphere(altitude_m):
        """Standard troposphere hydrostatic/ideal-gas model; not actual weather."""
        h=float(altitude_m)
        if not np.isfinite(h) or h < -1e-6 or h > 11000:
            raise ValueError('The first environment model supports altitude 0..11000 m.')
        h=max(0.0,h)  # sub-micrometre roundoff at the ground only
        temperature=288.15-.0065*h
        gas_constant=287.05287
        pressure=101325.0*(temperature/288.15)**(9.80665/(gas_constant*.0065))
        density=pressure/(gas_constant*temperature)
        return Atmosphere(temperature,pressure,density,np.sqrt(1.4*gas_constant*temperature))

    def altitude(self,state,properties,station_body_m=None):
        position=state.position_body_m
        if station_body_m is not None:
            position=position+vector3(station_body_m)-properties.cg_station_body_m
        rotation=body_to_initial_body_matrix(state.quaternion_body_to_initial_body_wxyz)
        return float(self.ground_altitude_m+self.up_initial_body@(rotation@position))

    def gravity_body(self,state):
        rotation=body_to_initial_body_matrix(state.quaternion_body_to_initial_body_wxyz)
        return -self.gravity_mps2*(rotation.T@self.up_initial_body)

    def wind_body(self,state):
        rotation=body_to_initial_body_matrix(state.quaternion_body_to_initial_body_wxyz)
        return rotation.T@self.wind_initial_body_mps
