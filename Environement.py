"""Everything the vehicle interacts with but does not carry: atmosphere,
gravity, wind and the spherical ground. Dynamics asks this module; the module
never asks back. Evaluators are pure; queries outside the modeled envelope
raise, never extrapolate.

Frame convention (project-wide, see docs/env_issues.md): the ECEF origin is
the Earth centre (0,0,0), the launch site sits at (0,0,earth_radius) on +z,
and Earth rotation is not modeled: ECEF axes are treated as inertial. Shared
vectors are current Body components, Body +x toward the nose; quaternion is
Hamilton [w,x,y,z] mapping Body to ECEF.
"""
from dataclasses import dataclass
import numpy as np
from frames import body_to_ecef_matrix, vector3

# USSA 1976 / ISO 2533 constants; layers: (base altitude [m], base temperature [K], lapse rate [K/m]).
SEA_LEVEL_TEMPERATURE_K=288.15
SEA_LEVEL_PRESSURE_PA=101325.0
DRY_AIR_GAS_CONSTANT_JPKGK=287.05287
AIR_HEAT_CAPACITY_RATIO=1.4
STANDARD_GRAVITY_MPS2=9.80665
EARTH_MEAN_RADIUS_M=6371000.0
STANDARD_LAYERS=(
    (0.0,SEA_LEVEL_TEMPERATURE_K,-0.0065),
    (11000.0,216.65,0.0),
    (20000.0,216.65,0.001),
    (32000.0,228.65,0.0028),
    (47000.0,270.65,0.0),
    (51000.0,270.65,-0.0028),
    (71000.0,214.65,-0.002),
)
MESOPAUSE_ALTITUDE_M=84852.0
# Above the mesopause an isothermal tail decays to numerical vacuum: a
# documented near-vacuum placeholder to the 2000 km ceiling, not
# thermospheric physics (rationale: docs/env_issues.md).
ATMOSPHERE_CEILING_M=2000000.0
# Wind speed fraction of the surface reference wind: (altitude [m], scale).
# Stylized, documented profile: boundary-layer growth toward the ~11 km
# jet-stream core, decay through the stratosphere, exactly zero at the
# mesopause - there is no wind above the physically modeled atmosphere.
WIND_PROFILE=(
    (0.0,1.0),
    (100.0,1.2),
    (1000.0,1.45),
    (5000.0,1.75),
    (11000.0,2.0),
    (20000.0,1.1),
    (32000.0,0.5),
    (47000.0,0.15),
    (71000.0,0.02),
    (MESOPAUSE_ALTITUDE_M,0.0),
)


@dataclass(frozen=True)
class Atmosphere:
    """Air state at one point."""
    temperature_k: float
    pressure_pa: float
    density_kgm3: float
    sound_speed_mps: float


def _layer_pressure(base_pressure,base_temperature,lapse,height_above_base):
    """Hydrostatic pressure within one standard-atmosphere layer."""
    if lapse==0.0:
        return base_pressure*np.exp(-STANDARD_GRAVITY_MPS2*height_above_base
                                    /(DRY_AIR_GAS_CONSTANT_JPKGK*base_temperature))
    temperature=base_temperature+lapse*height_above_base
    return base_pressure*(temperature/base_temperature)**(-STANDARD_GRAVITY_MPS2
                                                          /(DRY_AIR_GAS_CONSTANT_JPKGK*lapse))


def _wind_scale(altitude_m):
    """Piecewise-linear WIND_PROFILE interpolation; zero above the mesopause."""
    knots=WIND_PROFILE
    height=max(0.0,float(altitude_m))
    if height >= knots[-1][0]:
        return float(knots[-1][1])
    for index in range(len(knots)-1):
        base,scale=knots[index]
        top,upper=knots[index+1]
        if height < top:
            fraction=(height-base)/(top-base)
            return float(scale+fraction*(upper-scale))
    return float(knots[-1][1])


class Environment:
    """Scenario configuration plus pure evaluators; validated once at construction."""

    def __init__(self,*,surface_gravity_mps2=STANDARD_GRAVITY_MPS2,
                 earth_radius_m=EARTH_MEAN_RADIUS_M,ground_altitude_m=0.0,
                 wind_ecef_mps=(0,0,0),wind_rotation_rad=0.0):
        gravity=float(surface_gravity_mps2)
        wind=vector3(wind_ecef_mps)
        wind_rotation=float(wind_rotation_rad)
        radius=float(earth_radius_m)
        if (not np.all(np.isfinite([gravity,ground_altitude_m,radius,wind_rotation]))
                or gravity < 1e-12 or not np.all(np.isfinite(wind))
                or radius <= 0):
            raise ValueError('Invalid surface gravity, Earth radius, ground altitude, '
                             'wind or wind rotation.')
        self.planet_center_ecef_m=np.zeros(3)
        self.launch_site_ecef_m=np.array([0.0,0.0,radius])
        self.earth_radius_m=radius
        # mu = g_surface * R^2 reproduces the configured gravity at the launch site.
        self.gravity_mu_m3ps2=gravity*radius*radius
        self.ground_altitude_m=float(ground_altitude_m)
        # The scenario wind is the surface reference wind rotated about the
        # launch-site vertical (ECEF +z), positive from ECEF +x toward +y: the
        # rotation spreads the wind over the horizontal plane so crossflow
        # (and not one lone coordinate) angles the vehicle.
        cosine,sine=np.cos(wind_rotation),np.sin(wind_rotation)
        self.wind_ecef_mps=np.array([cosine*wind[0]-sine*wind[1],
                                     sine*wind[0]+cosine*wind[1],wind[2]])
        self.wind_rotation_rad=wind_rotation

    @staticmethod
    def atmosphere(altitude_m):
        """USSA 1976 layers to the mesopause plus an isothermal tail to the
        ceiling; altitude 0..2000000 m, rejected outside."""
        h=float(altitude_m)
        if not np.isfinite(h) or h < -1e-6 or h > ATMOSPHERE_CEILING_M:
            raise ValueError('The environment model supports altitude 0..2000000 m.')
        h=max(0.0,h)
        index=len(STANDARD_LAYERS)-1
        for candidate in range(len(STANDARD_LAYERS)-1):
            if h < STANDARD_LAYERS[candidate+1][0]:
                index=candidate
                break
        pressure=SEA_LEVEL_PRESSURE_PA
        for layer in range(index):
            base,base_temperature,lapse=STANDARD_LAYERS[layer]
            next_base=STANDARD_LAYERS[layer+1][0]
            pressure=_layer_pressure(pressure,base_temperature,lapse,next_base-base)
        base,base_temperature,lapse=STANDARD_LAYERS[index]
        if h <= MESOPAUSE_ALTITUDE_M:
            temperature=base_temperature+lapse*(h-base)
            pressure=_layer_pressure(pressure,base_temperature,lapse,h-base)
        else:
            mesopause_temperature=base_temperature+lapse*(MESOPAUSE_ALTITUDE_M-base)
            mesopause_pressure=_layer_pressure(pressure,base_temperature,lapse,MESOPAUSE_ALTITUDE_M-base)
            temperature=mesopause_temperature
            pressure=mesopause_pressure*np.exp(-STANDARD_GRAVITY_MPS2*(h-MESOPAUSE_ALTITUDE_M)
                                                /(DRY_AIR_GAS_CONSTANT_JPKGK*mesopause_temperature))
        density=pressure/(DRY_AIR_GAS_CONSTANT_JPKGK*temperature)
        sound_speed=np.sqrt(AIR_HEAT_CAPACITY_RATIO*DRY_AIR_GAS_CONSTANT_JPKGK*temperature)
        return Atmosphere(float(temperature),float(pressure),float(density),float(sound_speed))

    def altitude(self,state,properties,station_body_m=None):
        """Height of the CG, or of a station measured from the CG, above the spherical ground."""
        position=state.position_body_m
        if station_body_m is not None:
            position=position+vector3(station_body_m)-properties.cg_station_body_m
        return self._altitude_above_ground(self._rotation_to_ecef(state)@position)

    def wind_body(self,state):
        """The altitude-scaled ECEF wind expressed in current Body components.

        The scenario wind is the surface reference wind rotated about the
        launch-site vertical, multiplied by the WIND_PROFILE fraction of the
        CG altitude: zero above the mesopause, the edge of the physically
        modeled atmosphere.
        """
        rotation=self._rotation_to_ecef(state)
        scale=_wind_scale(self._altitude_above_ground(rotation@state.position_body_m))
        if scale==0.0:
            return np.zeros(3)
        return rotation.T@(scale*self.wind_ecef_mps)

    def _altitude_above_ground(self,position_ecef):
        return float(np.linalg.norm(position_ecef)-self.earth_radius_m+self.ground_altitude_m)

    def gravity_body(self,state):
        """Inverse-square point-mass gravity toward the planet centre, in Body components."""
        rotation=self._rotation_to_ecef(state)
        position_ecef=rotation@state.position_body_m
        distance=np.linalg.norm(position_ecef)
        if distance < 1e-6:
            raise ValueError('Gravity is undefined at the planet centre.')
        magnitude=self.gravity_mu_m3ps2/distance**2
        return rotation.T@(-magnitude*position_ecef/distance)

    @staticmethod
    def _rotation_to_ecef(state):
        return body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
