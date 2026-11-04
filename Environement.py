"""Environment model covering all simulated flight phases of this project.

Everything the vehicle interacts with but does not carry: atmosphere, gravity,
wind and the ground reference. Dynamics code asks this module; the module never
asks back.

Envelope and assumptions, each replaceable without touching callers:
- Atmosphere: dry US Standard Atmosphere 1976 layers, 0..84852 m, extended by an
  explicit isothermal tail, 84852..120000 m. Local geometric height is used as
  geopotential height (A); the deviation grows with altitude and is documented
  in first_ascent_segment.md. Queries outside 0..120000 m are rejected, never
  extrapolated.
- Gravity: direction fixed along the initial up axis (flat-Earth convention).
  Magnitude is constant by default, or inverse-square from the configured
  surface value when an Earth radius is supplied.
- Wind: one fixed vector expressed in the initial Body basis; air above ~50 km
  is too thin for wind to matter at this fidelity.
- Ground: a single plane of constant altitude along the initial up axis.

Shared position, velocity and force vectors are current Body components; Body
+x points toward the nose. Quaternion is Hamilton [w,x,y,z] mapping current
Body to the fixed initial Body basis.
"""
from dataclasses import dataclass
import numpy as np
from frames import body_to_initial_body_matrix, vector3

# Standard-atmosphere constants (USSA 1976 / ISO 2533).
SEA_LEVEL_TEMPERATURE_K=288.15
SEA_LEVEL_PRESSURE_PA=101325.0
DRY_AIR_GAS_CONSTANT_JPKGK=287.05287
AIR_HEAT_CAPACITY_RATIO=1.4
STANDARD_GRAVITY_MPS2=9.80665

# USSA 1976 layers: (base altitude [m], base temperature [K], lapse rate [K/m]).
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
# The 120 km ceiling is a modeling limit, not an arbitrary knob:
# - The isothermal tail assumes constant temperature (187 K) and uniform air
#   composition; above ~100 km the thermosphere heats to 600..1000+ K and the
#   constituents separate by mass, so R and gamma are no longer valid there.
# - Between 120 km and any higher ceiling, drag and ambient-pressure thrust
#   are already negligible, so a taller table would have no consumers.
# - Beyond ~120 km trajectories become orbital, which this simulation's
#   flat-Earth frame cannot represent (fixed gravity direction, no Earth
#   rotation, altitude above a ground plane). Answering there would produce
#   plausible numbers from an invalid model; per team convention, unmodeled
#   regimes are rejected instead. Orbital phases require a new frame model
#   and environment class, decided at team level.
ATMOSPHERE_CEILING_M=120000.0


@dataclass(frozen=True)
class Atmosphere:
    """Air state at one point; the only multi-number output of this module."""
    temperature_k: float
    pressure_pa: float
    density_kgm3: float
    sound_speed_mps: float


class TroposphereEnvironment:
    """Scenario configuration plus evaluators; constructed once per run.

    All validation happens here at construction, so an invalid scenario fails
    before flight, not mid-integration. Evaluators stay pure: they read state,
    return values, and keep no memory between calls. The class name predates
    the layered extension and is kept for import stability.
    """

    def __init__(self,*,gravity_mps2=STANDARD_GRAVITY_MPS2,gravity_earth_radius_m=None,
                 ground_altitude_m=0.0,up_initial_body=(1,0,0),
                 wind_initial_body_mps=(0,0,0)):
        up=vector3(up_initial_body)
        wind=vector3(wind_initial_body_mps)
        radius=None if gravity_earth_radius_m is None else float(gravity_earth_radius_m)
        if (not np.all(np.isfinite([gravity_mps2,ground_altitude_m]))
                or gravity_mps2 < 0 or np.linalg.norm(up) < 1e-12
                or not np.all(np.isfinite(wind))
                or (radius is not None and (not np.isfinite(radius) or radius <= 0))):
            raise ValueError('Invalid gravity, ground altitude or initial vertical.')
        self.gravity_mps2=float(gravity_mps2)
        self.gravity_earth_radius_m=radius
        self.ground_altitude_m=float(ground_altitude_m)
        self.up_initial_body=up/np.linalg.norm(up)
        self.wind_initial_body_mps=wind.copy()

    @staticmethod
    def atmosphere(altitude_m):
        """USSA 1976 layered standard atmosphere plus an explicit isothermal tail.

        Layers 0..84852 m reproduce the published table (evidence: USSA 1976 /
        ISO 2533, geopotential height). Above the mesopause the air is modeled
        as isothermal at the mesopause temperature with hydrostatic exponential
        decay (B); at those densities the choice only refines where drag is
        already negligible. Altitude outside 0..120000 m is unmodeled and
        rejected; the clamp only absorbs sub-micrometre roundoff at the ground.
        """
        h=float(altitude_m)
        if not np.isfinite(h) or h < -1e-6 or h > ATMOSPHERE_CEILING_M:
            raise ValueError('The environment model supports altitude 0..120000 m.')
        h=max(0.0,h)
        # Layer containing h; pressure bases are walked from sea level.
        index=len(STANDARD_LAYERS)-1
        for candidate in range(len(STANDARD_LAYERS)-1):
            if h < STANDARD_LAYERS[candidate+1][0]:
                index=candidate
                break
        pressure=SEA_LEVEL_PRESSURE_PA
        for layer in range(index):
            base,base_temperature,lapse=STANDARD_LAYERS[layer]
            next_base=STANDARD_LAYERS[layer+1][0]
            pressure=TroposphereEnvironment._layer_pressure(pressure,base_temperature,lapse,next_base-base)
        base,base_temperature,lapse=STANDARD_LAYERS[index]
        if h <= MESOPAUSE_ALTITUDE_M:
            temperature=base_temperature+lapse*(h-base)
            pressure=TroposphereEnvironment._layer_pressure(pressure,base_temperature,lapse,h-base)
        else:
            # Isothermal tail continuing the top tabulated layer.
            mesopause_temperature=base_temperature+lapse*(MESOPAUSE_ALTITUDE_M-base)
            mesopause_pressure=TroposphereEnvironment._layer_pressure(pressure,base_temperature,lapse,MESOPAUSE_ALTITUDE_M-base)
            temperature=mesopause_temperature
            pressure=mesopause_pressure*np.exp(-STANDARD_GRAVITY_MPS2*(h-MESOPAUSE_ALTITUDE_M)
                                                /(DRY_AIR_GAS_CONSTANT_JPKGK*mesopause_temperature))
        density=pressure/(DRY_AIR_GAS_CONSTANT_JPKGK*temperature)
        sound_speed=np.sqrt(AIR_HEAT_CAPACITY_RATIO*DRY_AIR_GAS_CONSTANT_JPKGK*temperature)
        return Atmosphere(float(temperature),float(pressure),float(density),float(sound_speed))

    def altitude(self,state,properties,station_body_m=None):
        """Height of the CG, or of an optional structural station, above the ground plane."""
        position=self._station_position_body(state,properties,station_body_m)
        return float(self._height_along_up(state,position))

    def gravity_body(self,state):
        """Gravity in current Body components; fixed direction along the up axis.

        Magnitude is the configured surface value, or that value scaled
        inverse-square with the CG height along the up axis when an Earth
        radius is configured.
        """
        rotation=self._rotation_to_initial_body(state)
        height=self.ground_altitude_m+self.up_initial_body@(rotation@state.position_body_m)
        magnitude=self.gravity_mps2
        if self.gravity_earth_radius_m is not None:
            magnitude*= (self.gravity_earth_radius_m/(self.gravity_earth_radius_m+height))**2
        return -magnitude*(rotation.T@self.up_initial_body)

    def wind_body(self,state):
        """Fixed initial-basis wind expressed in current Body components."""
        rotation=self._rotation_to_initial_body(state)
        return rotation.T@self.wind_initial_body_mps

    @staticmethod
    def _layer_pressure(base_pressure,base_temperature,lapse,height_above_base):
        """Hydrostatic pressure within one standard-atmosphere layer."""
        if lapse==0.0:
            return base_pressure*np.exp(-STANDARD_GRAVITY_MPS2*height_above_base
                                        /(DRY_AIR_GAS_CONSTANT_JPKGK*base_temperature))
        temperature=base_temperature+lapse*height_above_base
        return base_pressure*(temperature/base_temperature)**(-STANDARD_GRAVITY_MPS2
                                                              /(DRY_AIR_GAS_CONSTANT_JPKGK*lapse))

    def _height_along_up(self,state,position_body_m):
        rotation=self._rotation_to_initial_body(state)
        return self.ground_altitude_m+self.up_initial_body@(rotation@position_body_m)

    @staticmethod
    def _rotation_to_initial_body(state):
        return body_to_initial_body_matrix(state.quaternion_body_to_initial_body_wxyz)

    @staticmethod
    def _station_position_body(state,properties,station_body_m):
        """Position of the CG, or of an optional station measured from the CG."""
        position=state.position_body_m
        if station_body_m is not None:
            position=position+vector3(station_body_m)-properties.cg_station_body_m
        return position
