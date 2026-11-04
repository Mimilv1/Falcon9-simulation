"""Environment model covering all simulated flight phases of this project.

Everything the vehicle interacts with but does not carry: atmosphere, gravity,
wind and the ground reference. Dynamics code asks this module; the module never
asks back.

Envelope and assumptions, each replaceable without touching callers:
- Atmosphere: dry US Standard Atmosphere 1976 layers, 0..84852 m, extended by an
  explicit isothermal tail to the 2000 km model ceiling. Local geometric height is
  used as geopotential height (A); the deviation grows with altitude and is
  documented in first_ascent_segment.md. Above ~100 km the tail is a documented
  near-vacuum placeholder, not thermospheric physics. Queries outside
  0..2000000 m are rejected, never extrapolated.
- Gravity: a vector. The surface gravity magnitude at the launch site plus
  the planet radius define a spherical planet whose centre is the ECEF origin
  (0,0,0); the launch site sits at (0,0,earth_radius) on the +z axis. Earth is
  assimilated to a point of constant mass at that centre: gravity is mu/d^2
  along the normalized vector from the centre to the vehicle, so its direction
  follows the vehicle around the planet.
- Wind: one fixed vector expressed in ECEF; air above ~50 km
  is too thin for wind to matter at this fidelity.
- Ground: a sphere concentric with the planet; the launch site sits
  ground_altitude_m above it.

Shared position, velocity and force vectors are current Body components; Body
+x points toward the nose. Quaternion is Hamilton [w,x,y,z] mapping current
Body to ECEF. Earth rotation is not modeled in v0: the ECEF axes are treated
as non-rotating (inertial).
"""
from dataclasses import dataclass
import numpy as np
from frames import body_to_ecef_matrix, vector3

# Standard-atmosphere constants (USSA 1976 / ISO 2533).
SEA_LEVEL_TEMPERATURE_K=288.15
SEA_LEVEL_PRESSURE_PA=101325.0
DRY_AIR_GAS_CONSTANT_JPKGK=287.05287
AIR_HEAT_CAPACITY_RATIO=1.4
STANDARD_GRAVITY_MPS2=9.80665
EARTH_MEAN_RADIUS_M=6371000.0

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
# The 2000 km ceiling is the mission envelope, and two properties make it safe:
# - Below ~100 km the table is USSA 1976 evidence; above the mesopause the
#   isothermal tail decays to numerical vacuum (density < 1e-150 by 2000 km),
#   so drag and ambient-pressure thrust read as exactly negligible there
#   instead of wrong-but-plausible values. The tail is a placeholder, not
#   thermospheric physics: thermosphere temperature and composition above
#   ~100 km are not modeled.
# - The former 120 km rejection guarded the flat-Earth frame, in which high
#   arcs produced plausible-but-wrong trajectories. Gravity and altitude are
#   now spherical (see gravity_body), so suborbital arcs to 2000 km are
#   represented honestly. The ECEF datum is adopted project-wide (origin at
#   the Earth centre, launch site at (0,0,earth_radius)); Earth rotation and
#   orbital phases remain out of scope. Per team convention, unmodeled
#   regimes are rejected.
ATMOSPHERE_CEILING_M=2000000.0


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

    def __init__(self,*,surface_gravity_mps2=STANDARD_GRAVITY_MPS2,
                 earth_radius_m=EARTH_MEAN_RADIUS_M,ground_altitude_m=0.0,
                 wind_ecef_mps=(0,0,0)):
        gravity=float(surface_gravity_mps2)
        wind=vector3(wind_ecef_mps)
        radius=float(earth_radius_m)
        if (not np.all(np.isfinite([ground_altitude_m,radius]))
                or gravity < 1e-12 or not np.all(np.isfinite(wind))
                or radius <= 0):
            raise ValueError('Invalid surface gravity, Earth radius, ground altitude or wind.')
        # ECEF convention: the planet centre is the origin (0,0,0) and the
        # launch site sits at (0,0,earth_radius); surface gravity there points
        # toward the centre, along -z.
        self.planet_center_ecef_m=np.zeros(3)
        self.launch_site_ecef_m=np.array([0.0,0.0,radius])
        self.surface_gravity_ecef_mps2=np.array([0.0,0.0,-gravity])
        self.earth_radius_m=radius
        # mu = g_surface * R^2 reproduces the configured magnitude at the launch site.
        self.gravity_mu_m3ps2=gravity*radius*radius
        self.ground_altitude_m=float(ground_altitude_m)
        self.wind_ecef_mps=wind.copy()

    @staticmethod
    def atmosphere(altitude_m):
        """USSA 1976 layered standard atmosphere plus an explicit isothermal tail.

        Layers 0..84852 m reproduce the published table (evidence: USSA 1976 /
        ISO 2533, geopotential height). Above the mesopause the air is modeled
        as isothermal at the mesopause temperature with hydrostatic exponential
        decay (B); at those densities the choice only refines where drag is
        already negligible. Altitude outside 0..2000000 m is unmodeled and
        rejected; the clamp only absorbs sub-micrometre roundoff at the ground.
        """
        h=float(altitude_m)
        if not np.isfinite(h) or h < -1e-6 or h > ATMOSPHERE_CEILING_M:
            raise ValueError('The environment model supports altitude 0..2000000 m.')
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
        """Height of the CG, or of an optional structural station, above the spherical ground."""
        position=self._station_position_body(state,properties,station_body_m)
        return float(self._height_above_ground(state,position))

    def gravity_body(self,state):
        """Gravity in current Body components; Newtonian point-mass attraction.

        Earth is assimilated to a point of constant mass at the ECEF origin
        (no gravity field computation): gravity is mu/d^2 along the normalized
        vector from the centre to the vehicle, where d is the distance between
        the two and mu = g_surface * R^2 is the point-mass parameter fixed
        by the configured surface gravity magnitude at the launch site
        (0,0,earth_radius).
        """
        position_ecef=self._position_ecef(state,state.position_body_m)
        distance=np.linalg.norm(position_ecef)
        if distance < 1e-6:
            raise ValueError('Gravity is undefined at the planet centre.')
        # Normalized vector from the planet centre to the vehicle.
        direction=position_ecef/distance
        magnitude=self.gravity_mu_m3ps2/distance**2
        return self._rotation_to_ecef(state).T@(-magnitude*direction)

    def wind_body(self,state):
        """Fixed ECEF wind expressed in current Body components."""
        return self._rotation_to_ecef(state).T@self.wind_ecef_mps

    @staticmethod
    def _layer_pressure(base_pressure,base_temperature,lapse,height_above_base):
        """Hydrostatic pressure within one standard-atmosphere layer."""
        if lapse==0.0:
            return base_pressure*np.exp(-STANDARD_GRAVITY_MPS2*height_above_base
                                        /(DRY_AIR_GAS_CONSTANT_JPKGK*base_temperature))
        temperature=base_temperature+lapse*height_above_base
        return base_pressure*(temperature/base_temperature)**(-STANDARD_GRAVITY_MPS2
                                                              /(DRY_AIR_GAS_CONSTANT_JPKGK*lapse))

    def _height_above_ground(self,state,position_body_m):
        """Radial distance from the planet centre minus the ground-sphere radius."""
        position_ecef=self._position_ecef(state,position_body_m)
        return np.linalg.norm(position_ecef)-self.earth_radius_m+self.ground_altitude_m

    def _position_ecef(self,state,position_body_m):
        """Position vector in ECEF; the planet centre is the origin."""
        return self._rotation_to_ecef(state)@position_body_m

    @staticmethod
    def _rotation_to_ecef(state):
        return body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)

    @staticmethod
    def _station_position_body(state,properties,station_body_m):
        """Position of the CG, or of an optional station measured from the CG."""
        position=state.position_body_m
        if station_body_m is not None:
            position=position+vector3(station_body_m)-properties.cg_station_body_m
        return position
