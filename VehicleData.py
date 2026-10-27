"""Team objects: Body +x nose; all shared vectors use current Body axes."""
from copy import deepcopy
from dataclasses import dataclass, field, replace
import numpy as np
from frames import body_to_initial_body_matrix, normalize_quaternion, vector3
from mass_properties import combine_components, fixed_cylinder
from vehicle_parameters import Falcon9Parameters
from vehicle_mass_model import VehicleModel, Configuration


def quat_rotate(q, v):
    """Rotate current Body components into the fixed initial Body basis."""
    return body_to_initial_body_matrix(q) @ vector3(v)


@dataclass(frozen=True)
class EngineSpec:
    thrust_vac: float
    isp_vac: float
    isp_sl: float | None = None
    exit_area: float | None = None
    throttle_min: float | None = None
    gimbal_max: float | None = None
    dry_mass: float | None = None  # already included in stage dry mass


@dataclass
class Engine:
    spec: EngineSpec
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))  # stage-base Body coordinates
    running: bool = False
    throttle: float = 0.0
    gimbal_pitch: float = 0.0
    gimbal_yaw: float = 0.0
    ignition_count: int = 0


@dataclass(frozen=True)
class TankSpec:
    name: str
    initial_mass_kg: float
    density_kgm3: float
    reference_temperature_k: float
    inner_radius_m: float
    bottom_axial_m: float
    section_length_m: float


@dataclass
class Stage:
    dry_mass: float
    propellant_capacity: float
    propellant_mass: float
    length: float
    diameter: float
    cg_dry: float | None = None
    inertia_dry: np.ndarray | None = None  # own dry CG, Body x/y/z
    engines: list = field(default_factory=list)
    offset: float = 0.0  # axial stage-base station relative to Vehicle base
    name: str = 'stage'
    tank_specs: tuple[TankSpec, ...] = ()
    dry_components: tuple = ()
    liquid_profile: object | None = field(default=None, repr=False)

    def mass_properties(self, propellant_mass=None, consumption_kgps=0.0,consumption_rate_kgps2=0.0):
        remaining = self.propellant_mass if propellant_mass is None else float(propellant_mass)
        if (not np.all(np.isfinite([remaining, consumption_kgps,consumption_rate_kgps2]))
                or not 0 <= remaining <= self.propellant_capacity+1e-7 or consumption_kgps < 0):
            raise ValueError('Invalid remaining propellant or consumption rate.')
        if self.dry_components:
            if not np.isclose(sum(c.mass_kg for c in self.dry_components),self.dry_mass,atol=1e-7,rtol=0):
                raise ValueError('Rebuild dry components after changing the dry-mass budget.')
            components = list(self.dry_components)
            dry = combine_components(components)
            if (not np.isclose(dry.cg_station_body_m[0],self.cg_dry)
                    or not np.allclose(dry.inertia_body_kgm2,self.inertia_dry)):
                raise ValueError('Dry CG/inertia snapshots must agree with the component budget.')
        else:
            component = fixed_cylinder(self.name+'_dry',self.dry_mass,0,self.cg_dry,self.length,self.diameter/2)
            if self.inertia_dry is not None:
                inertia = np.asarray(self.inertia_dry,dtype=float)
                if inertia.shape != (3,3) or not np.all(np.isfinite(inertia)) or not np.allclose(inertia,inertia.T):
                    raise ValueError('Dry inertia must be a finite symmetric 3x3 matrix.')
                component = replace(component,inertia_at_center_kgm2=inertia.copy())
            components = [component]
        if self.liquid_profile is not None:
            components.extend(self.liquid_profile.stage1_liquid_components(remaining,consumption_kgps,consumption_rate_kgps2))
        else:
            components.append(fixed_cylinder(self.name+'_propellant',remaining,-consumption_kgps,
                                              self.length/2,self.length,self.diameter/2,
                                              mass_acceleration_kgps2=-consumption_rate_kgps2))
        return combine_components(components)

    def tank_inventories(self):
        """Derive liquid inventories from the single authoritative total state."""
        if self.liquid_profile is None:
            return {}
        masses = self.liquid_profile.stage1_inventory(self.propellant_mass)
        return {spec.name:mass for spec,mass in zip(self.tank_specs,masses)}


@dataclass
class Vehicle:
    stages: list = field(default_factory=list)
    fairing_mass: float = 0.0
    payload_mass: float = 0.0
    fairing_attached: bool = True
    t: float = 0.0
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))  # initial origin -> structural base, CURRENT Body m
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))  # base velocity relative to initial origin, CURRENT Body m/s
    quaternion: np.ndarray = field(default_factory=lambda: np.array([1.,0.,0.,0.]))
    omega: np.ndarray = field(default_factory=lambda: np.zeros(3))  # Body rad/s
    parameters: object | None = field(default=None, repr=False)
    datum_shift_m: float = 0.0  # current base station relative to original datum

    def __post_init__(self):
        self.position = vector3(self.position).copy()
        self.velocity = vector3(self.velocity).copy()
        self.omega = vector3(self.omega).copy()
        self.quaternion = normalize_quaternion(self.quaternion)

    def separate(self, index=1):
        """No impulse; preserve positions/velocities of corresponding structural points."""
        if not 0 < index < len(self.stages):
            raise ValueError('index out of range')
        cut = self.stages[index].offset
        r_body = np.array([cut,0.0,0.0])
        lower,upper = deepcopy(self),deepcopy(self)
        lower.stages = lower.stages[:index]
        lower.fairing_mass = lower.payload_mass = 0.0
        lower.fairing_attached = False
        upper.stages = upper.stages[index:]
        for stage in upper.stages:
            stage.offset -= cut
        upper.datum_shift_m += cut
        upper.position = self.position+r_body
        upper.velocity = self.velocity+np.cross(self.omega,r_body)
        return lower,upper

    def getmass(self):
        fairing = self.fairing_mass if self.fairing_attached else 0.0
        return sum(s.dry_mass+s.propellant_mass for s in self.stages)+fairing+self.payload_mass

    def get_mass_properties(self, *, propellant_by_stage=None, consumption_by_stage=None,
                            consumption_derivative_by_stage=None):
        propellant_by_stage = propellant_by_stage or {}
        consumption_by_stage = consumption_by_stage or {}
        consumption_derivative_by_stage = consumption_derivative_by_stage or {}
        components = []
        for stage in self.stages:
            properties = stage.mass_properties(propellant_by_stage.get(stage.name),consumption_by_stage.get(stage.name,0.0),
                                                consumption_derivative_by_stage.get(stage.name,0.0))
            components.extend(replace(c,center_axial_m=c.center_axial_m+stage.offset) for c in properties.components)
        if (self.fairing_attached and self.fairing_mass) or self.payload_mass:
            if self.parameters is None:
                raise ValueError('Fairing/payload geometry must be supplied explicitly.')
            p = self.parameters
            if self.fairing_attached and self.fairing_mass:
                components.append(fixed_cylinder('fairing',self.fairing_mass,0,p.fairing_cg_station-self.datum_shift_m,
                                                  p.fairing_height,p.fairing_diameter/2))
            if self.payload_mass:
                components.append(fixed_cylinder('payload',self.payload_mass,0,p.payload_cg_station-self.datum_shift_m,
                                                  p.payload_length,p.payload_radius))
        return combine_components(components)

    def to_state(self, *, consumption_by_stage=None,consumption_derivative_by_stage=None):
        """Base -> geometric CG, both expressed in current Body, including migration."""
        from dynamics import State
        properties = self.get_mass_properties(consumption_by_stage=consumption_by_stage,
                                              consumption_derivative_by_stage=consumption_derivative_by_stage)
        arm = properties.cg_station_body_m
        cg_velocity_body = self.velocity+np.cross(self.omega,arm)+properties.cg_rate_body_mps
        prop = {stage.name:stage.propellant_mass for stage in self.stages}
        return State(self.position+arm,cg_velocity_body,self.quaternion,
                     self.omega,[prop.get('stage1',0),prop.get('stage2',0)])

    def apply_state(self,state,time_s, *, consumption_by_stage=None,consumption_derivative_by_stage=None):
        """Synchronize an accepted A state; do not call this for RK trial states."""
        prop = dict(zip(('stage1','stage2'),state.propellant_kg))
        properties = self.get_mass_properties(propellant_by_stage=prop,consumption_by_stage=consumption_by_stage,
                                              consumption_derivative_by_stage=consumption_derivative_by_stage)
        arm = properties.cg_station_body_m
        self.position = state.position_body_m-arm
        self.velocity = state.velocity_body_mps-np.cross(state.omega_body_radps,arm)-properties.cg_rate_body_mps
        self.quaternion = state.quaternion_body_to_initial_body_wxyz.copy()
        self.omega = state.omega_body_radps.copy()
        for stage in self.stages: stage.propellant_mass = float(prop[stage.name])
        self.t = float(time_s)

    def mass_model(self):
        return VehicleMassAdapter(self)


class VehicleMassAdapter:
    """Evaluate RK4 trial inventory without mutating shared Stage/Vehicle state."""
    def __init__(self, vehicle):
        self.vehicle = vehicle
        self.mass_model_name = ('dual_tank' if any(s.liquid_profile is not None for s in vehicle.stages)
                                else 'uniform_cylinder')

    def propellant_bounds(self,configuration=None):
        if configuration is not None:
            raise ValueError('Use configuration=None with attached Vehicle stages.')
        lower,upper = np.zeros(2),np.zeros(2)
        for stage in self.vehicle.stages:
            index = {'stage1':0,'stage2':1}[stage.name]
            upper[index] = stage.propellant_capacity
            if stage.liquid_profile is not None:
                lower[index] = stage.liquid_profile.stage1_burn_limits()['minimum_remaining_kg']
        return lower,upper

    def mass_properties(self, configuration, propellant_kg, consumption_kgps=(0,0),consumption_rate_kgps2=(0,0)):
        if configuration is not None:
            raise ValueError('Use configuration=None; attached stages are defined by this Vehicle.')
        values = np.asarray(propellant_kg,dtype=float)
        rates = np.asarray(consumption_kgps,dtype=float)
        changes = np.asarray(consumption_rate_kgps2,dtype=float)
        if (values.shape != (2,) or rates.shape != (2,) or changes.shape != (2,)
                or not np.all(np.isfinite(np.r_[values,rates,changes])) or np.any(values < 0) or np.any(rates < 0)):
            raise ValueError('Expected nonnegative finite [stage1, stage2] inventories and rates.')
        attached = {stage.name for stage in self.vehicle.stages}
        for name,mass,rate,change in zip(('stage1','stage2'),values,rates,changes):
            if name not in attached and (mass > 1e-7 or rate > 0 or change != 0):
                raise ValueError('Detached stages cannot contribute propellant or consumption.')
        prop = dict(zip(('stage1','stage2'),values))
        rate = dict(zip(('stage1','stage2'),rates))
        change = dict(zip(('stage1','stage2'),changes))
        return self.vehicle.get_mass_properties(propellant_by_stage=prop,consumption_by_stage=rate,
                                                consumption_derivative_by_stage=change)


def build_falcon9(parameters=None):
    """Create fresh team objects from one evidence profile; do not run a scenario."""
    p = parameters or Falcon9Parameters()
    profile = VehicleModel(inputs=p,mass_model='dual_tank')
    dry_components = profile.stage1_dry_components()
    dry = combine_components(dry_components)
    spec1 = EngineSpec(p.stage1_engine_thrust_vac,p.stage1_isp_vac,isp_sl=p.stage1_isp_sl,
                       exit_area=p.stage1_engine_exit_area,throttle_min=p.stage1_throttle_min,
                       gimbal_max=p.stage1_gimbal_max,dry_mass=p.stage1_engine_dry_mass)
    spec2 = EngineSpec(p.stage2_engine_thrust_vac,p.stage2_isp_vac,isp_sl=None,
                       exit_area=p.stage2_engine_exit_area,throttle_min=p.stage2_throttle_min,
                       gimbal_max=p.stage2_gimbal_max,dry_mass=p.stage2_engine_dry_mass)
    engines = [Engine(spec1,np.array([p.stage1_engine_pivot_axial,0,0]))]
    engines.extend(Engine(spec1,np.array([p.stage1_engine_pivot_axial,
                                         p.stage1_engine_ring_radius*np.cos(a),p.stage1_engine_ring_radius*np.sin(a)]))
                   for a in np.deg2rad(np.arange(0,360,45)))
    tanks = tuple(TankSpec(fluid,getattr(p,f'stage1_{fluid}_initial'),getattr(p,f'stage1_{fluid}_density'),
                            getattr(p,f'stage1_{fluid}_temperature'),getattr(p,f'stage1_{fluid}_tank_inner_radius'),
                            getattr(p,f'stage1_{fluid}_tank_bottom'),getattr(p,f'stage1_{fluid}_tank_length'))
                  for fluid in ('rp1','lox'))
    stage1 = Stage(p.stage1_dry_mass,p.stage1_propellant,p.stage1_propellant,p.stage1_length,p.body_diameter,
                   cg_dry=float(dry.cg_station_body_m[0]),inertia_dry=dry.inertia_body_kgm2.copy(),engines=engines,
                   name='stage1',tank_specs=tanks,dry_components=dry_components,liquid_profile=profile)
    dry2 = fixed_cylinder('stage2_dry',p.stage2_dry_mass,0,p.stage2_length/2,p.stage2_length,p.radius)
    stage2 = Stage(p.stage2_dry_mass,p.stage2_propellant,p.stage2_propellant,p.stage2_length,p.body_diameter,
                   cg_dry=p.stage2_length/2,inertia_dry=dry2.inertia_at_center_kgm2.copy(),
                   engines=[Engine(spec2,np.array([p.stage2_engine_pivot_axial,0,0]))],
                   offset=p.stage1_length,name='stage2')
    return Vehicle([stage1,stage2],p.fairing_mass,p.payload_mass,parameters=p,
                   quaternion=np.array(p.quaternion_body_to_initial_body_wxyz))


# Preserve the team's global entry points. The factory provides independent objects.
falcon9 = build_falcon9()
stage1,stage2 = falcon9.stages
stage1_engines = stage1.engines
merlin_1d,merlin_vac = stage1.engines[0].spec,stage2.engines[0].spec
