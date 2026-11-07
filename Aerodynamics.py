"""B aero model: axial drag, crossflow normal force, station moments, rotary
damping — one model entry per vehicle configuration.

Physical decomposition (all forces in current Body axes, +x toward the nose):
total aero force = axial force + lateral (normal) force.
  axial:   -qbar*S_ref*CA0(M)*cos(a)*|cos(a)| along Body +x. CA0 is the adopted
           axial-drag curve of the configuration's primary flow direction; the
           cos^2 scaling is the axial-dynamic-pressure share and |cos| flips the
           direction for engine-first (retro) flight.
  lateral: potential (slender-body) term qbar*S_ref*CN_alpha*sin(a)*|cos(a)|
           plus a crossflow term qbar*CDc*S_side*sin(a)^2 per body side
           segment. Both act along -(lateral relative-velocity direction).
Total angle a in [0,pi]: atan2(|v_rel_yz|, v_rel_x Norm); v=0 is a rest case.
  sin(a)=0: pure axial flight; no lateral force, no station moment.
Moments about the current CG:
  potential force at the configuration's windward-end station: the adopted
  potential_center_station_m for nose-first flow, mirrored about the CG x for
  retro flow (declared approximation), clamped to [0, reference_length].
  crossflow force of each segment at the segment side-area centroid
  (overlapping segments are incremental, e.g. grid fins on top of the body).
  rotary damping: M = c * qbar*S_ref*L_ref^2/(2*V) * omega about each axis,
  with per-configuration pitch (yaw, by symmetry) and roll coefficients.
Approximation classes follow the project evidence convention: only the
full_stack axial curve is sourced (Pezent 2024 readout, 0 deg AoA); every other
coefficient is a declared estimate, documented in data/aero_configurations.json.
Unmodeled regimes are rejected, never extrapolated: Mach outside the adopted
curve raises; stations outside the configuration raise at construction.
"""
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np


class MachDragCurve:
    """Piecewise-linear figure readout; reject unmodeled Mach numbers."""
    def __init__(self,mach,coefficients):
        self.mach=np.asarray(mach,dtype=float)
        self.coefficients=np.asarray(coefficients,dtype=float)
        if (self.mach.ndim!=1 or len(self.mach)<2 or self.mach.shape!=self.coefficients.shape
            or not np.all(np.isfinite(self.mach)) or not np.all(np.isfinite(self.coefficients))
            or np.any(np.diff(self.mach)<=0) or self.mach[0]<0 or np.any(self.coefficients<0)):
            raise ValueError('Invalid Mach/drag-coefficient table.')

    @classmethod
    def from_json(cls,path,reference_area_m2):
        record=json.loads(Path(path).read_text(encoding='utf-8'))
        if not np.isclose(record['reference_area_m2'],reference_area_m2,rtol=0,atol=1e-9):
            raise ValueError('Cd curve and reference-area conventions must match.')
        return cls(record['mach'],record['drag_coefficient'])

    def __call__(self,mach):
        if not np.isfinite(mach) or mach<self.mach[0] or mach>self.mach[-1]:
            raise ValueError('Mach number is outside the adopted drag curve.')
        return float(np.interp(mach,self.mach,self.coefficients))


@dataclass(frozen=True)
class AeroOutput:
    force_body_n: np.ndarray
    relative_speed_mps: float
    mach: float
    dynamic_pressure_pa: float
    angle_of_attack_rad: float
    drag_coefficient: float
    normal_force_coefficient: float = 0.0
    free_moment_body_nm: np.ndarray | None = None


class DragOnlyAerodynamics:
    def __init__(self,reference_area_m2,drag_coefficient):
        if not np.isfinite(reference_area_m2) or reference_area_m2 <= 0:
            raise ValueError('Invalid reference area or drag coefficient.')
        if not callable(drag_coefficient) and (not np.isfinite(drag_coefficient) or drag_coefficient<0):
            raise ValueError('Invalid reference area or drag coefficient.')
        self.reference_area_m2=float(reference_area_m2)
        self.drag_coefficient=drag_coefficient

    def evaluate(self,state,properties,atmosphere,wind_body_mps):
        # Velocity of the fixed-body point coincident with CG at this instant,
        # not geometric-CG migration through the structure.
        migration=properties.cg_rate_body_mps
        relative=state.velocity_body_mps-migration-wind_body_mps
        speed=float(np.linalg.norm(relative))
        mach=speed/atmosphere.sound_speed_mps
        cd=float(self.drag_coefficient(mach) if callable(self.drag_coefficient) else self.drag_coefficient)
        if not np.isfinite(cd) or cd<0:
            raise ValueError('Drag coefficient must be finite and nonnegative.')
        dynamic_pressure=.5*atmosphere.density_kgm3*speed*speed
        drag=(-dynamic_pressure*self.reference_area_m2*cd*relative/speed
              if speed > 1e-10 else np.zeros(3))
        aoa=float(np.arctan2(np.linalg.norm(relative[1:]),relative[0])) if speed > 1e-10 else 0.0
        return AeroOutput(drag,speed,mach,dynamic_pressure,aoa,cd)


@dataclass(frozen=True)
class SideSegment:
    """Incremental projected side area with its own crossflow drag.

    Stations use the same structural datum as the phase's mass properties.
    Overlap with body segments is intentional: width_m is the ADDED projected
    width beyond the body (grid fins, deployed legs).
    """
    bottom_m: float
    top_m: float
    width_m: float
    crossflow_drag_coefficient: float

    @property
    def side_area_m2(self):
        return float(self.width_m*(self.top_m-self.bottom_m))

    @property
    def centroid_m(self):
        return float(0.5*(self.bottom_m+self.top_m))


@dataclass(frozen=True)
class AeroConfiguration:
    """One vehicle shape; stations in the phase's structural datum (base=0)."""
    name: str
    reference_area_m2: float
    reference_length_m: float
    axial_drag: MachDragCurve
    potential_normal_slope_per_rad: float
    potential_center_station_m: float
    side_segments: tuple[SideSegment, ...]
    pitch_damping_coefficient: float
    roll_damping_coefficient: float = 0.0

    def __post_init__(self):
        segments=tuple(self.side_segments)
        if (not isinstance(self.name,str) or not self.name
                or not np.all(np.isfinite([self.reference_area_m2,self.reference_length_m,
                                           self.potential_normal_slope_per_rad,
                                           self.potential_center_station_m,
                                           self.pitch_damping_coefficient,
                                           self.roll_damping_coefficient]))
                or self.reference_area_m2 <= 0 or self.reference_length_m <= 0
                or not 0 <= self.potential_normal_slope_per_rad
                or not 0 <= self.potential_center_station_m <= self.reference_length_m
                or self.pitch_damping_coefficient > 0 or self.roll_damping_coefficient > 0):
            raise ValueError('Invalid aerodynamic configuration parameters.')
        if not isinstance(self.axial_drag,MachDragCurve):
            raise ValueError('The axial drag law must be a MachDragCurve.')
        if self.axial_drag.mach[0] > 0:
            raise ValueError('The axial drag curve must cover Mach 0 for the rest case.')
        if not segments:
            raise ValueError('A configuration needs at least one side segment.')
        for segment in segments:
            if (not isinstance(segment,SideSegment)
                    or not np.all(np.isfinite([segment.bottom_m,segment.top_m,
                                               segment.width_m,
                                               segment.crossflow_drag_coefficient]))
                    or not 0 <= segment.bottom_m < segment.top_m <= self.reference_length_m
                    or segment.width_m <= 0
                    or not 0 <= segment.crossflow_drag_coefficient):
                raise ValueError('Invalid configuration side segment.')
        object.__setattr__(self,'side_segments',segments)


class Aerodynamics:
    """Multi-configuration aerodynamic model; the provider selects the shape."""

    def __init__(self,configurations):
        adopted=dict(configurations)
        if not adopted:
            raise ValueError('At least one aerodynamic configuration is required.')
        for name,configuration in adopted.items():
            if not isinstance(configuration,AeroConfiguration) or configuration.name!=name:
                raise ValueError('Configurations must be AeroConfiguration objects keyed by name.')
        self.configurations=adopted

    @classmethod
    def from_json(cls,path):
        """Load declared configuration records; axial curves by file or inline.

        The JSON's own evidence/class/uncertainty notes are the model
        documentation; this loader only checks the numeric contract.
        """
        source=Path(path)
        record=json.loads(source.read_text(encoding='utf-8'))
        base_directory=source.resolve().parent
        configurations={}
        for name,spec in record['configurations'].items():
            curve_record=spec['axial_drag_curve']
            if 'file' in curve_record:
                axial=MachDragCurve.from_json(base_directory/curve_record['file'],
                                              spec['reference_area_m2'])
            else:
                axial=MachDragCurve(curve_record['mach'],curve_record['drag_coefficient'])
            segments=tuple(SideSegment(segment['bottom_m'],segment['top_m'],
                                      segment['width_m'],
                                      segment['crossflow_drag_coefficient'])
                          for segment in spec['side_segments'])
            configurations[name]=AeroConfiguration(
                name=name,
                reference_area_m2=spec['reference_area_m2'],
                reference_length_m=spec['reference_length_m'],
                axial_drag=axial,
                potential_normal_slope_per_rad=spec['potential_normal_slope_per_rad'],
                potential_center_station_m=spec['potential_center_station_m'],
                side_segments=segments,
                pitch_damping_coefficient=spec['pitch_damping_coefficient'],
                roll_damping_coefficient=spec.get('roll_damping_coefficient',0.0),
            )
        return cls(configurations)

    def evaluate(self,state,properties,atmosphere,wind_body_mps,configuration_name):
        configuration=self.configurations.get(configuration_name)
        if configuration is None:
            raise ValueError('Unknown aerodynamic configuration: '+str(configuration_name))
        # Velocity of the fixed-body point coincident with CG at this instant,
        # not geometric-CG migration through the structure.
        migration=properties.cg_rate_body_mps
        relative=state.velocity_body_mps-migration-wind_body_mps
        speed=float(np.linalg.norm(relative))
        reference_area=configuration.reference_area_m2
        if speed<=1e-10:
            return AeroOutput(np.zeros(3),0.0,0.0,0.0,0.0,
                              float(configuration.axial_drag(0.0)),0.0,np.zeros(3))
        direction=relative/speed
        cosine=float(direction[0])
        sine=float(np.linalg.norm(direction[1:]))
        total_angle=float(np.arctan2(sine,cosine))
        mach=speed/atmosphere.sound_speed_mps
        axial_coefficient=float(configuration.axial_drag(mach))
        dynamic_pressure=.5*atmosphere.density_kgm3*speed*speed
        force=np.zeros(3)
        moment=np.zeros(3)
        # Axial force: axial-dynamic-pressure share, directed against the
        # axial relative motion (nose-first or engine-first).
        force[0]-=dynamic_pressure*reference_area*axial_coefficient*cosine*abs(cosine)
        normal_force=0.0
        cg=properties.cg_station_body_m
        if sine>0.0:
            side=np.zeros(3)
            side[1:]=-direction[1:]/sine

            def arm(station_m):
                return np.array([station_m,0.0,0.0])-cg

            potential=dynamic_pressure*reference_area*configuration.potential_normal_slope_per_rad*sine*abs(cosine)
            force+=potential*side
            # Windward end: adopted station nose-first; mirrored about the CG
            # for retro flow; clamped to the modeled body length.
            station=(configuration.potential_center_station_m if cosine>=0.0
                      else 2.0*cg[0]-configuration.potential_center_station_m)
            station=min(max(station,0.0),configuration.reference_length_m)
            moment+=np.cross(arm(station),potential*side)
            normal_force+=potential
            for segment in configuration.side_segments:
                crossflow=(dynamic_pressure*segment.crossflow_drag_coefficient
                           *segment.side_area_m2*sine*sine)
                force+=crossflow*side
                moment+=np.cross(arm(segment.centroid_m),crossflow*side)
                normal_force+=crossflow
        # Rotary damping about the CG; zero when the reference coefficients
        # are zero (declared unmodeled) or the dynamic pressure vanishes.
        damping_scale=dynamic_pressure*reference_area*configuration.reference_length_m**2/(2.0*speed)
        omega=state.omega_body_radps
        moment+=damping_scale*np.array([
            configuration.roll_damping_coefficient*omega[0],
            configuration.pitch_damping_coefficient*omega[1],
            configuration.pitch_damping_coefficient*omega[2]])
        return AeroOutput(force,speed,mach,dynamic_pressure,total_angle,axial_coefficient,
                          normal_force/(dynamic_pressure*reference_area),moment)
