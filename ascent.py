"""Temporary A/B load provider; mission conditions are supplied by the team."""
from dataclasses import dataclass
import numpy as np
from propulsion import steady_mass_flow,evaluate_steady_engine,scale_shutdown_output
from load_interface import AppliedForce,assemble_loads


@dataclass(frozen=True)
class AscentEvaluation:
    loads: object
    properties: object
    atmosphere: object
    aero: object
    thrust_n: float
    cg_altitude_m: float
    base_altitude_m: float


class SteadyPoweredAscent:
    """Steady engine commands; B can replace this provider without changing A."""
    def __init__(self,vehicle,environment,aerodynamics, *, nozzle_exit_axial_m=0.0,
                 max_angle_of_attack_rad=None,command_envelope=None):
        self.vehicle=vehicle
        self.model=vehicle.mass_model()
        self.environment=environment
        self.aerodynamics=aerodynamics
        self.nozzle_exit_axial_m=float(nozzle_exit_axial_m)
        self.command_envelope=command_envelope
        self.max_angle_of_attack_rad=max_angle_of_attack_rad
        if max_angle_of_attack_rad is not None and (not np.isfinite(max_angle_of_attack_rad) or max_angle_of_attack_rad<=0):
            raise ValueError('A configured angle-of-attack check limit must be positive.')
        self.stage1=next(s for s in vehicle.stages if s.name=='stage1')
        self._commands=tuple((e.running,e.throttle,e.gimbal_pitch,e.gimbal_yaw) for e in self.stage1.engines)
        self.consumption_kgps=np.array([sum(steady_mass_flow(e,vehicle.parameters.G0) for e in self.stage1.engines),0.])

    def flow_at(self,time_s):
        amplitude,derivative=(1.0,0.0) if self.command_envelope is None else self.command_envelope(time_s)
        return self.consumption_kgps*amplitude,self.consumption_kgps*derivative

    def evaluate(self,time_s,state):
        if self._commands != tuple((e.running,e.throttle,e.gimbal_pitch,e.gimbal_yaw) for e in self.stage1.engines):
            raise ValueError('Steady commands changed: construct a new properly bounded phase provider.')
        consumption,consumption_derivative=self.flow_at(time_s)
        properties=self.model.mass_properties(None,state.propellant_kg,consumption,
                                              consumption_rate_kgps2=consumption_derivative)
        altitude=self.environment.altitude(state,properties)
        atmosphere=self.environment.atmosphere(altitude)
        forces=[];ports=[];thrust=0.0
        for engine in self.stage1.engines:
            exit_station=np.array([self.stage1.offset+self.nozzle_exit_axial_m,engine.position[1],engine.position[2]])
            exit_height=self.environment.altitude(state,properties,exit_station)
            ambient=self.environment.atmosphere(exit_height).pressure_pa
            result=evaluate_steady_engine(engine,ambient,stage_index=0,stage_offset_m=self.stage1.offset,
                                           nozzle_exit_axial_m=self.nozzle_exit_axial_m,g0=self.vehicle.parameters.G0)
            if self.command_envelope is not None:
                result=scale_shutdown_output(result,*self.command_envelope(time_s))
            forces.append(result.applied_force);ports.append(result.exhaust_port);thrust+=result.thrust_n
        aero=self.aerodynamics.evaluate(state,properties,atmosphere,self.environment.wind_body(state))
        if self.max_angle_of_attack_rad is not None and aero.relative_speed_mps>1e-6 and aero.angle_of_attack_rad>self.max_angle_of_attack_rad:
            raise ValueError('Attitude check exceeded its declared angle-of-attack limit.')
        # Drag-through-CG is an explicit first-model assumption: no aero moment.
        forces.append(AppliedForce(aero.force_body_n,properties.cg_station_body_m))
        loads=assemble_loads(properties,forces,gravity_body_mps2=self.environment.gravity_body(state),
                              consumption_kgps=consumption,exhaust_ports=ports)
        return AscentEvaluation(loads,properties,atmosphere,aero,thrust,altitude,
                                 self.environment.altitude(state,properties,[0,0,0]))

    def __call__(self,time_s,state):
        return self.evaluate(time_s,state).loads
