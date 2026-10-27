"""A+B connection for an early, vertical powered-ascent segment."""
from dataclasses import dataclass
from pathlib import Path
import csv,hashlib,json,subprocess
import numpy as np
from VehicleData import build_falcon9
from Environement import TroposphereEnvironment
from Aerodynamics import DragOnlyAerodynamics,MachDragCurve
from propulsion import steady_mass_flow,evaluate_steady_engine
from load_interface import AppliedForce,assemble_loads
from dynamics import acceleration_truth
from propagation import propagate


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
    def __init__(self,vehicle,environment,aerodynamics, *, nozzle_exit_axial_m=0.0):
        self.vehicle=vehicle
        self.model=vehicle.mass_model()
        self.environment=environment
        self.aerodynamics=aerodynamics
        self.nozzle_exit_axial_m=float(nozzle_exit_axial_m)
        self.stage1=next(s for s in vehicle.stages if s.name=='stage1')
        self._commands=tuple((e.running,e.throttle,e.gimbal_pitch,e.gimbal_yaw) for e in self.stage1.engines)
        self.consumption_kgps=np.array([sum(steady_mass_flow(e,vehicle.parameters.G0) for e in self.stage1.engines),0.])

    def evaluate(self,time_s,state):
        if self._commands != tuple((e.running,e.throttle,e.gimbal_pitch,e.gimbal_yaw) for e in self.stage1.engines):
            raise ValueError('Steady commands changed: construct a new properly bounded phase provider.')
        properties=self.model.mass_properties(None,state.propellant_kg,self.consumption_kgps)
        altitude=self.environment.altitude(state,properties)
        atmosphere=self.environment.atmosphere(altitude)
        forces=[];ports=[];thrust=0.0
        for engine in self.stage1.engines:
            exit_station=np.array([self.stage1.offset+self.nozzle_exit_axial_m,engine.position[1],engine.position[2]])
            exit_height=self.environment.altitude(state,properties,exit_station)
            ambient=self.environment.atmosphere(exit_height).pressure_pa
            result=evaluate_steady_engine(engine,ambient,stage_index=0,stage_offset_m=self.stage1.offset,
                                           nozzle_exit_axial_m=self.nozzle_exit_axial_m,g0=self.vehicle.parameters.G0)
            forces.append(result.applied_force);ports.append(result.exhaust_port);thrust+=result.thrust_n
        aero=self.aerodynamics.evaluate(state,properties,atmosphere,self.environment.wind_body(state))
        # Drag-through-CG is an explicit first-model assumption: no aero moment.
        forces.append(AppliedForce(aero.force_body_n,properties.cg_station_body_m))
        loads=assemble_loads(properties,forces,gravity_body_mps2=self.environment.gravity_body(state),
                              consumption_kgps=self.consumption_kgps,exhaust_ports=ports)
        return AscentEvaluation(loads,properties,atmosphere,aero,thrust,altitude,
                                 self.environment.altitude(state,properties,[0,0,0]))

    def __call__(self,time_s,state):
        return self.evaluate(time_s,state).loads


def load_config(path):
    with open(path,encoding='utf-8') as f: return json.load(f)


def run_ascent(config, *, step_s=None):
    if config['phase']!='early_vertical_powered_ascent':
        raise ValueError('Only the first powered-ascent segment is implemented.')
    vehicle=build_falcon9()
    throttle=float(config['stage1_throttle'])
    if not vehicle.stages[0].engines[0].spec.throttle_min <= throttle <= 1:
        raise ValueError('The powered-ascent segment requires steady running engines within throttle limits.')
    for engine in vehicle.stages[0].engines:
        engine.running=True;engine.throttle=throttle
    vehicle.position=np.asarray(config['initial_base_position_body_m'],dtype=float)
    vehicle.velocity=np.asarray(config['initial_base_velocity_body_mps'],dtype=float)
    vehicle.omega=np.asarray(config['initial_omega_body_radps'],dtype=float)
    environment=TroposphereEnvironment(gravity_mps2=config['gravity_mps2'],
        ground_altitude_m=config['ground_altitude_m'],up_initial_body=config['up_initial_body'],
        wind_initial_body_mps=config['wind_initial_body_mps'])
    curve=MachDragCurve.from_json(Path(__file__).parent/config['drag_curve_file'],config['aero_reference_area_m2'])
    aero=DragOnlyAerodynamics(config['aero_reference_area_m2'],curve)
    provider=SteadyPoweredAscent(vehicle,environment,aero,nozzle_exit_axial_m=config['nozzle_exit_axial_m'])
    initial=vehicle.to_state(consumption_by_stage={'stage1':provider.consumption_kgps[0]})
    used_step=float(config['step_s'] if step_s is None else step_s)
    result=propagate(0,initial,float(config['duration_s']),used_step,
                     provider.model,None,provider)
    rows=[]
    sensor_station=np.asarray(config['imu_station_body_m'],dtype=float)
    for t,state in zip(result.times_s,result.states):
        sample=provider.evaluate(t,state);p=sample.properties
        sensor=acceleration_truth(state,sample.loads,p,sensor_station_body_m=sensor_station)
        sensor_altitude=environment.altitude(state,p,sensor_station)
        pressure=environment.atmosphere(sensor_altitude).pressure_pa
        inventory=vehicle.stages[0].liquid_profile.stage1_inventory(state.propellant_kg[0])
        row={'time_s':float(t),'phase':config['phase'],'total_mass_kg':p.mass_kg,
             'stage1_propellant_kg':float(state.propellant_kg[0]),'rp1_kg':inventory[0],'lox_kg':inventory[1],
             'cg_station_x_m':float(p.cg_station_body_m[0]),'base_altitude_m':sample.base_altitude_m,
             'cg_altitude_m':sample.cg_altitude_m,'airspeed_mps':sample.aero.relative_speed_mps,
             'mach':sample.aero.mach,'drag_coefficient':sample.aero.drag_coefficient,
             'dynamic_pressure_pa':sample.aero.dynamic_pressure_pa,
             'angle_of_attack_rad':sample.aero.angle_of_attack_rad,'thrust_n':sample.thrust_n,
             'drag_n':float(np.linalg.norm(sample.aero.force_body_n)),
             'pressure_imu_pa':pressure,'density_cg_kgm3':sample.atmosphere.density_kgm3,
             'temperature_cg_k':sample.atmosphere.temperature_k}
        for prefix,axes,values in [
            ('position_body_m','xyz',state.position_body_m),('velocity_body_mps','xyz',state.velocity_body_mps),
            ('quaternion_body_to_initial','wxyz',state.quaternion_body_to_initial_body_wxyz),
            ('gyro_ideal_body_radps','xyz',state.omega_body_radps),
            ('accelerometer_ideal_body_mps2','xyz',sensor['specific_force_body_mps2']),
            ('acceleration_body_mps2','xyz',sensor['acceleration_body_mps2']),
            ('gravity_body_mps2','xyz',sample.loads.gravity_body_mps2),
            ('inertia_body_kgm2','xyz',np.diag(p.inertia_body_kgm2))]:
            row.update({f'{prefix}_{axis}':float(value) for axis,value in zip(axes,values)})
        rows.append(row)
    final_loads=provider(result.times_s[-1],result.states[-1])
    vehicle.apply_state(result.states[-1],result.times_s[-1],
                         consumption_by_stage={'stage1':final_loads.consumption_kgps[0]})
    summary={'phase':config['phase'],'stop_reason':result.stop_reason,'samples':len(rows),
             'step_s':used_step,'initial':rows[0],'final':rows[-1],
             'maximum_dynamic_pressure_pa':max(r['dynamic_pressure_pa'] for r in rows),
             'maximum_mach':max(r['mach'] for r in rows),
             'maximum_quaternion_norm_error':max(abs(np.linalg.norm(s.quaternion_body_to_initial_body_wxyz)-1) for s in result.states),
             'sensor_note':'Ideal truth only. C adds sensor bias/noise/sampling; the estimator must not consume continuous truth.',
             'scope':'First 30 s/selected duration of vertical ascent, not all flight phases or real-flight validation.',
             'assumptions':config['provenance']}
    return rows,summary


def write_outputs(rows,summary,config,output_dir):
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    config={**config,'step_s':summary['step_s']}
    with (output_dir/'truth_for_c.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output_dir/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    (output_dir/'run_config.json').write_text(json.dumps(config,indent=2,ensure_ascii=False),encoding='utf-8')
    source=Path(__file__).parent/'data'/'vehicle_parameters.json'
    (output_dir/'vehicle_sources.json').write_text(source.read_text(encoding='utf-8'),encoding='utf-8')
    curve=Path(__file__).parent/config['drag_curve_file']
    (output_dir/'drag_curve_sources.json').write_bytes(curve.read_bytes())
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=Path(__file__).parent,capture_output=True,text=True)
    status=subprocess.run(['git','status','--porcelain'],cwd=Path(__file__).parent,capture_output=True,text=True)
    source_paths=list(Path(__file__).parent.glob('*.py'))+[source,curve]
    manifest={'code_revision':revision.stdout.strip() if revision.returncode==0 else None,
              'working_tree_dirty':bool(status.stdout.strip()) if status.returncode==0 else None,
              'source_sha256':{str(p.relative_to(Path(__file__).parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source_paths)},
              'coordinate_frame':'current Body, +x nose; geographic transforms deferred',
              'imu_station_body_m':config['imu_station_body_m'],'phase':config['phase']}
    (output_dir/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
