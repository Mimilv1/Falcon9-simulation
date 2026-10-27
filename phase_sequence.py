"""A shutdown/coast/separation verification sequence, with two tracked branches."""
from pathlib import Path
import csv,json
import numpy as np
from VehicleData import build_falcon9
from Environement import TroposphereEnvironment
from Aerodynamics import DragOnlyAerodynamics,MachDragCurve
from ascent import SteadyPoweredAscent,AscentEvaluation,sample_truth,write_outputs
from dynamics import Loads
from engine_schedule import SmoothShutdown
from prescribed_loads import force_doublet_from_config
from propagation import propagate
from staging import separate_coasting_vehicle


class GravityOnlyCoast:
    """Frozen liquid geometry, constant gravity, no drag: verification scope only."""
    def __init__(self,vehicle,environment):
        if any(e.running or e.throttle!=0 for s in vehicle.stages for e in s.engines):
            raise ValueError('Gravity-only coast requires off engine commands.')
        self.vehicle=vehicle;self.model=vehicle.mass_model();self.environment=environment
        self.aerodynamics=DragOnlyAerodynamics(1.,0.)

    def evaluate(self,time_s,state):
        p=self.model.mass_properties(None,state.propellant_kg)
        altitude=self.environment.altitude(state,p)
        atmosphere=self.environment.atmosphere(altitude)
        aero=self.aerodynamics.evaluate(state,p,atmosphere,self.environment.wind_body(state))
        loads=Loads([0,0,0],[0,0,0],self.environment.gravity_body(state),[0,0])
        return AscentEvaluation(loads,p,atmosphere,aero,0.,altitude,
                                 self.environment.altitude(state,p,[0,0,0]),np.zeros(3))

    def __call__(self,time_s,state):return self.evaluate(time_s,state).loads


def run_sequence(config, *, step_s=None):
    if config['phase']!='shutdown_separation_check':raise ValueError('Unknown sequence.')
    shutdown=SmoothShutdown(config['shutdown_start_s'],config['shutdown_duration_s'])
    separation=float(config['separation_time_s']);end=float(config['duration_s'])
    if not np.all(np.isfinite([separation,end])) or not shutdown.end_s<separation<end:
        raise ValueError('Require shutdown end < separation < simulation end.')
    if config['post_separation_model']!='gravity_only_frozen_liquids':
        raise ValueError('Separate aerodynamic/re-settling models have not been provided.')
    vehicle=build_falcon9()
    throttle=float(config['stage1_throttle'])
    for engine in vehicle.stages[0].engines:engine.running=True;engine.throttle=throttle
    vehicle.position=np.asarray(config['initial_base_position_body_m'],dtype=float)
    vehicle.velocity=np.asarray(config['initial_base_velocity_body_mps'],dtype=float)
    vehicle.omega=np.asarray(config['initial_omega_body_radps'],dtype=float)
    env=TroposphereEnvironment(gravity_mps2=config['gravity_mps2'],ground_altitude_m=config['ground_altitude_m'],
                              up_initial_body=config['up_initial_body'],wind_initial_body_mps=config['wind_initial_body_mps'])
    curve=MachDragCurve.from_json(Path(__file__).parent/config['drag_curve_file'],config['aero_reference_area_m2'])
    prescribed=tuple(force_doublet_from_config(r) for r in config.get('prescribed_loads',[]))
    if any(p.start_s+p.duration_s>shutdown.end_s for p in prescribed):
        raise ValueError('Prescribed attitude loads must finish before shutdown ends.')
    provider=SteadyPoweredAscent(vehicle,env,DragOnlyAerodynamics(config['aero_reference_area_m2'],curve),
        nozzle_exit_axial_m=config['nozzle_exit_axial_m'],prescribed_loads=prescribed,
        max_angle_of_attack_rad=np.deg2rad(config['max_angle_of_attack_deg']),command_envelope=shutdown)
    q,dq=provider.flow_at(0)
    initial=vehicle.to_state(consumption_by_stage={'stage1':q[0]},consumption_derivative_by_stage={'stage1':dq[0]})
    step=float(config['step_s'] if step_s is None else step_s)
    result=propagate(0,initial,separation,step,provider.model,None,provider,
        breakpoints_s=[*shutdown.breakpoints_s,*(t for p in prescribed for t in p.breakpoints_s)])
    def phase(t):
        return 'powered_ascent' if t<shutdown.start_s else ('shutdown' if t<shutdown.end_s else 'coast_attached')
    rows=[sample_truth(vehicle,provider,t,s,config,phase=phase(t),event_side='pre' if t==separation else 'regular')
          for t,s in zip(result.times_s,result.states)]
    if result.stop_reason!='end_time':
        raise RuntimeError('Propellant limit reached before scheduled separation; reschedule shutdown first.')
    separation_loads=provider(separation,result.states[-1])
    # Commands change only after accepted propagation, never during RK trials.
    for stage in vehicle.stages:
        for engine in stage.engines:engine.running=False;engine.throttle=0.
    split=separate_coasting_vehicle(vehicle,result.states[-1],separation,separation_loads)
    final_by_vehicle={'stack':rows[-1]}
    for name,child,state in [('booster',split.lower,split.lower_state),('upper_stack',split.upper,split.upper_state)]:
        coast=GravityOnlyCoast(child,env)
        branch_cfg=dict(config)
        if name=='upper_stack':
            branch_cfg.update(imu_station_body_m=config['upper_imu_station_body_m'],comparison_imu_station_body_m=None,
                              imu_id='upper_base_imu')
        branch=propagate(separation,state,end,step,coast.model,None,coast)
        branch_rows=[sample_truth(child,coast,t,s,branch_cfg,phase='coast_separated',vehicle_id=name,
                                 event_side='post' if t==separation else 'regular')
                     for t,s in zip(branch.times_s,branch.states)]
        rows.extend(branch_rows);final_by_vehicle[name]=branch_rows[-1]
        child.apply_state(branch.states[-1],branch.times_s[-1])
    rows.sort(key=lambda row:(row['time_s'],0 if row['vehicle_id']=='stack' else 1,row['vehicle_id']))
    cutoff=next(row for row in rows if row['vehicle_id']=='stack' and row['time_s']==shutdown.end_s)
    summary={'phase':config['phase'],'stop_reason':'end_time','samples':len(rows),'step_s':step,
             'initial':rows[0],'final':final_by_vehicle['booster'],'final_by_vehicle':final_by_vehicle,
             'shutdown_complete':cutoff,'separation_audit':split.audit,
             'vehicle_sensor_stations':{
                 'stack':{'imu_id':config.get('imu_id','original_base_imu'),'own_base_station_m':config['imu_station_body_m']},
                 'booster':{'imu_id':config.get('imu_id','original_base_imu'),'own_base_station_m':config['imu_station_body_m']},
                 'upper_stack':{'imu_id':'upper_base_imu','own_base_station_m':config['upper_imu_station_body_m']}},
             'events':[{'time_s':shutdown.start_s,'event':'shutdown_start'},
                       {'time_s':shutdown.end_s,'event':'shutdown_complete'},split.audit],
             'assumptions':config['provenance'],
             'sensor_note':'Keep vehicle_id streams separate. Booster retains the original base IMU; upper_stack has an independent virtual IMU at its own base.',
             'scope':'A-class event verification; arbitrary early shutdown/separation timing, frozen liquids and gravity-only child coasts. Not a Falcon 9 mission.'}
    return rows,summary


def write_sequence_outputs(rows,summary,config,output_dir):
    write_outputs(rows,summary,config,output_dir)
    output=Path(output_dir)
    for vehicle_id in ('stack','booster','upper_stack'):
        selected=[r for r in rows if r['vehicle_id']==vehicle_id]
        fields=list(dict.fromkeys(key for row in selected for key in row))
        with (output/f'truth_{vehicle_id}.csv').open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(selected)
    (output/'events.json').write_text(json.dumps(summary['events'],indent=2)+'\n',encoding='utf-8')
