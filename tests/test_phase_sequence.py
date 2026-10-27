"""Analytic shutdown inventory and zero-impulse separation invariants."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from engine_schedule import SmoothShutdown
from propulsion import evaluate_steady_engine,scale_shutdown_output
from VehicleData import build_falcon9
from dynamics import Loads
from staging import separate_coasting_vehicle
from phase_sequence import run_sequence,write_sequence_outputs
from mass_flow import ExhaustPort
from ascent import load_config

ROOT=Path(__file__).resolve().parents[1]


class PhaseSequenceChecks(unittest.TestCase):
    def test_shutdown_derivative_endpoints_and_integral(self):
        schedule=SmoothShutdown(3,2)
        self.assertEqual(schedule(3),(1.,0.));self.assertEqual(schedule(5),(0.,0.))
        for t in (3.2,3.8,4.3,4.9):
            h=1e-5
            self.assertAlmostEqual((schedule(t+h)[0]-schedule(t-h)[0])/(2*h),schedule(t)[1],places=8)
        t=np.linspace(3,5,1001)
        self.assertAlmostEqual(np.trapezoid([schedule(x)[0] for x in t],t),1.,places=10)
        with self.assertRaises(ValueError):SmoothShutdown(1,0)

    def test_shutdown_scales_pressure_corrected_thrust_and_flow_derivative(self):
        engine=build_falcon9().stages[0].engines[0];engine.running=True;engine.throttle=1
        steady=evaluate_steady_engine(engine,101325,stage_index=0,stage_offset_m=0,nozzle_exit_axial_m=0)
        half=scale_shutdown_output(steady,.5,-.2)
        self.assertAlmostEqual(half.thrust_n,.5*845099)
        self.assertAlmostEqual(half.exhaust_port.mass_flow_rate_kgps2,-.2*steady.mass_flow_kgps)
        off=scale_shutdown_output(steady,0,0)
        self.assertEqual(off.thrust_n,0);self.assertEqual(off.mass_flow_kgps,0)

    def test_rotating_separation_preserves_momentum_and_fixed_structural_sensor(self):
        v=build_falcon9();v.position=np.array([1000.,20.,30.]);v.velocity=np.array([100.,2.,3.])
        v.quaternion=np.array([.6,.2,-.3,.7]);v.omega=np.array([.01,.02,-.03])
        state=v.to_state();loads=Loads([0,0,0],[0,0,0],[-9.8,0,0],[0,0])
        split=separate_coasting_vehicle(v,state,12,loads)
        self.assertLess(split.audit['linear_relative_residual'],1e-12)
        self.assertLess(split.audit['angular_relative_residual'],1e-12)
        self.assertLess(split.audit['kinetic_energy_relative_residual'],1e-12)
        np.testing.assert_allclose(split.lower.position,v.position,atol=1e-12)
        np.testing.assert_allclose(split.lower.velocity,v.velocity,atol=1e-12)
        self.assertEqual(split.upper_state.propellant_kg[0],0)
        self.assertEqual(split.lower_state.propellant_kg[1],0)
        # Stage 2 material points inherit omega cross cut, not the same base velocity.
        expected=v.velocity+np.cross(v.omega,[42.6,0,0])
        np.testing.assert_allclose(split.upper.velocity,expected)
        np.testing.assert_allclose(split.lower.quaternion,split.upper.quaternion)
        self.assertEqual(len(v.stages),2)  # map does not mutate the caller

    def test_separation_rejects_running_engine_or_nonzero_outflow(self):
        v=build_falcon9();s=v.to_state();zero=Loads([0,0,0],[0,0,0],[0,0,0],[0,0])
        v.stages[0].engines[0].running=True
        with self.assertRaisesRegex(ValueError,'commands off'):separate_coasting_vehicle(v,s,0,zero)
        v.stages[0].engines[0].running=False
        flowing=Loads([1,0,0],[0,0,0],[0,0,0],[1,0],(ExhaustPort(0,[0,0,0],1),))
        with self.assertRaisesRegex(ValueError,'zero mass flow'):separate_coasting_vehicle(v,s,0,flowing)

    @staticmethod
    def short_config():
        c=load_config(ROOT/'config'/'shutdown_separation_check.json')
        c.update(shutdown_start_s=.5,shutdown_duration_s=1.,separation_time_s=2.,duration_s=2.5,
                 initial_base_position_body_m=[1000.,0,0],initial_base_velocity_body_mps=[100.,0,0],
                 initial_omega_body_radps=[.001,.002,-.003],prescribed_loads=[])
        return c

    def test_sequence_conserves_inventory_and_writes_separate_vehicle_streams(self):
        cfg=self.short_config();rows,summary=run_sequence(cfg,step_s=.05)
        flow=build_falcon9().parameters.stage1_mass_flow_max
        cutoff=summary['shutdown_complete']
        self.assertAlmostEqual(cutoff['stage1_propellant_kg'],410900-flow,places=6)
        self.assertEqual(cutoff['stage1_consumption_kgps'],0)
        self.assertEqual(cutoff['stage1_consumption_derivative_kgps2'],0)
        self.assertEqual(cutoff['thrust_n'],0)
        for name in ('booster','upper_stack'):
            final=summary['final_by_vehicle'][name]
            self.assertEqual(final['time_s'],2.5)
            self.assertEqual(final['thrust_n'],0);self.assertEqual(final['drag_n'],0)
        final=summary['final_by_vehicle']['booster']
        self.assertEqual(final['stage1_propellant_kg'],cutoff['stage1_propellant_kg'])
        self.assertEqual(final['stage2_propellant_kg'],0)
        upper=summary['final_by_vehicle']['upper_stack']
        self.assertEqual(upper['stage1_propellant_kg'],0)
        self.assertEqual(upper['stage2_propellant_kg'],107500)
        event_rows=[r for r in rows if r['time_s']==2]
        self.assertEqual({r['vehicle_id'] for r in event_rows},{'stack','booster','upper_stack'})
        self.assertEqual({r['event_side'] for r in event_rows},{'pre','post'})
        self.assertLess(summary['separation_audit']['angular_relative_residual'],1e-11)
        self.assertEqual(final['imu_id'],summary['initial']['imu_id'])
        self.assertNotEqual(upper['imu_id'],final['imu_id'])
        self.assertAlmostEqual(upper['imu_original_station_x_m'],42.6)
        with tempfile.TemporaryDirectory() as folder:
            write_sequence_outputs(rows,summary,cfg,folder)
            with (Path(folder)/'truth_for_c.csv').open(encoding='utf-8') as f:
                exported=list(csv.DictReader(f))
            self.assertEqual(len(exported),len(rows))
            self.assertEqual({r['vehicle_id'] for r in exported},{'stack','booster','upper_stack'})
            for name in ('stack','booster','upper_stack'):
                with (Path(folder)/f'truth_{name}.csv').open(encoding='utf-8') as f:
                    self.assertEqual({r['vehicle_id'] for r in csv.DictReader(f)},{name})
            events=json.loads((Path(folder)/'events.json').read_text())
            self.assertEqual(events[-1]['event'],'stage_separation')

    def test_shutdown_and_separation_times_must_be_ordered(self):
        cfg=self.short_config();cfg['separation_time_s']=1
        with self.assertRaisesRegex(ValueError,'shutdown end'):run_sequence(cfg)


if __name__=='__main__':unittest.main()
