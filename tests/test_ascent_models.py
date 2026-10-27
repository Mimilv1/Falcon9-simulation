"""First B models checked against standard values and physical limits."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Environement import TroposphereEnvironment
from propulsion import evaluate_steady_engine
from Aerodynamics import DragOnlyAerodynamics,MachDragCurve
from VehicleData import build_falcon9
from dynamics import State
from ascent import load_config,run_ascent,SteadyPoweredAscent


class AscentModelChecks(unittest.TestCase):
    def test_drag_table_interpolation_and_supported_domain(self):
        curve=MachDragCurve([0,1,2],[.1,1.1,.9])
        self.assertAlmostEqual(curve(.5),.6)
        with self.assertRaises(ValueError):curve(2.01)
        with self.assertRaises(ValueError):MachDragCurve([0,0],[.1,.2])
        source=Path(__file__).resolve().parents[1]/'data'/'drag_curve_pezent.json'
        with self.assertRaises(ValueError):MachDragCurve.from_json(source,10)
        actual=MachDragCurve.from_json(source,21.08)
        self.assertGreater(actual(1.06),1.1)
        self.assertLess(actual(.4),.15)

    def test_standard_atmosphere_reference_points(self):
        a=TroposphereEnvironment.atmosphere(0)
        self.assertEqual(a.temperature_k,288.15)
        self.assertEqual(a.pressure_pa,101325)
        self.assertAlmostEqual(a.density_kgm3,1.225,places=5)
        a=TroposphereEnvironment.atmosphere(1000)
        self.assertAlmostEqual(a.temperature_k,281.65)
        self.assertAlmostEqual(a.pressure_pa,89874.5705,delta=.05)
        self.assertAlmostEqual(a.density_kgm3,1.1116,delta=1e-4)

    def test_environment_rejects_unimplemented_altitudes(self):
        with self.assertRaises(ValueError):TroposphereEnvironment.atmosphere(11001)
        with self.assertRaises(ValueError):TroposphereEnvironment.atmosphere(-1)

    def test_gravity_changes_body_components_when_attitude_turns(self):
        state=State([0,0,0],[0,0,0],[np.sqrt(.5),0,0,np.sqrt(.5)],[0,0,0],[0,0])
        np.testing.assert_allclose(TroposphereEnvironment().gravity_body(state),[0,9.80665,0],atol=1e-12)

    def test_thrust_vacuum_limit_and_pressure_correction(self):
        e=build_falcon9().stages[0].engines[0];e.running=True;e.throttle=1
        args=dict(stage_index=0,stage_offset_m=0,nozzle_exit_axial_m=0)
        vac=evaluate_steady_engine(e,0,**args);sea=evaluate_steady_engine(e,101325,**args)
        self.assertEqual(vac.thrust_n,914000)
        self.assertAlmostEqual(vac.thrust_n-sea.thrust_n,.68*101325)
        self.assertEqual(vac.mass_flow_kgps,sea.mass_flow_kgps)

    def test_drag_opposes_relative_velocity_and_vanishes_at_rest(self):
        v=build_falcon9();p=v.get_mass_properties();a=TroposphereEnvironment.atmosphere(0)
        aero=DragOnlyAerodynamics(21.08,.8)
        state=State([100,0,0],[100,20,-5],[1,0,0,0],[0,0,0],[410900,107500])
        result=aero.evaluate(state,p,a,np.zeros(3))
        self.assertLess(np.dot(result.force_body_n,state.velocity_body_mps),0)
        rest=State([100,0,0],[0,0,0],[1,0,0,0],[0,0,0],[410900,107500])
        np.testing.assert_allclose(aero.evaluate(rest,p,a,np.zeros(3)).force_body_n,[0,0,0])

    def test_short_project_segment_produces_sensor_truth_and_climb(self):
        cfg=load_config(Path(__file__).resolve().parents[1]/'config'/'powered_ascent.json')
        cfg['duration_s']=1
        rows,summary=run_ascent(cfg,step_s=.1)
        self.assertEqual(summary['stop_reason'],'end_time')
        self.assertGreater(summary['final']['base_altitude_m'],0)
        self.assertGreater(summary['final']['thrust_n'],summary['initial']['thrust_n'])
        flow=build_falcon9().parameters.stage1_mass_flow_max
        self.assertAlmostEqual(summary['final']['total_mass_kg'],556300-flow,places=7)
        self.assertAlmostEqual(summary['initial']['pressure_imu_pa'],101325)
        self.assertGreater(summary['initial']['accelerometer_ideal_body_mps2_x'],9.80665)
        self.assertIn('quaternion_body_to_initial_w',rows[0])
        self.assertIn('gyro_ideal_body_radps_x',rows[0])
        self.assertLess(summary['maximum_quaternion_norm_error'],1e-12)


if __name__=='__main__':unittest.main()
