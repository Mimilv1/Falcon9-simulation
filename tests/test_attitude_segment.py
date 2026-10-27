"""Physical signs, sensor separation and temporal convergence of connected rotation."""
import unittest
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ascent import load_config,run_ascent,SteadyPoweredAscent
from prescribed_loads import SmoothForceDoublet
from VehicleData import build_falcon9
from Environement import TroposphereEnvironment
from Aerodynamics import DragOnlyAerodynamics
from dynamics import acceleration_truth

ROOT=Path(__file__).resolve().parents[1]


class AttitudeSegmentChecks(unittest.TestCase):
    def test_doublet_is_bounded_balanced_and_smooth_at_endpoints(self):
        pulse=SmoothForceDoublet(10,4,[0,10000,0],[1,0,0])
        self.assertEqual(pulse(9).force_body_n[1],0)
        self.assertEqual(pulse(14).force_body_n[1],0)
        self.assertAlmostEqual(pulse(11).force_body_n[1],10000)
        self.assertAlmostEqual(pulse(13).force_body_n[1],-10000)
        times=np.linspace(10,14,1001)
        impulse=np.trapezoid([pulse(t).force_body_n[1] for t in times],times)
        self.assertAlmostEqual(impulse,0,places=9)
        self.assertLess(np.linalg.norm(pulse(10+1e-5).force_body_n)/1e-5,1e-4)
        with self.assertRaises(ValueError):SmoothForceDoublet(0,0,[0,1,0],[0,0,0])

    def test_applied_force_station_gives_expected_rotation_signs(self):
        for axis,force in ((2,[0,10000,0]),(1,[0,0,10000])):
            v=build_falcon9()
            for e in v.stages[0].engines:e.running=True;e.throttle=1
            v.position=np.array([1000.,0,0]);v.velocity=np.array([100.,0,0])
            pulse=SmoothForceDoublet(0,4,force,[1,0,0])
            provider=SteadyPoweredAscent(v,TroposphereEnvironment(),DragOnlyAerodynamics(21.08,.1),
                                         prescribed_loads=[pulse])
            state=v.to_state(consumption_by_stage={'stage1':provider.consumption_kgps[0]})
            out=provider.evaluate(1,state)
            expected=np.cross(np.array([1,0,0])-out.properties.cg_station_body_m,force)
            np.testing.assert_allclose(out.loads.moment_body_nm,expected,atol=1e-8)
            sensor=acceleration_truth(state,out.loads,out.properties,[0,0,0])
            np.testing.assert_allclose(sensor['angular_acceleration_body_radps2'],
                                       expected/np.diag(out.properties.inertia_body_kgm2),atol=1e-14)
            if axis==2:self.assertLess(expected[axis],0)
            else:self.assertGreater(expected[axis],0)

    @staticmethod
    def short_config():
        cfg=load_config(ROOT/'config'/'ascent_attitude_check.json')
        cfg['duration_s']=3
        cfg['initial_base_position_body_m']=[1000,0,0]
        cfg['initial_base_velocity_body_mps']=[100,0,0]
        cfg['prescribed_loads'][0].update(start_s=.5,duration_s=1)
        return cfg

    def test_rotating_segment_sensor_lever_arm_and_gravity_sign(self):
        rows,summary=run_ascent(self.short_config(),step_s=.05)
        peak=max(rows,key=lambda r:abs(r['angular_acceleration_body_radps2_z']))
        omega=np.array([peak[f'gyro_ideal_body_radps_{a}'] for a in 'xyz'])
        alpha=np.array([peak[f'angular_acceleration_body_radps2_{a}'] for a in 'xyz'])
        measured=np.array([peak[f'comparison_accelerometer_ideal_body_mps2_{a}']-
                           peak[f'accelerometer_ideal_body_mps2_{a}'] for a in 'xyz'])
        arm=np.array([40.,0,0])
        np.testing.assert_allclose(measured,np.cross(alpha,arm)+np.cross(omega,np.cross(omega,arm)),atol=1e-12)
        self.assertGreater(np.linalg.norm(measured),.001)
        final=summary['final']
        self.assertLess(final['quaternion_body_to_initial_z'],0)
        self.assertLess(final['gravity_body_mps2_y'],0)
        self.assertGreater(summary['maximum_tilt_deg'],0)
        self.assertLess(summary['maximum_quaternion_norm_error'],1e-12)
        for boundary in [.5,1,1.5]:
            self.assertTrue(any(abs(r['time_s']-boundary)<1e-12 for r in rows))

    def test_force_reversal_mirrors_attitude_and_step_refinement_converges(self):
        cfg=self.short_config()
        _,coarse=run_ascent(cfg,step_s=.1)
        _,fine=run_ascent(cfg,step_s=.05)
        _,reference=run_ascent(cfg,step_s=.025)
        key='quaternion_body_to_initial_z'
        coarse_error=abs(coarse['final'][key]-reference['final'][key])
        fine_error=abs(fine['final'][key]-reference['final'][key])
        self.assertLess(fine_error,coarse_error/8)
        self.assertLess(fine_error,1e-8)
        cfg['prescribed_loads'][0]['peak_force_body_n'][1]*=-1
        _,mirror=run_ascent(cfg,step_s=.05)
        self.assertAlmostEqual(mirror['final']['quaternion_body_to_initial_z'],
                               -fine['final']['quaternion_body_to_initial_z'],places=11)
        self.assertAlmostEqual(mirror['final']['base_altitude_m'],fine['final']['base_altitude_m'],places=9)

    def test_selected_angle_of_attack_guard_is_enforced(self):
        cfg=self.short_config();cfg['max_angle_of_attack_deg']=1e-8
        with self.assertRaisesRegex(ValueError,'angle-of-attack'):
            run_ascent(cfg,step_s=.1)


if __name__=='__main__':unittest.main()
