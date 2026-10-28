"""Component checks for shutdown envelopes and no-impulse separation."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from engine_schedule import SmoothShutdown
from propulsion import evaluate_steady_engine,scale_shutdown_output
from VehicleData import build_falcon9
from dynamics import Loads
from staging import separate_coasting_vehicle
from mass_flow import ExhaustPort


class TransitionMathChecks(unittest.TestCase):
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


if __name__=="__main__":
    unittest.main()
