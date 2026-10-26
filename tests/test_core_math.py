"""Analytical core checks; synthetic fixtures, no vehicle-data import."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mass_properties import fixed_cylinder, combine_components
from frames import body_to_ned_matrix
from load_interface import AppliedForce, assemble_loads
from dynamics import State, Loads, rk4_step


class FixedModel:
    mass_model_name = 'fixed_diagnostic'

    def __init__(self):
        self.properties = combine_components([fixed_cylinder('synthetic', 1000, 0, 20, 40, 1)])

    def mass_properties(self, configuration, propellant_kg, consumption_kgps=(0, 0)):
        return self.properties


class CoreMathChecks(unittest.TestCase):
    def test_parallel_axis_mass_combination(self):
        p = combine_components([fixed_cylinder('a',100,0,0,0,.1),
                                fixed_cylinder('b',100,0,2,0,.1)])
        self.assertEqual(p.mass_kg, 200)
        np.testing.assert_allclose(p.cg_station_body_m, [0,0,1])
        np.testing.assert_allclose(p.inertia_body_kgm2, np.diag([200.5,200.5,1]))

    def test_base_side_force_has_negative_y_moment(self):
        p = FixedModel().properties
        loads = assemble_loads(p,[AppliedForce([1000,0,0],[0,0,0])],gravity_ned_mps2=[0,0,0])
        np.testing.assert_allclose(loads.moment_body_nm,[0,-20000,0])

    def test_vertical_body_axes_map_to_ned(self):
        r = body_to_ned_matrix([0,np.sqrt(.5),np.sqrt(.5),0])
        np.testing.assert_allclose(r,[[0,1,0],[1,0,0],[0,0,-1]],atol=1e-14)

    def test_free_fall(self):
        model = FixedModel()
        state = State([0,0,-100],[0,0,0],[0,np.sqrt(.5),np.sqrt(.5),0],[0,0,0],[0,0])
        provider = lambda t,s: Loads([0,0,0],[0,0,0],[0,0,10],[0,0])
        result = rk4_step(0,state,1,model,None,provider)
        np.testing.assert_allclose(result.position_ned_m,[0,0,-95],atol=1e-12)
        np.testing.assert_allclose(body_to_ned_matrix(result.quaternion_body_to_ned_wxyz)
                                   @ result.velocity_body_mps,[0,0,10],atol=1e-12)

    def test_constant_moment_rotation(self):
        model = FixedModel()
        p = model.properties
        state = State([0,0,-100],[0,0,0],[0,np.sqrt(.5),np.sqrt(.5),0],[0,0,0],[0,0])
        provider = lambda t,s: assemble_loads(p,[AppliedForce([1000,0,0],[0,0,0])],gravity_ned_mps2=[0,0,0])
        alpha = -20000/p.inertia_body_kgm2[1,1]
        for k in range(200): state = rk4_step(k*.01,state,.01,model,None,provider)
        nose = body_to_ned_matrix(state.quaternion_body_to_ned_wxyz)[:,2]
        tilt = np.arctan2(np.linalg.norm(nose[:2]),-nose[2])
        self.assertAlmostEqual(tilt,abs(alpha)*2,places=10)
        self.assertAlmostEqual(state.omega_body_radps[1],alpha*2,places=10)


if __name__ == '__main__': unittest.main()
