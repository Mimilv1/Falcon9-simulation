"""Team-object/source-profile contracts and staging physics."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from VehicleData import build_falcon9, VehicleModel, Configuration
from frames import body_to_enu_matrix
from dynamics import Loads, rk4_step


class VehicleIntegrationChecks(unittest.TestCase):
    def test_factory_and_engine_specs(self):
        v = build_falcon9()
        self.assertEqual(len(v.stages[0].engines),9)
        self.assertIs(v.stages[0].engines[0].spec,v.stages[0].engines[1].spec)
        self.assertGreater(v.stages[0].engines[0].spec.isp_sl,0)
        self.assertIsNone(v.stages[1].engines[0].spec.isp_sl)
        self.assertEqual(v.stages[1].engines[0].spec.thrust_vac,935000)

    def test_total_mass_without_duplicate_engines(self):
        v = build_falcon9()
        self.assertEqual(v.getmass(),556300)
        self.assertEqual(v.get_mass_properties().mass_kg,v.getmass())

    def test_dry_budget_and_axis_order(self):
        stage = build_falcon9().stages[0]
        self.assertAlmostEqual(stage.cg_dry,17.52107927927928)
        self.assertLess(stage.inertia_dry[0,0],stage.inertia_dry[1,1])
        self.assertEqual(stage.inertia_dry[1,1],stage.inertia_dry[2,2])

    def test_current_objects_agree_with_source_profile(self):
        v = build_falcon9()
        p = v.get_mass_properties()
        m = VehicleModel(mass_model='dual_tank')
        expected = m.mass_properties(Configuration.FULL_STACK,[410900,107500])
        np.testing.assert_allclose(p.cg_station_body_m,expected.cg_station_body_m)
        np.testing.assert_allclose(p.inertia_body_kgm2,expected.inertia_body_kgm2)

    def test_limiting_tank_and_capacity(self):
        s = build_falcon9().stages[0]
        minimum = s.liquid_profile.stage1_burn_limits()['minimum_remaining_kg']
        s.propellant_mass = minimum
        inventories = s.tank_inventories()
        self.assertAlmostEqual(inventories['rp1'],1202.1276595745)
        self.assertAlmostEqual(inventories['lox'],0)
        s.mass_properties()
        s.propellant_mass = minimum-1
        with self.assertRaises(ValueError): s.mass_properties()

    def test_separation_preserves_structural_points_and_mass(self):
        v = build_falcon9()
        v.position=np.array([10.,20.,30.]); v.velocity=np.array([1.,2.,3.])
        v.omega=np.array([.02,-.03,.04]); v.quaternion=np.array([.6,.2,-.3,.7])
        rotation=body_to_enu_matrix(v.quaternion)
        cut=v.stages[1].offset
        local_point=np.array([.5,0,0]); original_point=local_point+np.array([cut,0,0])
        lower,upper=v.separate()
        np.testing.assert_allclose(v.position+rotation@original_point,upper.position+rotation@local_point)
        np.testing.assert_allclose(v.velocity+rotation@np.cross(v.omega,original_point),
                                   upper.velocity+rotation@np.cross(upper.omega,local_point))
        self.assertEqual(lower.getmass()+upper.getmass(),v.getmass())
        original_moment=v.get_mass_properties().mass_kg*v.get_mass_properties().cg_station_body_m
        split_moment=(lower.getmass()*lower.get_mass_properties().cg_station_body_m
                      +upper.getmass()*(upper.get_mass_properties().cg_station_body_m+np.array([cut,0,0])))
        np.testing.assert_allclose(original_moment,split_moment)

    def test_trial_evaluation_does_not_mutate_shared_state(self):
        v=build_falcon9()
        adapter=v.mass_model()
        props=adapter.mass_properties(None,[200000,107500])
        self.assertEqual(v.stages[0].propellant_mass,410900)
        self.assertAlmostEqual(props.mass_kg,v.getmass()-(410900-200000))

    def test_base_to_cg_mapping_and_no_dual_tank_powered_claim(self):
        v=build_falcon9()
        state=v.to_state()
        self.assertGreater(state.position_enu_m[2],0)
        np.testing.assert_allclose(state.position_enu_m[:2],[0,0])
        provider=lambda t,s: Loads([0,0,0],[0,0,0],[0,0,-9.80665],[1,0])
        with self.assertRaisesRegex(ValueError,'open-system'):
            rk4_step(0,state,.01,v.mass_model(),None,provider)


if __name__ == '__main__': unittest.main()
