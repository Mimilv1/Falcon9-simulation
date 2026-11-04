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
from ascent import SteadyPoweredAscent


class AscentModelChecks(unittest.TestCase):
    def test_shared_load_provider_without_a_flight_scenario(self):
        vehicle=build_falcon9()
        for engine in vehicle.stages[0].engines:
            engine.running=True
            engine.throttle=1.0
        environment=TroposphereEnvironment()
        curve=MachDragCurve.from_json(Path(__file__).resolve().parents[1]/'data'/'drag_curve_pezent.json',21.08)
        provider=SteadyPoweredAscent(vehicle,environment,DragOnlyAerodynamics(21.08,curve))
        state=vehicle.to_state(consumption_by_stage={'stage1':provider.consumption_kgps[0]})
        loads=provider(0,state)
        np.testing.assert_allclose(loads.force_body_n,[9*845099,0,0],atol=1e-8)
        np.testing.assert_allclose(loads.consumption_kgps,provider.consumption_kgps)
        self.assertEqual(vehicle.stages[0].propellant_mass,410900)

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

    def test_environment_rejects_unmodeled_altitudes(self):
        with self.assertRaises(ValueError):TroposphereEnvironment.atmosphere(120001)
        with self.assertRaises(ValueError):TroposphereEnvironment.atmosphere(-1)

    def test_layered_atmosphere_checkpoint_values(self):
        # USSA 1976 geopotential-height checkpoints above the troposphere.
        checkpoints={
            11000:(216.65,22632.06,0.36392),
            20000:(216.65,5474.89,0.088035),
            32000:(228.65,868.02,0.013225),
            47000:(270.65,110.91,0.001428),
            51000:(270.65,66.939,0.000862),
        }
        for altitude,(temperature,pressure,density) in checkpoints.items():
            sample=TroposphereEnvironment.atmosphere(altitude)
            self.assertAlmostEqual(sample.temperature_k,temperature,delta=.01)
            self.assertAlmostEqual(sample.pressure_pa,pressure,delta=2e-3*pressure)
            self.assertAlmostEqual(sample.density_kgm3,density,delta=3e-3*density)

    def test_atmosphere_is_continuous_at_layer_breakpoints(self):
        for altitude in (11000,20000,32000,47000,51000,71000,84852):
            below=TroposphereEnvironment.atmosphere(altitude-.5)
            above=TroposphereEnvironment.atmosphere(altitude+.5)
            self.assertAlmostEqual(below.temperature_k,above.temperature_k,delta=.01)
            self.assertAlmostEqual(below.pressure_pa,above.pressure_pa,delta=1e-3*below.pressure_pa)

    def test_gravity_changes_body_components_when_attitude_turns(self):
        state=State([0,0,0],[0,0,0],[np.sqrt(.5),0,0,np.sqrt(.5)],[0,0,0],[0,0])
        np.testing.assert_allclose(TroposphereEnvironment().gravity_body(state),[0,9.80665,0],atol=1e-12)

    def test_inverse_square_gravity_falloff_option(self):
        radius=6371008.8
        env=TroposphereEnvironment(gravity_mps2=9.80665,gravity_earth_radius_m=radius)
        at_ground=State([0,0,0],[0,0,0],[1,0,0,0],[0,0,0],[0,0])
        np.testing.assert_allclose(env.gravity_body(at_ground),[-9.80665,0,0],atol=1e-12)
        at_radius=State([radius,0,0],[0,0,0],[1,0,0,0],[0,0,0],[0,0])
        np.testing.assert_allclose(env.gravity_body(at_radius),[-9.80665/4,0,0],atol=1e-12)
        with self.assertRaises(ValueError):
            TroposphereEnvironment(gravity_earth_radius_m=-1.0)

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



if __name__=='__main__':unittest.main()
