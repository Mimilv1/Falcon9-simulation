import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Environement import (Environment,Atmosphere,MESOPAUSE_ALTITUDE_M,
                          STANDARD_LAYERS,EARTH_MEAN_RADIUS_M,STANDARD_GRAVITY_MPS2)
from frames import (LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,body_to_ecef_matrix,
                    launch_position_body_m,normalize_quaternion,quaternion_product)
from dynamics import State

R=EARTH_MEAN_RADIUS_M
G0=STANDARD_GRAVITY_MPS2


def state_at(position_body,quaternion=None,velocity=(0,0,0)):
    return State(np.asarray(position_body,dtype=float),velocity,
                 quaternion if quaternion is not None else LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,
                 [0,0,0],[0,0])


def state_at_ecef(position_ecef,quaternion=None):
    quaternion=quaternion if quaternion is not None else LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ
    return state_at(body_to_ecef_matrix(quaternion).T@np.asarray(position_ecef,dtype=float),quaternion)


class AtmosphereChecks(unittest.TestCase):
    def test_sea_level_reference_point(self):
        a=Environment.atmosphere(0)
        self.assertAlmostEqual(a.temperature_k,288.15,places=9)
        self.assertAlmostEqual(a.pressure_pa,101325.0,places=6)
        self.assertAlmostEqual(a.density_kgm3,1.2250,places=6)
        self.assertAlmostEqual(a.sound_speed_mps,340.294,places=3)

    def test_published_ussa_1976_checkpoint_values(self):


        expected=((5000,255.65,54019.9,7.3612e-1),
                  (11000,216.65,22632.1,3.6392e-1),
                  (20000,216.65,5474.9,8.8035e-2),
                  (32000,228.65,868.02,1.3225e-2),
                  (47000,270.65,110.906,1.4275e-3),
                  (51000,270.65,66.939,8.6160e-4),
                  (71000,214.65,3.9564,6.4211e-5),
                  (80000,196.65,0.88627,1.5700e-5))
        for altitude,temperature,pressure,density in expected:
            a=Environment.atmosphere(altitude)
            self.assertAlmostEqual(a.temperature_k,temperature,places=2,msg=f'T at {altitude} m')


            np.testing.assert_allclose(a.pressure_pa,pressure,rtol=5e-5,err_msg=f'p at {altitude} m')
            np.testing.assert_allclose(a.density_kgm3,density,rtol=5e-5,err_msg=f'rho at {altitude} m')
            self.assertAlmostEqual(a.sound_speed_mps,
                                   np.sqrt(1.4*287.05287*a.temperature_k),places=6)

    def test_mesopause_and_tail_start(self):
        a=Environment.atmosphere(MESOPAUSE_ALTITUDE_M)
        self.assertAlmostEqual(a.temperature_k,186.946,places=3)
        self.assertAlmostEqual(a.pressure_pa,0.37338,places=5)

        self.assertAlmostEqual(Environment.atmosphere(MESOPAUSE_ALTITUDE_M+1).temperature_k,
                               a.temperature_k,places=12)
        self.assertAlmostEqual(Environment.atmosphere(2000000).temperature_k,
                               a.temperature_k,places=12)

    def test_tail_decays_exponentially_with_scale_height(self):


        g0,gas,tail_temperature=G0,287.05287,186.946
        ratio=Environment.atmosphere(95000).pressure_pa/Environment.atmosphere(90000).pressure_pa
        self.assertAlmostEqual(ratio,np.exp(-g0*5000/(gas*tail_temperature)),places=12)

    def test_tail_reaches_the_ceiling_as_numerical_vacuum(self):
        for altitude in (200000,1000000,2000000):
            sample=Environment.atmosphere(altitude)
            self.assertLess(sample.density_kgm3,1e-12,msg=f'density at {altitude} m')
        self.assertEqual(Environment.atmosphere(2000000).density_kgm3,
                         Environment.atmosphere(2000000).pressure_pa
                         /(287.05287*186.946))

    def test_rejects_altitudes_outside_the_modeled_envelope(self):
        for bad in (-1,2000001,float('nan'),float('inf'),-5e-5,-2e6):
            with self.assertRaises(ValueError,msg=f'altitude {bad}'):
                Environment.atmosphere(bad)

        for good in (0,2000000,-5e-7):
            Environment.atmosphere(good)

    def test_layer_boundaries_are_continuous(self):
        for boundary in [layer[0] for layer in STANDARD_LAYERS]+[MESOPAUSE_ALTITUDE_M]:
            below=Environment.atmosphere(boundary-1e-6)
            above=Environment.atmosphere(boundary+1e-6)
            np.testing.assert_allclose(below.pressure_pa,above.pressure_pa,rtol=1e-9)
            np.testing.assert_allclose(below.temperature_k,above.temperature_k,rtol=1e-9)
            np.testing.assert_allclose(below.density_kgm3,above.density_kgm3,rtol=1e-9)

    def test_pressure_and_density_decrease_monotonically(self):
        grid=np.unique(np.r_[np.linspace(0,MESOPAUSE_ALTITUDE_M,300),
                              np.linspace(MESOPAUSE_ALTITUDE_M,2000000,200)])
        for i in range(1,len(grid)):
            lower=Environment.atmosphere(grid[i-1])
            upper=Environment.atmosphere(grid[i])
            self.assertLess(upper.pressure_pa,lower.pressure_pa,
                            msg=f'pressure not decreasing at {grid[i]} m')
            self.assertLess(upper.density_kgm3,lower.density_kgm3,
                            msg=f'density not decreasing at {grid[i]} m')

    def test_hydrostatic_balance_dp_dh_equals_minus_g_rho(self):

        for altitude in (3000.0,15000.0,40000.0):
            dh=1.0
            mid=Environment.atmosphere(altitude)
            slope=(Environment.atmosphere(altitude+dh).pressure_pa
                   -Environment.atmosphere(altitude-dh).pressure_pa)/(2*dh)
            self.assertAlmostEqual(slope,-G0*mid.density_kgm3,places=4,
                                   msg=f'hydrostatic balance at {altitude} m')

    def test_outputs_are_finite_and_positive_everywhere(self):
        for altitude in np.linspace(0,2000000,421):
            a=Environment.atmosphere(altitude)
            for value in (a.temperature_k,a.pressure_pa,a.density_kgm3,a.sound_speed_mps):
                self.assertTrue(np.isfinite(value) and value>0,msg=f'at {altitude} m')

    def test_evaluator_is_pure_and_static(self):
        first=Environment.atmosphere(12345)
        second=Environment.atmosphere(12345)
        self.assertEqual(first,second)
        self.assertIsInstance(first,Atmosphere)


class GravityChecks(unittest.TestCase):
    def test_points_down_along_minus_body_x_at_the_launch_site(self):
        env=Environment()
        pad=state_at(launch_position_body_m(R))
        np.testing.assert_allclose(env.gravity_body(pad),[-G0,0,0],atol=1e-12)

    def test_respects_configured_surface_gravity(self):
        env=Environment(surface_gravity_mps2=3.7,earth_radius_m=1e6)
        pad=state_at(launch_position_body_m(1e6))
        np.testing.assert_allclose(env.gravity_body(pad),[-3.7,0,0],atol=1e-12)

    def test_inverse_square_falloff_with_distance_from_the_centre(self):
        env=Environment()
        near=state_at(launch_position_body_m(R,R))
        far=state_at(launch_position_body_m(R,2*R))
        np.testing.assert_allclose(env.gravity_body(near),[-G0/4,0,0],atol=1e-12)
        np.testing.assert_allclose(env.gravity_body(far),[-G0/9,0,0],atol=1e-12)

    def test_magnitude_consistent_with_altitude_above_the_site(self):
        env=Environment()
        height=500000.0
        state=state_at(launch_position_body_m(R,height))
        expected=G0*(R/(R+height))**2
        self.assertAlmostEqual(np.linalg.norm(env.gravity_body(state)),expected,places=9)

    def test_direction_and_magnitude_survive_the_body_round_trip(self):

        env=Environment()
        position_ecef=np.array([5e6,-3e6,7e6])
        quaternion=normalize_quaternion([0.6,0.3,-0.2,0.7])
        vector=env.gravity_body(state_at_ecef(position_ecef,quaternion))
        rebuilt=body_to_ecef_matrix(quaternion)@vector
        distance=np.linalg.norm(position_ecef)
        expected=-G0*R*R/distance**2*position_ecef/distance
        np.testing.assert_allclose(rebuilt,expected,rtol=1e-12)

    def test_tilts_toward_the_centre_downrange(self):
        env=Environment()
        downrange_ecef=np.array([R,0.0,R])
        vector=env.gravity_body(state_at_ecef(downrange_ecef))
        expected=-G0/2*downrange_ecef/np.linalg.norm(downrange_ecef)
        np.testing.assert_allclose(body_to_ecef_matrix(
            LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ)@vector,expected,atol=1e-12)

    def test_components_rotate_with_attitude_at_a_fixed_position(self):
        env=Environment(wind_ecef_mps=(0,0,0))
        yawed=quaternion_product(LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,[np.sqrt(.5),0,0,np.sqrt(.5)])
        site_body=body_to_ecef_matrix(yawed).T@np.array([0,0,R])
        vector=env.gravity_body(state_at(site_body,yawed))
        np.testing.assert_allclose(vector,[0,G0,0],atol=1e-12)

    def test_rejects_the_planet_centre(self):
        env=Environment()
        with self.assertRaises(ValueError):
            env.gravity_body(state_at([0,0,0]))


class WindChecks(unittest.TestCase):
    def test_default_is_calm(self):
        env=Environment()
        np.testing.assert_array_equal(env.wind_ecef_mps,[0,0,0])
        np.testing.assert_array_equal(env.wind_body(state_at(launch_position_body_m(R))),[0,0,0])

    def test_pad_wind_components_follow_the_launch_frame_mapping(self):


        env=Environment(wind_ecef_mps=(3,-4,0))
        np.testing.assert_allclose(env.wind_body(state_at(launch_position_body_m(R))),
                                   [0,-4,-3],atol=1e-12)
        env=Environment(wind_ecef_mps=(0,0,5))
        np.testing.assert_allclose(env.wind_body(state_at(launch_position_body_m(R))),
                                   [5,0,0],atol=1e-12)

    def test_wind_rotates_with_attitude_and_keeps_its_magnitude(self):
        env=Environment(wind_ecef_mps=(3,-4,0))
        yawed=quaternion_product(LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,[np.sqrt(.5),0,0,np.sqrt(.5)])
        pad=env.wind_body(state_at(launch_position_body_m(R)))
        turned=env.wind_body(state_at(launch_position_body_m(R),yawed))
        self.assertFalse(np.allclose(pad,turned))
        self.assertAlmostEqual(np.linalg.norm(pad),np.linalg.norm(turned),places=12)

        for quaternion in (LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,yawed):
            vector=env.wind_body(state_at(launch_position_body_m(R),quaternion))
            np.testing.assert_allclose(body_to_ecef_matrix(quaternion)@vector,[3,-4,0],atol=1e-12)

    def test_wind_ignores_position_velocity_and_time(self):
        env=Environment(wind_ecef_mps=(7,1,-2))
        pad=env.wind_body(state_at(launch_position_body_m(R)))
        moving=env.wind_body(state_at(launch_position_body_m(R),velocity=(100,20,-5)))
        coasting=env.wind_body(state_at(launch_position_body_m(R,500000)))
        np.testing.assert_array_equal(pad,moving)
        np.testing.assert_array_equal(pad,coasting)

    def test_constructor_stores_an_independent_copy(self):
        source=np.array([1,2,3])
        env=Environment(wind_ecef_mps=source)
        source[:]=0
        np.testing.assert_array_equal(env.wind_ecef_mps,[1,2,3])


class AltitudeChecks(unittest.TestCase):
    def test_zero_at_the_launch_site(self):
        env=Environment()
        self.assertAlmostEqual(env.altitude(state_at(launch_position_body_m(R)),None),0.0,delta=1e-6)

    def test_ground_altitude_shifts_the_reference_sphere(self):
        env=Environment(ground_altitude_m=100.0)
        self.assertAlmostEqual(env.altitude(state_at(launch_position_body_m(R)),None),100.0,delta=1e-6)
        self.assertAlmostEqual(env.altitude(state_at(launch_position_body_m(R,-100)),None),0.0,delta=1e-6)

    def test_measures_height_above_the_curved_surface(self):
        env=Environment()
        height=500.0
        for position_ecef in ((0,0,R+height),
                              (np.array([R+height,0,R+height])/np.sqrt(2)),
                              (0,R+height,0)):
            state=state_at_ecef(position_ecef)
            self.assertAlmostEqual(env.altitude(state,None),height,places=6,
                                   msg=f'altitude at {position_ecef}')
        downrange_ecef=np.array([R,0.0,R])
        state=state_at_ecef(downrange_ecef)
        self.assertAlmostEqual(env.altitude(state,None),R*(np.sqrt(2)-1),places=6)
        below=state_at(launch_position_body_m(R,-50))
        self.assertAlmostEqual(env.altitude(below,None),-50.0,delta=1e-6)

    def test_independent_of_attitude(self):
        env=Environment()
        position_ecef=np.array([2e6,-4e6,R+1234.0])
        yawed=quaternion_product(LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,[np.sqrt(.5),0,0,np.sqrt(.5)])
        readings=[env.altitude(state_at_ecef(position_ecef,q),None)
                   for q in (LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ,yawed,
                             normalize_quaternion([0.6,0.3,-0.2,0.7]))]
        for reading in readings[1:]:
            self.assertAlmostEqual(reading,readings[0],places=6)

    def test_station_reading_is_offset_from_the_cg(self):
        env=Environment()
        pad=state_at(launch_position_body_m(R))
        self.assertAlmostEqual(env.altitude(pad,SimpleNamespace(cg_station_body_m=np.zeros(3)),[10,0,0]),
                               10.0,delta=1e-6)

        self.assertAlmostEqual(env.altitude(pad,SimpleNamespace(cg_station_body_m=np.array([4.,0,0])),[10,0,0]),
                               6.0,delta=1e-6)

        lateral=env.altitude(pad,SimpleNamespace(cg_station_body_m=np.zeros(3)),[0,0,10])
        self.assertAlmostEqual(lateral,np.sqrt(R*R+100)-R,delta=1e-6)

        self.assertAlmostEqual(env.altitude(pad,None),0.0,delta=1e-6)


class ConstructionChecks(unittest.TestCase):
    def test_default_geometry_and_point_mass_parameter(self):
        env=Environment()
        np.testing.assert_array_equal(env.planet_center_ecef_m,[0,0,0])
        np.testing.assert_allclose(env.launch_site_ecef_m,[0,0,R])
        self.assertEqual(env.earth_radius_m,R)
        self.assertAlmostEqual(env.gravity_mu_m3ps2,G0*R*R,places=3)

    def test_rejects_invalid_scenarios(self):
        cases=[dict(earth_radius_m=-1),dict(earth_radius_m=0),
               dict(surface_gravity_mps2=0.0),dict(surface_gravity_mps2=float('nan')),
               dict(wind_ecef_mps=(1,float('nan'),0)),dict(wind_ecef_mps=(float('inf'),0,0)),
               dict(ground_altitude_m=float('nan')),dict(earth_radius_m=float('nan'))]
        for case in cases:
            with self.assertRaises(ValueError,msg=f'{case}'):
                Environment(**case)

    def test_parameters_are_keyword_only(self):
        with self.assertRaises(TypeError):
            Environment(9.80665)

    def test_wind_vector_must_have_three_finite_components(self):
        with self.assertRaises(ValueError):
            Environment(wind_ecef_mps=(1,2))
        with self.assertRaises(ValueError):
            Environment(wind_ecef_mps='wind')


class SimulatorSanityChecks(unittest.TestCase):

    def test_max_q_regime_produces_kilopascal_dynamic_pressure(self):

        air=Environment.atmosphere(12000)
        q=.5*air.density_kgm3*250**2
        self.assertGreater(q,5e3)
        self.assertLess(q,3e4)

    def test_gravity_and_atmosphere_cover_the_boostback_arc(self):


        env=Environment()
        apogee=state_at(launch_position_body_m(R,150000))
        self.assertAlmostEqual(np.linalg.norm(env.gravity_body(apogee)),
                                G0*(R/(R+150000))**2,places=6)
        self.assertGreater(Environment.atmosphere(120000).density_kgm3,1e-9)
        self.assertLess(Environment.atmosphere(150000).density_kgm3,1e-9)

    def test_drag_inputs_are_self_consistent_at_one_flight_point(self):


        env=Environment(wind_ecef_mps=(-12,3,0))
        quaternion=normalize_quaternion([0.9,0.1,0.2,0.3])
        position_ecef=np.array([3e5,7e5,R+21000])
        state=state_at_ecef(position_ecef,quaternion)
        air=Environment.atmosphere(env.altitude(state,None))

        self.assertAlmostEqual(env.altitude(state,None),
                               np.linalg.norm(position_ecef)-R,places=6)
        self.assertAlmostEqual(np.linalg.norm(env.gravity_body(state)),
                               G0*R*R/np.linalg.norm(position_ecef)**2,places=8)
        self.assertAlmostEqual(np.linalg.norm(env.wind_body(state)),
                               np.linalg.norm([-12,3,0]),places=9)
        self.assertGreater(air.sound_speed_mps,0)


if __name__=='__main__':
    unittest.main()
