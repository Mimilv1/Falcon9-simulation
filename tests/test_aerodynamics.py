"""Multi-configuration aerodynamics tests: force/moment conventions, retro
flow, damping, configuration data and provider integration."""
import unittest
from pathlib import Path
import numpy as np

from Aerodynamics import (Aerodynamics,AeroConfiguration,SideSegment,GridFinControl,
                           DragOnlyAerodynamics,MachDragCurve)
from Environement import Environment
from dynamics import State
from mass_properties import MassProperties
from ascent import SteadyPoweredAscent
from VehicleData import build_falcon9

ATMOSPHERE=Environment.atmosphere(0.0)
DATA=Path(__file__).resolve().parents[1]/'data'
CONFIGURATIONS=DATA/'aero_configurations.json'
DRAG_CURVE=DATA/'drag_curve_pezent.json'


def properties(cg_x_m=21.0):
    return MassProperties(mass_kg=25000.0,cg_station_body_m=np.array([cg_x_m,0.0,0.0]),
                          inertia_body_kgm2=np.eye(3),inertia_rate_body_kgm2ps=np.zeros((3,3)),
                          cg_rate_body_mps=np.zeros(3))


def state(velocity_body,omega=(0.0,0.0,0.0)):
    return State(np.array([7000000.0,0.0,0.0]),np.asarray(velocity_body,dtype=float),
                 np.array([1.0,0.0,0.0,0.0]),np.asarray(omega,dtype=float),[400000.0,100000.0])


def booster_configuration(name='booster',crossflow=1.1,pitch=-20.0,roll=-2.0,
                         center=38.6,slope=2.0):
    return AeroConfiguration(name=name,reference_area_m2=10.52,reference_length_m=42.6,
        axial_drag=MachDragCurve([0.0,1.0,2.0],[0.42,0.94,0.6]),
        potential_normal_slope_per_rad=slope,potential_center_station_m=center,
        side_segments=(SideSegment(0.0,42.6,3.66,crossflow),),
        pitch_damping_coefficient=pitch,roll_damping_coefficient=roll)


class TestAxialEquivalence(unittest.TestCase):
    def test_zero_alpha_matches_drag_only(self):
        curve=MachDragCurve([0.0,1.0,2.0],[0.42,0.94,0.6])
        model=Aerodynamics({'booster':booster_configuration()})
        reference=DragOnlyAerodynamics(10.52,curve)
        for speed in (240.0,600.0):
            result=model.evaluate(state([speed,0.0,0.0]),properties(),ATMOSPHERE,
                                  np.zeros(3),'booster')
            expected=reference.evaluate(state([speed,0.0,0.0]),properties(),ATMOSPHERE,
                                         np.zeros(3))
            np.testing.assert_allclose(result.force_body_n,expected.force_body_n)
            self.assertEqual(result.drag_coefficient,expected.drag_coefficient)
            self.assertEqual(result.angle_of_attack_rad,0.0)
            np.testing.assert_allclose(result.free_moment_body_nm,np.zeros(3),atol=1e-12)


class TestLateralForceAndMoment(unittest.TestCase):
    def test_normal_force_opposes_lateral_relative_velocity(self):
        model=Aerodynamics({'booster':booster_configuration()})
        result=model.evaluate(state([200.0,30.0,-40.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'booster')
        self.assertGreater(result.normal_force_coefficient,0.0)
        self.assertLess(np.dot(result.force_body_n[1:],np.array([30.0,-40.0])),0.0)
        self.assertGreater(result.angle_of_attack_rad,0.0)

    def test_station_moments_are_transverse_and_force_free(self):
        model=Aerodynamics({'booster':booster_configuration()})
        result=model.evaluate(state([200.0,-15.0,25.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'booster')
        moment=result.free_moment_body_nm
        self.assertAlmostEqual(moment[0],0.0,places=9)
        tolerance=1e-7*np.linalg.norm(moment)*np.linalg.norm(result.force_body_n)+1e-9
        self.assertLess(abs(float(np.dot(moment,result.force_body_n))),tolerance)

    def test_crossflow_only_at_ninety_degrees(self):
        body=booster_configuration(name='body',crossflow=1.1)
        with_fins=AeroConfiguration(name='fins',reference_area_m2=10.52,
            reference_length_m=42.6,axial_drag=body.axial_drag,
            potential_normal_slope_per_rad=2.0,potential_center_station_m=38.6,
            side_segments=(SideSegment(0.0,42.6,3.66,1.1),
                           SideSegment(38.6,40.6,2.4,1.0)),
            pitch_damping_coefficient=-20.0)
        model=Aerodynamics({'body':body,'fins':with_fins})
        flight=state([0.0,150.0,0.0])
        plain=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),'body')
        finned=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),'fins')
        expected=1.1*3.66*42.6/10.52
        self.assertAlmostEqual(plain.normal_force_coefficient,expected,places=9)
        self.assertAlmostEqual(plain.angle_of_attack_rad,np.pi/2,places=9)
        self.assertGreater(finned.normal_force_coefficient,
                           plain.normal_force_coefficient)
        self.assertTrue(np.any(np.not_equal(finned.free_moment_body_nm,
                                             plain.free_moment_body_nm)))


class TestRetroFlow(unittest.TestCase):
    def test_engine_first_axial_drag_opposes_motion(self):
        model=Aerodynamics({'booster':booster_configuration()})
        result=model.evaluate(state([-200.0,0.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'booster')
        self.assertGreater(result.force_body_n[0],0.0)
        self.assertLess(float(np.dot(result.force_body_n,np.array([-200.0,0.0,0.0]))),0.0)
        self.assertAlmostEqual(result.angle_of_attack_rad,np.pi,places=12)

    def test_retro_potential_center_is_mirrored_about_cg(self):
        configuration=booster_configuration(crossflow=0.0)
        model=Aerodynamics({'booster':configuration})
        relative=np.array([-200.0,0.0,12.0])
        result=model.evaluate(state(relative),properties(),ATMOSPHERE,
                              np.zeros(3),'booster')
        speed=float(np.linalg.norm(relative))
        direction=relative/speed
        cosine=float(direction[0]); sine=float(np.linalg.norm(direction[1:]))
        dynamic_pressure=.5*ATMOSPHERE.density_kgm3*speed*speed
        magnitude=dynamic_pressure*10.52*2.0*sine*abs(cosine)
        side=np.r_[0.0,-direction[1:]/sine]
        expected=np.cross(np.array([2.0*21.0-38.6-21.0,0.0,0.0]),magnitude*side)
        np.testing.assert_allclose(result.free_moment_body_nm,expected,rtol=1e-12)
        # The mirrored flow with the mirrored station yields the same lateral
        # moment: the mirror maps arm and force direction coherently.
        nose_direction=-relative/speed
        nose_side=np.r_[0.0,-nose_direction[1:]/sine]
        nose_expected=np.cross(np.array([38.6-21.0,0.0,0.0]),magnitude*nose_side)
        nose_first=model.evaluate(state(-relative),properties(),ATMOSPHERE,
                                   np.zeros(3),'booster')
        np.testing.assert_allclose(nose_first.free_moment_body_nm,nose_expected,rtol=1e-12)
        np.testing.assert_allclose(nose_first.free_moment_body_nm,
                                   result.free_moment_body_nm,rtol=1e-12)


class TestRotaryDamping(unittest.TestCase):
    def test_damping_follows_declared_coefficients(self):
        configuration=AeroConfiguration(name='axial',reference_area_m2=10.52,
            reference_length_m=42.6,axial_drag=MachDragCurve([0.0,1.0,2.0],[0.42,0.94,0.6]),
            potential_normal_slope_per_rad=0.0,potential_center_station_m=38.6,
            side_segments=(SideSegment(0.0,42.6,3.66,0.0),),
            pitch_damping_coefficient=-20.0,roll_damping_coefficient=-2.0)
        model=Aerodynamics({'axial':configuration})
        flight=state([200.0,0.0,0.0],omega=[0.1,0.5,-0.2])
        result=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),'axial')
        speed=200.0
        dynamic_pressure=.5*ATMOSPHERE.density_kgm3*speed*speed
        scale=dynamic_pressure*10.52*42.6**2/(2.0*speed)
        expected=scale*np.array([-2.0*0.1,-20.0*0.5,-20.0*(-0.2)])
        # With slope and crossflow zero, the entire free moment is damping.
        np.testing.assert_allclose(result.free_moment_body_nm,expected,rtol=1e-12)

    def test_damping_opposes_rotation_signs(self):
        model=Aerodynamics({'booster':booster_configuration()})
        result=model.evaluate(state([200.0,0.0,0.0],omega=[0.0,0.4,0.0]),
                              properties(),ATMOSPHERE,np.zeros(3),'booster')
        self.assertLess(result.free_moment_body_nm[1],0.0)


class TestDomainAndRest(unittest.TestCase):
    def test_rest_returns_zero_loads_and_zero_moment(self):
        model=Aerodynamics({'booster':booster_configuration()})
        result=model.evaluate(state([0.0,0.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'booster')
        np.testing.assert_allclose(result.force_body_n,np.zeros(3))
        np.testing.assert_allclose(result.free_moment_body_nm,np.zeros(3))
        self.assertEqual(result.drag_coefficient,0.42)

    def test_unmodeled_mach_is_rejected(self):
        model=Aerodynamics({'booster':booster_configuration()})
        with self.assertRaises(ValueError):
            model.evaluate(state([3000.0,0.0,0.0]),properties(),ATMOSPHERE,
                           np.zeros(3),'booster')

    def test_unknown_configuration_is_rejected(self):
        model=Aerodynamics({'booster':booster_configuration()})
        with self.assertRaises(ValueError):
            model.evaluate(state([200.0,0.0,0.0]),properties(),ATMOSPHERE,
                           np.zeros(3),'not_a_configuration')
        with self.assertRaises(ValueError):
            Aerodynamics({'bad':AeroConfiguration(name='mismatch',reference_area_m2=10.52,
                reference_length_m=42.6,axial_drag=MachDragCurve([0.0,1.0],[0.4,0.4]),
                potential_normal_slope_per_rad=2.0,potential_center_station_m=38.6,
                side_segments=(SideSegment(0.0,42.6,3.66,1.1),),
                pitch_damping_coefficient=-20.0)})


class TestConfigurationData(unittest.TestCase):
    def test_declared_configuration_file_loads_all_phases(self):
        model=Aerodynamics.from_json(CONFIGURATIONS)
        self.assertEqual(set(model.configurations),
                         {'full_stack','booster_alone','booster_landing','second_stack'})
        for name in model.configurations:
            result=model.evaluate(state([200.0,10.0,-5.0]),properties(),ATMOSPHERE,
                                  np.zeros(3),name)
            self.assertGreater(result.dynamic_pressure_pa,0.0)
        rest=model.evaluate(state([0.0,0.0,0.0]),properties(),ATMOSPHERE,
                            np.zeros(3),'booster_landing')
        np.testing.assert_allclose(rest.free_moment_body_nm,np.zeros(3))

    def test_full_stack_axial_curve_is_the_sourced_evidence(self):
        model=Aerodynamics.from_json(CONFIGURATIONS)
        expected=MachDragCurve.from_json(DRAG_CURVE,21.08)
        self.assertEqual(model.configurations['full_stack'].reference_area_m2,21.08)
        self.assertAlmostEqual(model.configurations['full_stack'].axial_drag(1.5),
                               expected(1.5),places=12)

    def test_landing_configuration_has_incremental_fin_and_leg_segments(self):
        configuration=Aerodynamics.from_json(CONFIGURATIONS).configurations['booster_landing']
        self.assertEqual(len(configuration.side_segments),3)
        self.assertGreater(sum(s.side_area_m2 for s in configuration.side_segments),
                            Aerodynamics.from_json(CONFIGURATIONS)
                            .configurations['booster_alone'].side_segments[0].side_area_m2)


class TestPhysicalInvariants(unittest.TestCase):
    """Model-level invariants the propagation core and load assembly rely on."""

    def test_aerodynamic_power_is_never_positive(self):
        """Every force term must extract energy from the relative motion."""
        model=Aerodynamics({'booster':booster_configuration()})
        for speed in (80.0,200.0,670.0):
            for angle in np.linspace(0.0,np.pi,9):
                direction=np.array([np.cos(angle),np.sin(angle)*np.cos(0.7),
                                    np.sin(angle)*np.sin(0.7)])
                relative=speed*direction
                result=model.evaluate(state(relative),properties(),ATMOSPHERE,
                                      np.zeros(3),'booster')
                tolerance=1e-9*abs(np.dot(relative,result.force_body_n)) \
                          +1e-6*np.linalg.norm(relative)*np.linalg.norm(result.force_body_n)
                self.assertLessEqual(float(np.dot(relative,result.force_body_n)),
                                     tolerance)

    def test_damping_power_is_never_positive(self):
        """Rotary damping with negative coefficients must oppose rotation."""
        model=Aerodynamics({'booster':booster_configuration(roll=-2.0)})
        omega=[0.3,-0.4,0.5]
        result=model.evaluate(state([300.0,0.0,0.0],omega=omega),
                              properties(),ATMOSPHERE,np.zeros(3),'booster')
        self.assertLess(float(np.dot(result.free_moment_body_nm,omega)),0.0)

    def test_wind_enters_only_through_relative_velocity(self):
        model=Aerodynamics({'booster':booster_configuration()})
        wind=np.array([3.0,-2.0,5.0])
        flight=state([260.0,-20.0,15.0])
        with_wind=model.evaluate(flight,properties(),ATMOSPHERE,wind,'booster')
        shifted=model.evaluate(state([257.0,-18.0,10.0]),properties(),ATMOSPHERE,
                               np.zeros(3),'booster')
        np.testing.assert_allclose(with_wind.force_body_n,shifted.force_body_n)
        np.testing.assert_allclose(with_wind.free_moment_body_nm,shifted.free_moment_body_nm)
        self.assertAlmostEqual(with_wind.relative_speed_mps,shifted.relative_speed_mps)
        self.assertAlmostEqual(with_wind.mach,shifted.mach)

    def test_force_field_is_continuous_across_ninety_degrees(self):
        model=Aerodynamics({'booster':booster_configuration()})
        import math
        half=1e-6
        before=model.evaluate(state([half,150.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'booster')
        after=model.evaluate(state([-half,150.0,0.0]),properties(),ATMOSPHERE,
                             np.zeros(3),'booster')
        np.testing.assert_allclose(before.force_body_n,after.force_body_n,atol=1e-3)
        np.testing.assert_allclose(before.free_moment_body_nm,after.free_moment_body_nm,
                                   atol=1e-1)
        self.assertAlmostEqual(before.angle_of_attack_rad,np.pi/2,places=5)
        self.assertAlmostEqual(after.angle_of_attack_rad,np.pi/2,places=5)

    def test_outputs_scale_with_dynamic_pressure(self):
        model=Aerodynamics({'booster':booster_configuration()})
        slow=model.evaluate(state([170.0,12.0,0.0]),properties(),ATMOSPHERE,
                            np.zeros(3),'booster')
        fast=model.evaluate(state([340.0,24.0,0.0]),properties(),ATMOSPHERE,
                            np.zeros(3),'booster')
        self.assertAlmostEqual(fast.mach,2.0*slow.mach)
        self.assertAlmostEqual(fast.dynamic_pressure_pa,4.0*slow.dynamic_pressure_pa)
        self.assertAlmostEqual(fast.normal_force_coefficient,slow.normal_force_coefficient,
                               places=6)

    def test_reports_are_consistent_with_stated_formulas(self):
        """force and moment reproduce the documented decomposition term by term."""
        configuration=booster_configuration(crossflow=1.1)
        model=Aerodynamics({'booster':configuration})
        relative=np.array([250.0,0.0,30.0])
        cg=properties(21.0)
        result=model.evaluate(state(relative),cg,ATMOSPHERE,np.zeros(3),'booster')
        speed=float(np.linalg.norm(relative))
        direction=relative/speed
        cosine=float(direction[0]); sine=float(np.linalg.norm(direction[1:]))
        qbar=.5*ATMOSPHERE.density_kgm3*speed*speed
        axial=-qbar*10.52*configuration.axial_drag(speed/ATMOSPHERE.sound_speed_mps)\
              *cosine*abs(cosine)
        potential=qbar*10.52*configuration.potential_normal_slope_per_rad*sine*abs(cosine)
        segment=configuration.side_segments[0]
        crossflow=qbar*segment.crossflow_drag_coefficient*segment.side_area_m2*sine*sine
        side=np.r_[0.0,-direction[1:]/sine]
        expected_force=np.array([axial,0.0,0.0])+(potential+crossflow)*side
        expected_moment=(np.cross(np.array([38.6-21.0,0.0,0.0]),potential*side)
                         +np.cross(np.array([segment.centroid_m-21.0,0.0,0.0]),
                                   crossflow*side))
        np.testing.assert_allclose(result.force_body_n,expected_force,rtol=1e-12)
        np.testing.assert_allclose(result.free_moment_body_nm,expected_moment,rtol=1e-12)
        self.assertAlmostEqual(result.normal_force_coefficient,
                                (potential+crossflow)/(qbar*10.52),places=12)


class TestHypersonicExtension(unittest.TestCase):
    """Mach 7..10 knots are declared estimates; the domain ends there."""

    def setUp(self):
        self.model=Aerodynamics.from_json(CONFIGURATIONS)
        self.atmosphere=Environment.atmosphere(80000.0)

    def test_all_configurations_evaluate_across_mach_seven(self):
        for name in self.model.configurations:
            for mach in (7.01,8.0,9.99):
                speed=mach*self.atmosphere.sound_speed_mps
                result=self.model.evaluate(state([speed,0.0,0.0]),properties(),
                                          self.atmosphere,np.zeros(3),name)
                self.assertAlmostEqual(result.mach,mach,places=6)

    def test_extension_rejects_mach_beyond_ten(self):
        for name in self.model.configurations:
            with self.assertRaises(ValueError):
                self.model.evaluate(state([10.5*self.atmosphere.sound_speed_mps,0.0,0.0]),
                                    properties(),self.atmosphere,np.zeros(3),name)

    def test_full_stack_extension_joins_the_sourced_knot(self):
        model=self.model
        self.assertAlmostEqual(model.configurations['full_stack'].axial_drag(7.0),0.8131,
                               places=4)
        self.assertLess(model.configurations['full_stack'].axial_drag(10.0),
                        model.configurations['full_stack'].axial_drag(8.0))
        self.assertGreater(model.configurations['full_stack'].axial_drag(10.0),
                           0.7*model.configurations['full_stack'].axial_drag(7.0))

    def test_vertical_ascent_propagates_past_mach_seven(self):
        from propagation import propagate
        vehicle=build_falcon9()
        for engine in vehicle.stages[0].engines:
            engine.running=True
            engine.throttle=1.0
        provider=SteadyPoweredAscent(vehicle,Environment(ground_altitude_m=10.0,earth_rotation_rate_radps=0.0),
                                     self.model,aero_configuration='full_stack')
        initial=vehicle.to_state(consumption_by_stage={
            'stage1':provider.consumption_kgps[0]})
        result=propagate(0.0,initial,150.0,0.5,vehicle.mass_model(),None,provider)
        self.assertIn(result.stop_reason,('end_time','propellant_limit'))
        machs=[provider.evaluate(t,s).aero.mach
               for t,s in zip(result.times_s[::4],result.states[::4])]
        self.assertGreater(max(machs),7.0)
        self.assertLessEqual(max(machs),10.0+1e-9)


def fin_control_configuration():
    """Isolated fin control: no potential force, no crossflow, no damping."""
    return AeroConfiguration(name='fins',reference_area_m2=10.52,
        reference_length_m=42.6,axial_drag=MachDragCurve([0.0,1.0],[0.4,0.4]),
        potential_normal_slope_per_rad=0.0,potential_center_station_m=21.0,
        side_segments=(SideSegment(0.0,42.6,3.66,0.0),),
        pitch_damping_coefficient=0.0,roll_damping_coefficient=0.0,
        grid_fin_control=GridFinControl(fin_count=4,
            azimuths_rad=(0.0,np.pi/2,np.pi,3*np.pi/2),
            fin_area_m2=2.0,lift_slope_per_rad=1.0,pivot_radius_m=2.0,
            station_m=40.0,maximum_deflection_rad=0.3))


class TestGridFinControl(unittest.TestCase):
    """Fin increments: tangential force at the mount station, one moment each."""

    def test_single_fin_force_and_moment_at_mount_station(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        result=model.evaluate(state([300.0,0.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'fins',[0.2,0.0,0.0,0.0])
        dynamic_pressure=.5*ATMOSPHERE.density_kgm3*300.0**2
        magnitude=dynamic_pressure*2.0*1.0*0.2
        expected_force=np.array([0.0,0.0,magnitude])
        expected_moment=np.cross(np.array([40.0,2.0,0.0])-properties().cg_station_body_m,
                                 expected_force)
        self.assertAlmostEqual(result.force_body_n[2],magnitude,places=6)
        np.testing.assert_allclose(result.free_moment_body_nm,expected_moment,rtol=1e-12)
        self.assertAlmostEqual(result.grid_fin_control_force_n,magnitude,places=6)

    def test_equal_deflections_give_pure_roll_moment(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        result=model.evaluate(state([300.0,0.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'fins',[0.2,0.2,0.2,0.2])
        single=.5*ATMOSPHERE.density_kgm3*300.0**2*2.0*0.2
        np.testing.assert_allclose(result.force_body_n[1:],np.zeros(2),atol=1e-8)
        np.testing.assert_allclose(result.free_moment_body_nm,
                                   np.array([4*2.0*single,0.0,0.0]),
                                   rtol=1e-9,atol=1e-6)

    def test_opposite_pair_differential_gives_pure_pitch(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        result=model.evaluate(state([300.0,0.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'fins',[0.2,0.0,-0.2,0.0])
        single=.5*ATMOSPHERE.density_kgm3*300.0**2*2.0*0.2
        # The axial aero force stays untouched; only the fin increment steers.
        self.assertAlmostEqual(result.force_body_n[2],2.0*single,places=6)
        self.assertAlmostEqual(result.force_body_n[1],0.0,places=8)
        np.testing.assert_allclose(result.free_moment_body_nm,
                                   np.array([0.0,-2.0*(40.0-21.0)*single,0.0]),
                                   rtol=1e-9,atol=1e-6)

    def test_control_force_scales_with_dynamic_pressure(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        commands=[0.2,0.0,0.0,0.0]
        slow=model.evaluate(state([150.0,0.0,0.0]),properties(),ATMOSPHERE,
                            np.zeros(3),'fins',commands)
        fast=model.evaluate(state([300.0,0.0,0.0]),properties(),ATMOSPHERE,
                            np.zeros(3),'fins',commands)
        self.assertAlmostEqual(fast.grid_fin_control_force_n,
                                4.0*slow.grid_fin_control_force_n)

    def test_single_fin_across_all_declared_azimuths(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        cg=properties().cg_station_body_m
        dynamic_pressure=.5*ATMOSPHERE.density_kgm3*300.0**2
        for index,azimuth in enumerate((0.0,np.pi/2,np.pi,3*np.pi/2)):
            commands=[0.0,0.0,0.0,0.0]
            commands[index]=0.2
            result=model.evaluate(state([300.0,0.0,0.0]),properties(),ATMOSPHERE,
                                  np.zeros(3),'fins',commands)
            fin_force=(np.array([0.0,-np.sin(azimuth),np.cos(azimuth)])
                                 *dynamic_pressure*2.0*1.0*0.2)
            fin_station=np.array([40.0,2.0*np.cos(azimuth),2.0*np.sin(azimuth)])
            np.testing.assert_allclose(
                result.free_moment_body_nm,
                np.cross(fin_station-cg,fin_force),rtol=1e-9,atol=1e-6)
            np.testing.assert_allclose(
                result.force_body_n[1:],fin_force[1:],rtol=1e-9,atol=1e-6)

    def test_fin_control_reverses_in_retro_flow(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        commands=[0.2,0.0,0.0,0.0]
        nose=model.evaluate(state([300.0,0.0,0.0]),properties(),ATMOSPHERE,
                            np.zeros(3),'fins',commands)
        retro=model.evaluate(state([-300.0,0.0,0.0]),properties(),ATMOSPHERE,
                             np.zeros(3),'fins',commands)
        # Same speed, mirrored through-flow: the fin increment flips sign.
        self.assertAlmostEqual(retro.force_body_n[2],-nose.force_body_n[2],
                               places=6)
        np.testing.assert_allclose(retro.free_moment_body_nm,
                                   -nose.free_moment_body_nm,rtol=1e-9,atol=1e-6)
        self.assertAlmostEqual(retro.grid_fin_control_force_n,
                                nose.grid_fin_control_force_n,places=3)

    def test_fin_control_follows_the_axial_flow_share(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        commands=[0.2,0.0,0.0,0.0]
        result45=model.evaluate(state([200.0,200.0,0.0]),properties(),ATMOSPHERE,
                                np.zeros(3),'fins',commands)
        cos45=np.cos(np.pi/4)
        dynamic_pressure=.5*ATMOSPHERE.density_kgm3*(200.0**2+200.0**2)
        expected=dynamic_pressure*2.0*1.0*0.2*cos45*abs(cos45)
        self.assertAlmostEqual(result45.force_body_n[2],expected,places=6)
        result90=model.evaluate(state([0.0,150.0,0.0]),properties(),ATMOSPHERE,
                                np.zeros(3),'fins',commands)
        self.assertEqual(result90.grid_fin_control_force_n,0.0)

    def test_control_adds_exactly_to_uncommanded_flight(self):
        model=Aerodynamics.from_json(CONFIGURATIONS)
        flight=state([300.0,-20.0,0.0])
        commanded=[0.15,0.0,-0.15,0.0]
        uncommanded=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                                   'booster_landing')
        steered=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                               'booster_landing',commanded)
        control=model.configurations['booster_landing'].grid_fin_control
        direction=flight.velocity_body_mps/np.linalg.norm(flight.velocity_body_mps)
        share=direction[0]*abs(direction[0])
        dynamic_pressure=.5*ATMOSPHERE.density_kgm3*np.dot(flight.velocity_body_mps,
                                                           flight.velocity_body_mps)
        increment=0.0
        expected_force=np.zeros(3); expected_moment=np.zeros(3)
        cg=properties().cg_station_body_m
        for azimuth,deflection in zip(control.azimuths_rad,commanded):
            if deflection==0.0:
                continue
            magnitude=dynamic_pressure*control.fin_area_m2\
                      *control.lift_slope_per_rad*deflection*share
            fin_force=np.array([0.0,-np.sin(azimuth),np.cos(azimuth)])*magnitude
            fin_station=np.array([control.station_m,
                                    control.pivot_radius_m*np.cos(azimuth),
                                    control.pivot_radius_m*np.sin(azimuth)])
            expected_force+=fin_force
            expected_moment+=np.cross(fin_station-cg,fin_force)
        np.testing.assert_allclose(steered.force_body_n-uncommanded.force_body_n,
                                   expected_force,rtol=1e-9,atol=1e-6)
        np.testing.assert_allclose(steered.free_moment_body_nm
                                   -uncommanded.free_moment_body_nm,
                                   expected_moment,rtol=1e-9,atol=1e-6)

    def test_boundary_deflection_is_accepted_and_beyond_rejected(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        flight=state([300.0,0.0,0.0])
        model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                       'fins',[0.3,0.0,0.0,0.0])
        model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                       'fins',[-0.3,0.0,0.0,0.0])
        with self.assertRaises(ValueError):
            model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                            'fins',[0.3000001,0.0,0.0,0.0])

    def test_rest_with_valid_commands_returns_zero_loads(self):
        model=Aerodynamics({'fins':fin_control_configuration()})
        result=model.evaluate(state([0.0,0.0,0.0]),properties(),ATMOSPHERE,
                              np.zeros(3),'fins',[0.2,0.2,0.2,0.2])
        np.testing.assert_allclose(result.force_body_n,np.zeros(3))
        np.testing.assert_allclose(result.free_moment_body_nm,np.zeros(3))
        self.assertEqual(result.grid_fin_control_force_n,0.0)

    def test_roll_command_accelerates_roll_in_the_dynamics_core(self):
        from mass_properties import MassProperties
        from load_interface import AppliedForce,assemble_loads
        from dynamics import accelerations
        model=Aerodynamics({'fins':fin_control_configuration()})
        flight=state([300.0,0.0,0.0])
        mass=MassProperties(mass_kg=25000.0,
                            cg_station_body_m=np.array([21.0,0.0,0.0]),
                            inertia_body_kgm2=np.diag([100000.0,1e6,1e6]),
                            inertia_rate_body_kgm2ps=np.zeros((3,3)),
                            cg_rate_body_mps=np.zeros(3))
        result=model.evaluate(flight,mass,ATMOSPHERE,np.zeros(3),
                               'fins',[0.2,0.2,0.2,0.2])
        loads=assemble_loads(mass,[AppliedForce(result.force_body_n,
                                                mass.cg_station_body_m)],
                            gravity_body_mps2=np.zeros(3),
                            free_moment_body_nm=result.free_moment_body_nm)
        _,angular=accelerations(flight,loads,mass)
        self.assertGreater(angular[0],0.0)
        self.assertAlmostEqual(angular[0],
                               result.free_moment_body_nm[0]/100000.0,places=6)
        self.assertAlmostEqual(angular[1],0.0,places=6)
        self.assertAlmostEqual(angular[2],0.0,places=6)

    def test_command_contract_is_enforced(self):
        model=Aerodynamics.from_json(CONFIGURATIONS)
        flight=state([300.0,0.0,0.0])
        with self.assertRaises(ValueError):
            model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                           'booster_landing',[0.1,0.1,0.1])
        with self.assertRaises(ValueError):
            model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                           'booster_landing',[0.5,0.1,0.1,0.1])
        with self.assertRaises(ValueError):
            model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                           'booster_landing',[float('nan'),0.0,0.0,0.0])
        with self.assertRaises(ValueError):
            model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                           'booster_alone',[0.1,0.1,0.1,0.1])
        model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                       'booster_alone',[0.0,0.0,0.0,0.0])

    def test_zero_commands_match_the_undeflected_model(self):
        model=Aerodynamics.from_json(CONFIGURATIONS)
        flight=state([200.0,15.0,0.0])
        undeflected=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                                   'booster_landing')
        zeros=model.evaluate(flight,properties(),ATMOSPHERE,np.zeros(3),
                             'booster_landing',[0.0,0.0,0.0,0.0])
        np.testing.assert_array_equal(undeflected.force_body_n,zeros.force_body_n)
        np.testing.assert_array_equal(undeflected.free_moment_body_nm,
                                      zeros.free_moment_body_nm)

    def test_only_booster_landing_declares_fin_control(self):
        model=Aerodynamics.from_json(CONFIGURATIONS)
        for name in ('full_stack','booster_alone','second_stack'):
            self.assertIsNone(model.configurations[name].grid_fin_control)
        control=model.configurations['booster_landing'].grid_fin_control
        self.assertEqual(control.fin_count,4)
        self.assertEqual(len(control.azimuths_rad),4)

    def test_invalid_fin_declarations_are_rejected(self):
        with self.assertRaises(ValueError):
            GridFinControl(fin_count=0,azimuths_rad=(),fin_area_m2=1.0,
                           lift_slope_per_rad=1.0,pivot_radius_m=1.0,
                           station_m=39.0,maximum_deflection_rad=0.3)
        with self.assertRaises(ValueError):
            AeroConfiguration(name='bad',reference_area_m2=10.52,reference_length_m=42.6,
                axial_drag=MachDragCurve([0.0,1.0],[0.4,0.4]),
                potential_normal_slope_per_rad=0.0,potential_center_station_m=21.0,
                side_segments=(SideSegment(0.0,42.6,3.66,1.1),),
                pitch_damping_coefficient=-20.0,
                grid_fin_control=GridFinControl(fin_count=1,azimuths_rad=(0.0,),
                    fin_area_m2=1.0,lift_slope_per_rad=1.0,pivot_radius_m=1.0,
                    station_m=50.0,maximum_deflection_rad=0.3))


class TestProviderIntegration(unittest.TestCase):
    def test_steady_provider_uses_selected_configuration_and_free_moment(self):
        vehicle=build_falcon9()
        for engine in vehicle.stages[0].engines:
            engine.running=True
            engine.throttle=1.0
        model=Aerodynamics.from_json(CONFIGURATIONS)
        provider=SteadyPoweredAscent(vehicle,Environment(),model,aero_configuration='full_stack')
        initial=vehicle.to_state(consumption_by_stage={'stage1':provider.consumption_kgps[0]})
        evaluation=provider.evaluate(0.0,initial)
        self.assertAlmostEqual(evaluation.aero.relative_speed_mps,0.0)
        self.assertGreater(evaluation.thrust_n,0.0)
        np.testing.assert_allclose(evaluation.aero.free_moment_body_nm,np.zeros(3))

    def test_steady_provider_rejects_unknown_configuration(self):
        vehicle=build_falcon9()
        model=Aerodynamics.from_json(CONFIGURATIONS)
        with self.assertRaises(ValueError):
            SteadyPoweredAscent(vehicle,Environment(),model,aero_configuration='stage_three')

    def test_vertical_propagation_runs_with_the_multiconfiguration_model(self):
        from propagation import propagate
        vehicle=build_falcon9()
        for engine in vehicle.stages[0].engines:
            engine.running=True
            engine.throttle=1.0
        model=Aerodynamics.from_json(CONFIGURATIONS)
        # Pad exit nozzles sit at the site level; the strict environment
        # requires the site above the ground datum (ground_altitude_m>0).
        provider=SteadyPoweredAscent(vehicle,Environment(ground_altitude_m=10.0,earth_rotation_rate_radps=0.0),model,
                                     aero_configuration='full_stack')
        initial=vehicle.to_state(consumption_by_stage={
            'stage1':provider.consumption_kgps[0]})
        result=propagate(0.0,initial,2.0,0.25,vehicle.mass_model(),None,provider)
        final=result.states[-1]
        self.assertEqual(result.stop_reason,'end_time')
        self.assertGreater(final.velocity_body_mps[0],5.0)
        self.assertGreater(final.position_body_m[0],initial.position_body_m[0]+4.0)
        np.testing.assert_allclose(final.velocity_body_mps[1:],np.zeros(2),atol=1e-9)
        loads=provider(2.0,final)
        self.assertGreater(loads.force_body_n[0],0.0)


if __name__=='__main__':
    unittest.main()
