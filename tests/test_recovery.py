import unittest
import numpy as np
from dynamics import State
from Events import load_mission, MissionSequence
from VehicleData import build_falcon9
from recovery import RecoveryController, recovery_plan, smooth_command
from frames import body_to_ecef_matrix, quaternion_product


class RecoveryChecks(unittest.TestCase):
    def setup_booster(self, phase="entry_burn"):
        plan = recovery_plan(load_mission())
        v, _ = build_falcon9().separate()
        v.stages[0].propellant_mass = 30000
        v.position = np.array([6401000.0, 0, 0])
        v.velocity = np.array([-500.0, 20, 0])
        c = RecoveryController(plan)
        c.time = 300
        c.provider(v, phase, plan)
        flow, derivative = c.flow_at(300)
        s = v.to_state(
            consumption_by_stage={"stage1": flow[0]},
            consumption_derivative_by_stage={"stage1": derivative[0]},
        )
        return plan, v, c, s

    def test_throttle_ramp_has_matching_flow_derivative(self):
        for t in [0.1, 0.3, 0.7, 0.9]:
            y, d = smooth_command(t, 0, 1, 0, 0, 0.8)
            eps = 1e-5
            numeric = (
                smooth_command(t + eps, 0, 1, 0, 0, 0.8)[0]
                - smooth_command(t - eps, 0, 1, 0, 0, 0.8)[0]
            ) / (2 * eps)
            self.assertAlmostEqual(d, numeric, places=8)
        self.assertEqual(smooth_command(1, 0, 1, 0, 0, 0.8), (0.8, 0))

    def test_entry_ports_are_balanced_and_flow_matches_mass_properties(self):
        _, v, c, s = self.setup_booster()
        for t in [300, 300.3, 301]:
            e = c.evaluate(t, s)
            self.assertEqual(len(e.loads.exhaust_ports), 3)
            self.assertAlmostEqual(
                e.properties.mass_rate_kgps, -e.loads.consumption_kgps.sum()
            )
            np.testing.assert_allclose(
                sum(
                    p.mass_flow_kgps * p.station_body_m[1:]
                    for p in e.loads.exhaust_ports
                ),
                np.zeros(2),
                atol=1e-10,
            )
            self.assertTrue(np.all(np.isfinite(e.loads.moment_body_nm)))

    def test_shutdown_cannot_make_negative_flow_at_roundoff_boundary(self):
        _, v, c, s = self.setup_booster()
        c.ramp = (300, 1, 0.8, 0, 0)
        for t in np.linspace(300.9999999, 301.0000001, 100):
            self.assertGreaterEqual(c.flow_at(t)[0][0], 0)
            c.evaluate(t, s)

    def test_geometric_contact_is_not_cg_zero(self):
        _, v, c, s = self.setup_booster()
        v.position = np.array([6371001.0, 0, 0])
        v.velocity = np.zeros(3)
        s = v.to_state()
        props = v.get_mass_properties()
        self.assertLess(abs(c.contact_altitude(s, props)), 2e-6)
        self.assertGreater(c.environment.altitude(s, props), 10)

    def test_return_ignitions_use_state_conditions_with_source_markers_retained(self):
        plan = recovery_plan(load_mission())
        events = {e["name"]: e for e in plan["transitions"]}
        self.assertIsNone(events["entry_start"]["time_s"])
        self.assertEqual(events["entry_start"]["reference_time_s"], 392)
        self.assertIsNone(events["landing_start"]["time_s"])
        sequence = MissionSequence(plan)
        for e in sequence.transitions[:5]:
            sequence.accept(e.name, e.time_s or 200, dict.fromkeys(e.guards, True))
        self.assertIsNone(sequence.due(392))
        self.assertEqual(sequence.due(395, {"entry_start": True}).name, "entry_start")

    def test_force_inversion_respects_limits_instead_of_pointing_thrust_at_net_force(
        self,
    ):
        plan, v, c, state = self.setup_booster("entry_burn")
        v.stages[0].propellant_mass = 15000
        v.position = np.array([6372500.0, 0, 0])
        v.velocity = np.array([-100.0, 10, 0])
        c.time = 0
        c.state = v.to_state()
        c.properties = v.get_mass_properties()
        c.provider(v, "landing_burn", plan)
        state = v.to_state()
        props = v.get_mass_properties()
        p, velocity, up, *_ = c.geometry(state, props)
        desired = 5 * up + np.array([1.0, 0, 0])
        nose, throttle = c.solve_force(state, props, desired)
        self.assertTrue(v.stages[0].engines[0].spec.throttle_min <= throttle <= 1)
        self.assertLessEqual(np.degrees(np.arccos(np.clip(nose @ up, -1, 1))), 20.00001)
        self.assertTrue(np.all(np.isfinite(nose)))

    def test_lowest_leg_clearance_rate_matches_its_motion(self):
        _, vehicle, controller, _ = self.setup_booster()
        vehicle.omega = np.array([0, 0.03, 0.01])
        vehicle.quaternion = quaternion_product(
            vehicle.quaternion, np.array([np.cos(0.1), 0, np.sin(0.1), 0])
        )
        props = vehicle.get_mass_properties()
        rot = body_to_ecef_matrix(vehicle.quaternion)
        state = State(
            rot.T @ np.array([0.0, 0.0, 6401000.0]),
            rot.T @ np.array([20.0, 0.0, -500.0]),
            vehicle.quaternion,
            vehicle.omega,
            [vehicle.stages[0].propellant_mass, 0],
        )
        position, velocity, _ = controller.contact_state(state, props)
        rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        cg = rotation @ state.position_body_m
        cg_velocity = rotation @ state.velocity_body_mps
        dt = 0.001
        heights = []
        rate = np.linalg.norm(state.omega_body_radps)
        for sign in (-1, 1):
            delta = np.r_[
                np.cos(rate * dt / 2),
                sign * np.sin(rate * dt / 2) * state.omega_body_radps / rate,
            ]
            q = quaternion_product(state.quaternion_body_to_ecef_wxyz, delta)
            rot = body_to_ecef_matrix(q)
            trial = State(
                rot.T @ (cg + sign * dt * cg_velocity),
                rot.T @ cg_velocity,
                q,
                state.omega_body_radps,
                state.propellant_kg,
            )
            heights.append(controller.contact_altitude(trial, props))
        finite_difference = (heights[1] - heights[0]) / (2 * dt)
        self.assertAlmostEqual(
            finite_difference,
            float(position @ velocity / np.linalg.norm(position)),
            delta=2e-5,
        )

    def test_original_acceptance_limits_are_retained(self):
        self.assertEqual(
            load_mission()["landing_criteria"]["limits"],
            {
                "vertical_speed_mps": 2.0,
                "horizontal_speed_mps": 1.0,
                "tilt_deg": 5.0,
                "angular_rate_degps": 1.0,
                "target_error_m": 10.0,
            },
        )

    def test_allocated_vertical_force_matches_the_propulsion_and_aero_models(self):
        plan, vehicle, controller, _ = self.setup_booster()
        vehicle.stages[0].propellant_mass = 15000
        vehicle.position = np.array([6372500.0, 0, 0])
        vehicle.velocity = np.array([-100.0, 10, 0])
        original = vehicle.to_state()
        props = vehicle.get_mass_properties()
        controller.time = 0
        controller.state = original
        controller.properties = props
        controller.provider(vehicle, "landing_burn", plan)
        position, velocity, up, *_ = controller.geometry(original, props)
        desired = 5 * up + np.array([1.0, 0, 0])
        nose, throttle = controller.solve_force(original, props, desired)
        rot = body_to_ecef_matrix(original.quaternion_body_to_ecef_wxyz)
        first = rot[:, 0]
        axis = np.cross(first, nose)
        angle = np.arccos(np.clip(first @ nose, -1, 1))
        delta = np.r_[
            np.cos(angle / 2),
            axis / max(np.linalg.norm(axis), 1e-12) * np.sin(angle / 2),
        ]
        q = quaternion_product(delta, original.quaternion_body_to_ecef_wxyz)
        turned = body_to_ecef_matrix(q)
        state = State(
            turned.T @ (rot @ original.position_body_m),
            turned.T @ (rot @ original.velocity_body_mps),
            q,
            [0, 0, 0],
            original.propellant_kg,
        )
        air = controller.environment.atmosphere(
            controller.environment.altitude(state, props)
        )
        aero = controller.aero.evaluate(
            state,
            props,
            air,
            controller.environment.wind_body(state),
            "booster_landing",
        )
        from dataclasses import replace
        from propulsion import evaluate_steady_engine

        pressure = controller.environment.atmosphere(
            controller.environment.altitude(state, props, [0, 0, 0])
        ).pressure_pa
        output = evaluate_steady_engine(
            replace(vehicle.stages[0].engines[0], running=True, throttle=throttle),
            pressure,
            stage_index=0,
            stage_offset_m=0,
            nozzle_exit_axial_m=0,
        )
        actual = (
            turned
            @ (aero.force_body_n + output.applied_force.force_body_n)
            / props.mass_kg
            + controller.environment.gravity_ecef(rot @ original.position_body_m)
            + controller.environment.frame_acceleration_ecef(
                rot @ original.position_body_m, rot @ original.velocity_body_mps
            )
        )
        self.assertAlmostEqual(float(actual @ up), float(desired @ up), delta=0.02)

    def test_landing_uses_fin_force_only_after_virtual_tvc_authority_is_used(self):
        plan, vehicle, controller, state = self.setup_booster()
        controller.state = state
        controller.properties = vehicle.get_mass_properties()
        controller.provider(vehicle, "landing_burn", plan)
        evaluation = controller.evaluate(controller.time + 1, state)
        self.assertAlmostEqual(evaluation.aero.grid_fin_control_force_n, 0.0, places=6)
        self.assertEqual(plan["recovery"]["max_control_moment_nm"], 2000000)
        self.assertEqual(plan["recovery"]["max_angular_accel_radps2"], 0.12)
        self.assertEqual(plan["recovery"]["max_angular_rate_radps"], 0.22)

    def test_return_budget_fits_revised_meco_inventory(self):
        p = load_mission()["recovery"]
        total = sum(
            p[k]
            for k in [
                "boostback_budget_kg",
                "entry_budget_kg",
                "landing_budget_kg",
                "reserve_kg",
            ]
        )
        self.assertEqual(total, 66000)
        self.assertLess(total, 68889.374 - 1202.128)


if __name__ == "__main__":
    unittest.main()
