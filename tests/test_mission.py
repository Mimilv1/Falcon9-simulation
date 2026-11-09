from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import json
import unittest
from unittest.mock import patch
import numpy as np

from Events import (
    EventConditionError,
    MissionSequence,
    landing_verdict,
    load_mission,
    restart_conditions,
)
from VehicleData import build_falcon9
from dynamics import Loads
from mission import (
    default_provider,
    run_mission,
    terminal_plan,
    prepare_launch,
    structural_state,
    NavigationRun,
)


class MissionChecks(unittest.TestCase):
    def test_reference_clock_and_unknown_endpoints_are_not_mixed(self):
        plan = load_mission()
        events = {e["name"]: e["time_s"] for e in plan["transitions"]}
        self.assertEqual(
            [
                events[k]
                for k in (
                    "meco",
                    "stage_separation",
                    "boostback_start",
                    "entry_start",
                    "landing_start",
                )
            ],
            [141, 144, 161, 392, 453],
        )
        self.assertIsNone(events["boostback_end"])
        self.assertIsNone(events["entry_end"])
        self.assertIsNone(events["touchdown"])
        self.assertEqual(
            events["shutdown_start"], events["meco"] - plan["shutdown"]["duration_s"]
        )

    def test_sequence_rejects_skipped_early_and_unsafe_events(self):
        sequence = MissionSequence()
        with self.assertRaises(EventConditionError):
            sequence.accept("stage_separation", 144, {})
        with self.assertRaises(EventConditionError):
            sequence.accept("shutdown_start", 139, {"fuel_available": True})
        sequence.accept("shutdown_start", 140, {"fuel_available": True})
        with self.assertRaises(EventConditionError):
            sequence.accept("meco", 141, {"zero_flow": False})
        self.assertEqual(sequence.phase, "engine_shutdown")

    def test_restart_requires_declared_budget_and_ready_hardware(self):
        vehicle, _ = build_falcon9().separate()
        state = vehicle.to_state()
        loads = Loads([0, 0, 0], [0, 0, 0], [-9.80665, 0, 0], [0, 0])
        conditions = restart_conditions(vehicle, state, loads)
        self.assertFalse(conditions["fuel_budget_available"])
        self.assertFalse(conditions["attitude_ready"])
        self.assertFalse(conditions["propulsion_ready"])
        conditions = restart_conditions(
            vehicle,
            state,
            loads,
            required_propellant_kg=1000,
            attitude_ready=True,
            propulsion_ready=True,
        )
        self.assertTrue(conditions["fuel_budget_available"])
        with self.assertRaises(EventConditionError):
            MissionSequence().accept("boostback_start", 161, conditions)

    def test_ground_contact_is_not_automatically_a_successful_landing(self):
        limits = load_mission()["landing_criteria"]["limits"]
        metrics = {k: v / 2 for k, v in limits.items()}
        metrics.update(ground_contact=True, legs_deployed=True)
        self.assertTrue(landing_verdict(metrics, limits)["landing_success"])
        metrics["vertical_speed_mps"] = 20
        self.assertFalse(landing_verdict(metrics, limits)["landing_success"])
        self.assertFalse(
            landing_verdict({"ground_contact": True}, limits)["landing_success"]
        )

    def test_missing_burn_end_cannot_be_replaced_by_next_ignition_clock(self):
        sequence = MissionSequence()
        for event in sequence.transitions[:4]:
            sequence.accept(event.name, event.time_s, dict.fromkeys(event.guards, True))
        self.assertIsNone(sequence.due(392))
        self.assertEqual(sequence.next_clock_time, 392)
        with self.assertRaises(EventConditionError):
            sequence.accept("entry_start", 392, {})
        event = sequence.due(220, {"burn_complete": True})
        sequence.accept(event.name, 220, dict.fromkeys(event.guards, True))
        self.assertEqual(sequence.due(392).name, "entry_start")

    def shortened_plan(self):
        plan = deepcopy(load_mission())
        plan["shutdown"]["duration_s"] = 0.01
        plan["initial_state"]["structural_base_position_body_m"][0] += 10.0
        for event, time_s in zip(plan["transitions"][:3], (0.04, 0.05, 0.06)):
            event["time_s"] = time_s
        return plan

    def test_real_vehicle_sensor_pipeline_and_endpoint_separation(self):
        with TemporaryDirectory() as temp:
            result = run_mission(temp, until_s=0.06, plan=self.shortened_plan())
            self.assertEqual(result["status"], "requested_interval_complete")
            self.assertEqual(result["phase"], "booster_coast")
            self.assertFalse(result["full_mission_completed"])
            self.assertEqual(result["sensor_counts"]["imu"], 7)
            self.assertEqual(result["final_truth"]["stage2_propellant_kg"], 0)
            events = json.loads((Path(temp) / "events.json").read_text())
            split = events[-1]
            self.assertEqual(split["event"], "stage_separation")
            self.assertLess(
                split["conservation_audit"]["linear_relative_residual"], 1e-10
            )
            self.assertLess(
                split["conservation_audit"]["structural_position_continuity_error_m"],
                1e-8,
            )
            with (Path(temp) / "navigation.csv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 2)
            self.assertIsNone(result["plan"]["transitions"][4]["time_s"])

    def test_domain_limit_saves_last_accepted_state_without_claiming_completion(self):
        class LimitedProvider:
            def __init__(self, vehicle, phase, plan):
                self.base = default_provider(vehicle, phase, plan)

            def flow_at(self, time_s):
                return self.base.flow_at(time_s)

            def evaluate(self, time_s, state):
                if time_s > 0.012:
                    raise ValueError(
                        "The first environment model supports altitude 0..11000 m."
                    )
                return self.base.evaluate(time_s, state)

            def __call__(self, time_s, state):
                return self.evaluate(time_s, state).loads

        with TemporaryDirectory() as temp:
            result = run_mission(temp, until_s=0.02, provider_factory=LimitedProvider)
            self.assertEqual(result["status"], "model_domain_limit")
            self.assertAlmostEqual(result["actual_end_s"], 0.01)
            self.assertFalse(result["full_mission_completed"])
            self.assertTrue((Path(temp) / "summary.json").exists())

    def test_launch_base_is_at_rest_with_moving_cg_and_rotation_enabled(self):
        plan = load_mission()
        vehicle = prepare_launch(plan)
        provider = default_provider(vehicle, "powered_ascent", plan)
        flow, derivative = provider.flow_at(0)
        state = vehicle.to_state(
            consumption_by_stage={"stage1": flow[0]},
            consumption_derivative_by_stage={"stage1": derivative[0]},
        )
        ev = provider.evaluate(0, state)
        p, v = structural_state(state, ev.properties, [0, 0, 0])
        np.testing.assert_allclose(p, [0, 0, 6371000], atol=1e-8)
        np.testing.assert_allclose(v, np.zeros(3), atol=1e-12)
        self.assertGreater(np.linalg.norm(ev.loads.frame_acceleration_body_mps2), 0.01)
        self.assertAlmostEqual(
            ev.properties.cg_rate_body_mps[0], state.velocity_body_mps[0]
        )

    def test_navigation_is_not_reinitialized_at_separation(self):
        with TemporaryDirectory() as temp, patch(
            "mission.NavigationRun", wraps=NavigationRun
        ) as factory:
            result = run_mission(temp, until_s=0.08, plan=self.shortened_plan())
            self.assertEqual(factory.call_count, 1)
            self.assertEqual(result["sensor_counts"]["imu"], 9)
            self.assertEqual(result["phase"], "booster_coast")

    def test_terminal_case_is_independent_and_only_fires_central_engine(self):
        plan = terminal_plan()
        vehicle = prepare_launch(plan)
        self.assertEqual([s.name for s in vehicle.stages], ["stage1"])
        self.assertEqual(sum(e.running for e in vehicle.stages[0].engines), 1)
        self.assertEqual(vehicle.payload_mass, 0)
        with TemporaryDirectory() as temp:
            result = run_mission(temp, until_s=0.02, plan=plan)
            self.assertFalse(result["full_mission_completed"])
            self.assertEqual(result["phase"], "terminal_braking")
            self.assertAlmostEqual(
                result["initial_truth"]["base_altitude_m"], 1500, places=6
            )
            self.assertAlmostEqual(
                result["initial_truth"]["base_vertical_velocity_mps"], -100, places=8
            )
            self.assertGreater(
                result["final_truth"]["base_vertical_velocity_mps"], -100
            )
            self.assertGreater(result["lidar_updates"], 0)
            self.assertEqual(result["reference_comparison"]["count"], 0)

    def test_imu_undersampling_is_rejected(self):
        with TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "IMU period"):
                run_mission(temp, until_s=1, step_s=0.02)

    def test_runner_does_not_accept_a_fuel_ready_flag_without_a_budget(self):
        plan = self.shortened_plan()
        plan["transitions"][3]["time_s"] = 0.065

        def readiness(*args):
            return {
                "fuel_budget_available": True,
                "attitude_ready": True,
                "propulsion_ready": True,
            }

        with TemporaryDirectory() as temp:
            result = run_mission(
                temp, until_s=0.065, plan=plan, conditions_provider=readiness
            )
            self.assertEqual(result["status"], "event_conditions_unmet")
            self.assertIn("fuel_budget_available", result["detail"])
            self.assertEqual(result["phase"], "booster_coast")


if __name__ == "__main__":
    unittest.main()
