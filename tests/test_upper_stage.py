from tempfile import TemporaryDirectory
from pathlib import Path
import unittest
import numpy as np

from VehicleData import build_falcon9
from Events import load_mission
from upper_stage import UpperStageProvider, propagate_upper_stage
from frames import body_to_ecef_matrix


class UpperStageChecks(unittest.TestCase):
    def setup_upper(self):
        _, vehicle = build_falcon9().separate()
        vehicle.position[0] += 60000
        state = vehicle.to_state()
        return vehicle, state, load_mission()

    def test_smooth_ignition_and_cutoff_flow(self):
        vehicle, state, plan = self.setup_upper()
        provider = UpperStageProvider(vehicle, state, 144, plan)
        for t in (144, 152, 545, 600):
            self.assertEqual(provider.envelope(t), (0, 0))
        for t in (152.2, 152.5, 544.2, 544.7):
            value, derivative = provider.envelope(t)
            difference = (
                provider.envelope(t + 1e-5)[0] - provider.envelope(t - 1e-5)[0]
            ) / 2e-5
            self.assertAlmostEqual(derivative, difference, places=6)
            self.assertGreater(value, 0)

    def test_upper_flow_only_consumes_stage2_inventory(self):
        vehicle, state, plan = self.setup_upper()
        provider = UpperStageProvider(vehicle, state, 144, plan)
        loads, props, thrust, _ = provider.evaluate(160, state)
        self.assertEqual(loads.consumption_kgps[0], 0)
        self.assertGreater(loads.consumption_kgps[1], 0)
        self.assertEqual(loads.exhaust_ports[0].stage_index, 1)
        self.assertAlmostEqual(props.mass_rate_kgps, -loads.consumption_kgps[1])
        self.assertGreater(thrust, 0)

    def test_propagation_preserves_initial_separation_state_and_burns_fuel(self):
        vehicle, state, plan = self.setup_upper()
        initial_position = (
            body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
            @ state.position_body_m
        )
        with TemporaryDirectory() as temp:
            result = propagate_upper_stage(vehicle, state, 144, 154, plan, Path(temp))
        initial = result["initial"]
        np.testing.assert_allclose(
            [initial["cg_position_ecef_m_" + a] for a in "xyz"], initial_position
        )
        self.assertLess(result["final"]["stage2_propellant_kg"], state.propellant_kg[1])
        self.assertEqual(result["events"][-1]["event"], "stage2_ignition")
        self.assertTrue(result["finite_state"])


if __name__ == "__main__":
    unittest.main()
