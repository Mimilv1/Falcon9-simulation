import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mass_properties import fixed_cylinder, combine_components
from frames import body_to_ecef_matrix
from load_interface import AppliedForce, assemble_loads
from dynamics import State, Loads, rk4_step


class FixedModel:
    mass_model_name = "fixed_diagnostic"

    def __init__(self):
        self.properties = combine_components(
            [fixed_cylinder("synthetic", 1000, 0, 20, 40, 1)]
        )

    def mass_properties(self, configuration, propellant_kg, consumption_kgps=(0, 0)):
        return self.properties


class CoreMathChecks(unittest.TestCase):
    def test_rk4_convergence_with_fixed_sensor_independent_rotation(self):
        from report.validate_navigation import integrator_case
        result = integrator_case()
        for order in result['observed_orders']:
            self.assertTrue(3.8 <= order <= 4.2)

    def test_parallel_axis_mass_combination(self):
        p = combine_components(
            [
                fixed_cylinder("a", 100, 0, 0, 0, 0.1),
                fixed_cylinder("b", 100, 0, 2, 0, 0.1),
            ]
        )
        self.assertEqual(p.mass_kg, 200)
        np.testing.assert_allclose(p.cg_station_body_m, [1, 0, 0])
        np.testing.assert_allclose(p.inertia_body_kgm2, np.diag([1, 200.5, 200.5]))

    def test_base_side_force_has_negative_z_moment(self):
        p = FixedModel().properties
        loads = assemble_loads(
            p, [AppliedForce([0, 1000, 0], [0, 0, 0])], gravity_body_mps2=[0, 0, 0]
        )
        np.testing.assert_allclose(loads.moment_body_nm, [0, 0, -20000])

    def test_identity_quaternion_maps_body_axes_onto_ecef(self):
        np.testing.assert_allclose(body_to_ecef_matrix([1, 0, 0, 0]), np.eye(3))

    def test_free_fall(self):
        model = FixedModel()
        state = State([100, 0, 0], [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0, 0])
        provider = lambda t, s: Loads([0, 0, 0], [0, 0, 0], [-10, 0, 0], [0, 0])
        result = rk4_step(0, state, 1, model, None, provider)
        np.testing.assert_allclose(result.position_body_m, [95, 0, 0], atol=1e-12)
        np.testing.assert_allclose(
            body_to_ecef_matrix(result.quaternion_body_to_ecef_wxyz)
            @ result.velocity_body_mps,
            [-10, 0, 0],
            atol=1e-12,
        )

    def test_geocentric_trial_positions_do_not_penetrate_ground_when_body_rotates(self):
        model = FixedModel()
        radius = 6371000.0
        state = State([radius + 1, 0, 0], [0, 0, 0], [1, 0, 0, 0], [0, 0, 0.5], [0, 0])

        def provider(t, trial):
            self.assertLess(
                abs(np.linalg.norm(trial.position_body_m) - radius - 1), 1e-6
            )
            return Loads([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0])

        for k in range(20):
            state = rk4_step(k * 0.1, state, 0.1, model, None, provider)
        position = (
            body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
            @ state.position_body_m
        )
        np.testing.assert_allclose(position, [radius + 1, 0, 0], atol=1e-7)

    def test_constant_moment_rotation(self):
        model = FixedModel()
        p = model.properties
        state = State([100, 0, 0], [0, 0, 0], [1, 0, 0, 0], [0, 0, 0], [0, 0])
        provider = lambda t, s: assemble_loads(
            p, [AppliedForce([0, 1000, 0], [0, 0, 0])], gravity_body_mps2=[0, 0, 0]
        )
        alpha = -20000 / p.inertia_body_kgm2[2, 2]
        for k in range(200):
            state = rk4_step(k * 0.01, state, 0.01, model, None, provider)
        nose = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)[:, 0]
        tilt = np.arctan2(np.linalg.norm(nose[1:]), nose[0])
        self.assertAlmostEqual(tilt, abs(alpha) * 2, places=10)
        self.assertAlmostEqual(state.omega_body_radps[2], alpha * 2, places=10)

    def test_rotating_body_storage_preserves_straight_physical_motion(self):
        model = FixedModel()
        state = State([100, 0, 0], [10, 0, 0], [1, 0, 0, 0], [0, 0, 1], [0, 0])
        provider = lambda t, s: Loads([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0])
        for k in range(100):
            state = rk4_step(k * 0.01, state, 0.01, model, None, provider)
        rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        np.testing.assert_allclose(
            rotation @ state.position_body_m, [110, 0, 0], atol=2e-7
        )
        np.testing.assert_allclose(
            rotation @ state.velocity_body_mps, [10, 0, 0], atol=2e-8
        )


if __name__ == "__main__":
    unittest.main()
