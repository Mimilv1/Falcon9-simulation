from copy import deepcopy
import unittest
import numpy as np

from frames import (
    body_to_ecef_matrix, normalize_quaternion, quaternion_product,
    quaternion_conjugate, quaternion_from_rotation_vector,
    quaternion_to_rotation_vector, earth_rotation_quaternion,
)
from kalman import NavigationKF, _right_jacobian
from sensors import StarTrackerSpec


def error_between(reference, trial):
    parts = [trial.x - reference.x]
    if reference.attitude_enabled:
        parts += [quaternion_to_rotation_vector(quaternion_product(
            quaternion_conjugate(reference.q), trial.q)), trial.bg - reference.bg]
    if reference.calibration_enabled:
        parts += [trial.gyro_sf - reference.gyro_sf, trial.gyro_mis - reference.gyro_mis]
    if reference.accel_sf_enabled:
        parts += [trial.accel_sf - reference.accel_sf]
    if reference.gnss_bias_enabled:
        parts += [trial.gnss_bias - reference.gnss_bias]
    return np.concatenate(parts)


def perturb(filter_, index, value):
    if index < 9:
        filter_.x[index] += value
    elif index < 12:
        rv = np.zeros(3)
        rv[index - 9] = value
        filter_.q = quaternion_product(filter_.q, quaternion_from_rotation_vector(rv))
    elif index < 15:
        filter_.bg[index - 12] += value
    elif index < 18:
        filter_.gyro_sf[index - 15] += value
    elif index < 24:
        filter_.gyro_mis[index - 18] += value
    elif index < 27:
        filter_.accel_sf[index - 24] += value
    else:
        filter_.gnss_bias[index - 27] += value


class NavigationChecks(unittest.TestCase):
    def full_filter(self):
        k = NavigationKF(
            [30., -40., 1000.], [2., -3., 1.], np.ones(9),
            gravity_ecef=[0., 0., -9.81],
            q0=normalize_quaternion([.9, .1, -.2, .3]),
            estimate_gyro_calibration=True, estimate_accel_scale=True,
            gnss_bias_sigma_m=[.6, .6, 1.],
        )
        k.x[6:9] = [.03, -.01, .02]
        k.bg = np.array([.0004, -.0002, .0003])
        k.gyro_sf = np.array([.002, -.0015, .001])
        k.gyro_mis = np.array([.001, -.0008, .0006, .0012, -.001, .0005])
        k.accel_sf = np.array([.001, -.0005, .002])
        k.gnss_bias = np.array([.3, -.4, .2])
        k.prime_imu([1., -.2, 3.], [.03, -.02, .04])
        return k

    def test_exponential_right_jacobian_against_central_difference(self):
        for rv in [np.zeros(3), np.array([.4, -.2, .6])]:
            q = quaternion_from_rotation_vector(rv)
            numeric = np.zeros((3, 3))
            for i in range(3):
                step = np.eye(3)[i] * 1e-6
                values = [quaternion_to_rotation_vector(quaternion_product(
                    quaternion_conjugate(q), quaternion_from_rotation_vector(rv + sign * step)
                )) for sign in [-1, 1]]
                numeric[:, i] = (values[1] - values[0]) / 2e-6
            np.testing.assert_allclose(_right_jacobian(rv), numeric, atol=2e-9)

    def test_full_30_state_transition_against_finite_differences(self):
        initial = self.full_filter()
        nominal = deepcopy(initial)
        nominal.predict([1.1, -.3, 3.2], None, .2, gyro_meas=[.04, -.01, .03])
        numeric = np.zeros((30, 30))
        for i in range(30):
            step = .01 if i < 9 or i >= 27 else .001
            trials = []
            for sign in [-1, 1]:
                trial = deepcopy(initial)
                perturb(trial, i, sign * step)
                trial.predict([1.1, -.3, 3.2], None, .2, gyro_meas=[.04, -.01, .03])
                trials.append(error_between(nominal, trial))
            numeric[:, i] = (trials[1] - trials[0]) / (2 * step)
        np.testing.assert_allclose(nominal.last_transition, numeric, atol=3e-6, rtol=3e-5)

    def test_gnss_zero_latency_jacobian_including_calibration_and_colored_bias(self):
        initial = self.full_filter()
        arm = np.array([3., .3, -.2])
        raw = np.array([.03, -.02, .04])
        rot = initial.rotation
        pos = initial.x[:3] + initial.gnss_bias + rot @ arm
        vel = initial.x[3:6] + rot @ np.cross(initial.calibrated_omega(raw), arm)
        def observe(k):
            return k.update_gnss(pos, vel, None, raw, arm, np.ones(3), np.full(3, .1))
        baseline = deepcopy(initial)
        observe(baseline)
        numeric = np.zeros((6, 30))
        for i in range(30):
            step = 1e-5
            values = []
            for sign in [-1, 1]:
                trial = deepcopy(initial)
                perturb(trial, i, sign * step)
                values.append(observe(trial))
            numeric[:, i] = -(values[1] - values[0]) / (2 * step)
        np.testing.assert_allclose(baseline.last_measurement_jacobian, numeric, atol=2e-8)

    def test_lidar_position_and_attitude_jacobian(self):
        k = self.full_filter()
        radius = 6371000.
        k.x[:3] = [40., -20., radius + 1000.]
        k.q = quaternion_from_rotation_vector([.02, -np.pi / 2, .01])
        arm, axis = np.array([-37.5, .2, -.1]), np.array([-1., 0., 0.])
        rel = k.x[:3] + k.rotation @ arm
        direction = k.rotation @ axis
        b = rel @ direction
        measured = -b - np.sqrt(b*b - (rel @ rel - radius**2))
        def observe(filter_):
            return filter_.update_lidar(measured, None, arm, axis, np.zeros(3), radius, .1)
        baseline = deepcopy(k)
        observe(baseline)
        for i in [0, 1, 2, 9, 10, 11]:
            step = .01 if i < 3 else 1e-4
            values = []
            for sign in [-1, 1]:
                trial = deepcopy(k)
                perturb(trial, i, sign * step)
                values.append(observe(trial)[0])
            numeric = -(values[1] - values[0]) / (2*step)
            self.assertAlmostEqual(baseline.last_measurement_jacobian[0, i], numeric, delta=5e-4)

    def test_delayed_star_attitude_uses_rate_history_and_correct_frame(self):
        k = self.full_filter()
        q_start = k.q.copy()
        earth_rate = np.array([5e-5, 0., 4e-5])
        spec = StarTrackerSpec()
        for i in range(1, 11):
            k.predict([1., -.2, 3.], None, .01,
                      gyro_meas=[.03 + .08*i*.01, -.02, .04])
        measured = quaternion_product(quaternion_product(
            earth_rotation_quaternion(0., earth_rate), q_start), spec.q_body_from_sensor)
        residual = k.update_startracker(measured, 0., spec.q_body_from_sensor,
                       spec.sigma_sensor_rad, earth_rate, .1, gate_nis=16.27)
        self.assertIsNotNone(residual)
        self.assertLess(np.linalg.norm(residual), 1e-12)

    def test_delayed_star_calibration_jacobian(self):
        k = self.full_filter()
        q_start = k.q.copy()
        spec = StarTrackerSpec()
        for i in range(1, 11):
            k.predict([1., -.2, 3.], None, .01, gyro_meas=[.03, -.02 + .01*i, .04])
        measured = quaternion_product(q_start, spec.q_body_from_sensor)
        def observe(filter_):
            return filter_.update_startracker(measured, 0., spec.q_body_from_sensor,
                             spec.sigma_sensor_rad, np.zeros(3), .1, gate_nis=None)
        baseline = deepcopy(k)
        observe(baseline)
        for i in range(9, 24):
            step = 1e-5
            values = []
            for sign in [-1, 1]:
                trial = deepcopy(k)
                perturb(trial, i, sign * step)
                values.append(observe(trial))
            numeric = -(values[1] - values[0]) / (2*step)
            np.testing.assert_allclose(baseline.last_measurement_jacobian[:, i], numeric,
                                       atol=3e-7, rtol=2e-5)

    def test_star_gate_rejects_outlier_without_changing_state(self):
        k = self.full_filter()
        spec = StarTrackerSpec()
        valid = quaternion_product(k.q, spec.q_body_from_sensor)
        self.assertIsNotNone(k.update_startracker(valid, 0., spec.q_body_from_sensor,
                  spec.sigma_sensor_rad, np.zeros(3), 0., gate_nis=16.27))
        before, covariance = k.q.copy(), k.P.copy()
        outlier = quaternion_product(valid, quaternion_from_rotation_vector([.1, 0, 0]))
        self.assertIsNone(k.update_startracker(outlier, 0., spec.q_body_from_sensor,
                  spec.sigma_sensor_rad, np.zeros(3), 0., gate_nis=16.27))
        np.testing.assert_array_equal(k.q, before)
        np.testing.assert_array_equal(k.P, covariance)
        self.assertEqual(k.startracker_accepted, 1)
        self.assertEqual(k.startracker_rejected, 1)

    def test_colored_gnss_bias_has_stationary_covariance_and_exponential_mean(self):
        k = self.full_filter()
        initial = k.gnss_bias.copy()
        for _ in range(50):
            k.predict([0., 0., 0.], None, .2, gyro_meas=[0., 0., 0.])
        np.testing.assert_allclose(k.gnss_bias, initial*np.exp(-10/300), atol=1e-14)
        np.testing.assert_allclose(k.P[k.gnss_bias_slice, k.gnss_bias_slice],
                                   np.diag(k.gnss_bias_sigma**2), atol=2e-14)

    def test_prediction_rejects_invalid_interval_and_singular_accelerometer(self):
        k = self.full_filter()
        with self.assertRaises(ValueError):
            k.predict([0, 0, 0], None, 0, gyro_meas=[0, 0, 0])
        k.accel_sf[0] = -1
        with self.assertRaises(ValueError):
            k.predict([0, 0, 0], None, .01, gyro_meas=[0, 0, 0])


if __name__ == '__main__':
    unittest.main()
