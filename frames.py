import numpy as np


def vector3(value):
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ValueError("Expected three finite vector components.")
    return result


def normalize_quaternion(value):
    result = np.asarray(value, dtype=float)
    if result.shape != (4,) or not np.all(np.isfinite(result)):
        raise ValueError("Quaternion must contain four finite components [w,x,y,z].")
    norm = np.linalg.norm(result)
    if norm < 1e-12:
        raise ValueError("A zero quaternion does not define an attitude.")
    return result / norm


def quaternion_product(left, right):
    aw, ax, ay, az = np.asarray(left, dtype=float)
    bw, bx, by, bz = np.asarray(right, dtype=float)
    return np.array(
        [
            aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw,
        ]
    )


LAUNCH_TILT_RAD = -np.pi / 2
LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ = np.array(
    [np.cos(LAUNCH_TILT_RAD / 2), 0.0, np.sin(LAUNCH_TILT_RAD / 2), 0.0]
)


def launch_position_body_m(earth_radius_m, height_above_site_m=0.0):
    radius = float(earth_radius_m)
    height = float(height_above_site_m)
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(height):
        raise ValueError("Expected a positive finite Earth radius and a finite height.")
    return np.array([radius + height, 0.0, 0.0])


def body_to_ecef_matrix(quaternion_wxyz):
    w, x, y, z = normalize_quaternion(quaternion_wxyz)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
        ]
    )


def quaternion_rate_body_to_ecef(quaternion_wxyz, omega_body_radps):
    return 0.5 * quaternion_product(
        normalize_quaternion(quaternion_wxyz), np.r_[0.0, vector3(omega_body_radps)]
    )
