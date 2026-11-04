"""Body +x vector/quaternion operations; the fixed reference frame is ECEF.

ECEF convention (project-wide): the Earth centre is the origin (0,0,0) and the
vehicle launches from (0,0,earth_radius) on the +z axis. Earth rotation is not
modeled in v0: the ECEF axes are treated as non-rotating (inertial), so only
the datum changes relative to the former launch-point origin. Quaternions map
current Body components into ECEF. At launch the vehicle stands vertical,
Body +x toward the nose along ECEF +z; that fixed alignment is
LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ below.
"""

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
    return np.array([
        aw*bw - ax*bx - ay*by - az*bz,
        aw*bx + ax*bw + ay*bz - az*by,
        aw*by - ax*bz + ay*bw + az*bx,
        aw*bz + ax*by - ay*bx + az*bw,
    ])


# Fixed launch alignment: -90 degree rotation about the common Body/ECEF y axis
# maps Body +x (nose) onto ECEF +z (radially up at the launch site), Body y
# onto ECEF y and Body z onto ECEF x.
LAUNCH_TILT_RAD = -np.pi / 2
LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ = np.array([
    np.cos(LAUNCH_TILT_RAD / 2), 0.0, np.sin(LAUNCH_TILT_RAD / 2), 0.0
])


def launch_position_body_m(earth_radius_m, height_above_site_m=0.0):
    """ECEF launch point (0,0,R) plus a height h, in Body components at launch attitude.

    With the launch attitude, Body +x is ECEF +z, so the vector from the Earth
    centre to a point h above the launch site is [R+h, 0, 0] in Body components.
    """
    radius = float(earth_radius_m)
    height = float(height_above_site_m)
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(height):
        raise ValueError("Expected a positive finite Earth radius and a finite height.")
    return np.array([radius + height, 0.0, 0.0])


def body_to_ecef_matrix(quaternion_wxyz):
    w, x, y, z = normalize_quaternion(quaternion_wxyz)
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
        [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
        [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)],
    ])


def quaternion_rate_body_to_ecef(quaternion_wxyz, omega_body_radps):
    """Body/ECEF rate in Body components; Earth rotation is zero in v0."""
    return 0.5 * quaternion_product(
        normalize_quaternion(quaternion_wxyz), np.r_[0.0, vector3(omega_body_radps)]
    )
