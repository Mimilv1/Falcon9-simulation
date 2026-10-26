"""Right-handed Body (+x nose) and fixed launch-site ENU transformations.

Quaternions use Hamilton multiplication, scalar first, and rotate Body vectors
into ENU. ENU-to-ECEF is a separate geometric transformation, not a correction
for omitted Earth-rotation or Earth-curvature dynamics.
"""

from dataclasses import dataclass
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


def body_to_enu_matrix(quaternion_wxyz):
    w, x, y, z = normalize_quaternion(quaternion_wxyz)
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
        [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
        [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)],
    ])


def quaternion_rate_body_to_enu(quaternion_wxyz, omega_body_radps):
    """Body/ENU rate in Body components; Earth rotation is zero in v0."""
    return 0.5 * quaternion_product(
        normalize_quaternion(quaternion_wxyz), np.r_[0.0, vector3(omega_body_radps)]
    )


@dataclass(frozen=True)
class LaunchSite:
    latitude_deg: float
    longitude_deg: float
    ellipsoid_height_m: float = 0.0

    def __post_init__(self):
        if not np.all(np.isfinite([
            self.latitude_deg, self.longitude_deg, self.ellipsoid_height_m
        ])) or not -90 <= self.latitude_deg <= 90:
            raise ValueError("Invalid geodetic launch-site coordinates.")

    def origin_ecef_m(self):
        # WGS-84 geodetic origin; h is ellipsoidal height, not MSL altitude.
        a = 6378137.0
        flattening = 1.0 / 298.257223563
        e2 = flattening * (2.0 - flattening)
        latitude, longitude = np.deg2rad([self.latitude_deg, self.longitude_deg])
        n = a / np.sqrt(1.0 - e2 * np.sin(latitude)**2)
        h = self.ellipsoid_height_m
        return np.array([
            (n+h)*np.cos(latitude)*np.cos(longitude),
            (n+h)*np.cos(latitude)*np.sin(longitude),
            (n*(1.0-e2)+h)*np.sin(latitude),
        ])

    def enu_to_ecef_matrix(self):
        latitude, longitude = np.deg2rad([self.latitude_deg, self.longitude_deg])
        s, c = np.sin(latitude), np.cos(latitude)
        sl, cl = np.sin(longitude), np.cos(longitude)
        # Columns are east, north, up, expressed in ECEF axes.
        return np.array([[-sl, -s*cl, c*cl], [cl, -s*sl, c*sl], [0, c, s]])

    def position_enu_to_ecef(self, position_enu_m):
        return self.origin_ecef_m() + self.enu_to_ecef_matrix() @ vector3(position_enu_m)

    def position_ecef_to_enu(self, position_ecef_m):
        return self.enu_to_ecef_matrix().T @ (vector3(position_ecef_m)-self.origin_ecef_m())

    def vector_enu_to_ecef(self, vector_enu):
        """Free vectors have no origin offset; applies to Earth-relative velocity."""
        return self.enu_to_ecef_matrix() @ vector3(vector_enu)

    def attitude_body_to_ecef_matrix(self, quaternion_body_to_enu_wxyz):
        return self.enu_to_ecef_matrix() @ body_to_enu_matrix(quaternion_body_to_enu_wxyz)
