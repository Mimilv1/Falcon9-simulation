"""Axisymmetric component mass properties in the team's Body +z convention.

Liquid columns are a settled-fluid geometric approximation. This module does
not model slosh, exhaust angular-momentum flux, or rigid-body propagation.
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MassComponent:
    name: str
    mass_kg: float
    mass_rate_kgps: float
    center_z_m: float
    center_rate_mps: float
    inertia_at_center_kgm2: np.ndarray
    inertia_rate_at_center_kgm2ps: np.ndarray
    length_m: float
    radius_m: float


@dataclass(frozen=True)
class MassProperties:
    mass_kg: float
    cg_station_body_m: np.ndarray
    inertia_body_kgm2: np.ndarray
    inertia_rate_body_kgm2ps: np.ndarray
    mass_rate_kgps: float = 0.0
    cg_rate_body_mps: np.ndarray | None = None
    components: tuple[MassComponent, ...] = ()

    def arm_from_station(self, station_body_m):
        """Convert a structural-datum station to a force arm about current CG."""
        station = np.asarray(station_body_m, dtype=float)
        if station.shape != (3,) or not np.all(np.isfinite(station)):
            raise ValueError("A station must have three finite Body components.")
        return station - self.cg_station_body_m


def cylinder_inertia_per_kg(length_m, radius_m):
    length, radius = float(length_m), float(radius_m)
    if not np.all(np.isfinite([length, radius])) or length < 0 or radius <= 0:
        raise ValueError("Cylinder length must be nonnegative and radius positive.")
    return np.diag([(3*radius**2 + length**2)/12]*2 + [radius**2/2])


def fixed_cylinder(name, mass_kg, mass_rate_kgps, center_z_m, length_m, radius_m):
    values = [mass_kg, mass_rate_kgps, center_z_m, length_m, radius_m]
    if not np.all(np.isfinite(values)) or mass_kg < 0:
        raise ValueError("Component properties must be finite, with nonnegative mass.")
    shape = cylinder_inertia_per_kg(length_m, radius_m)
    return MassComponent(name, float(mass_kg), float(mass_rate_kgps),
                         float(center_z_m), 0.0, mass_kg*shape,
                         mass_rate_kgps*shape, float(length_m), float(radius_m))


def liquid_column(name, mass_kg, mass_rate_kgps, bottom_z_m, radius_m,
                  density_kgm3, max_height_m):
    """A uniform liquid cylinder draining from a fixed bottom along Body +z.

    Rates are derivatives of the geometric mass-property model, not additional
    forces or torques. Their use in open-system equations needs a separate
    momentum-flux derivation. At zero mass, an imposed outflow is a left limit.
    """
    values = [mass_kg, mass_rate_kgps, bottom_z_m, radius_m, density_kgm3, max_height_m]
    if (not np.all(np.isfinite(values)) or mass_kg < 0
            or min(radius_m, density_kgm3, max_height_m) <= 0):
        raise ValueError("Invalid liquid-column geometry or mass.")
    mass_per_height = density_kgm3*np.pi*radius_m**2
    height = mass_kg/mass_per_height
    if height > max_height_m + 1e-8:
        raise ValueError(f"{name} exceeds the assumed tank capacity.")
    height_rate = mass_rate_kgps/mass_per_height
    shape = cylinder_inertia_per_kg(height, radius_m)
    local_rate = mass_rate_kgps*shape
    local_rate += np.diag([mass_kg*height*height_rate/6]*2 + [0.0])
    return MassComponent(name, float(mass_kg), float(mass_rate_kgps),
                         bottom_z_m + height/2, height_rate/2,
                         mass_kg*shape, local_rate, height, radius_m)


def combine_components(components):
    """Combine inertias about the instantaneous CG using the parallel-axis theorem."""
    components = tuple(components)
    if not components:
        raise ValueError("At least one mass component is required.")
    mass = sum(c.mass_kg for c in components)
    mass_rate = sum(c.mass_rate_kgps for c in components)
    if mass <= 0:
        raise ValueError("Attached mass must be positive.")
    center = sum(c.mass_kg*c.center_z_m for c in components)/mass
    center_rate = (sum(c.mass_rate_kgps*c.center_z_m
                       + c.mass_kg*c.center_rate_mps for c in components)
                   - center*mass_rate)/mass
    inertia = np.zeros((3, 3))
    inertia_rate = np.zeros((3, 3))
    for c in components:
        distance = c.center_z_m - center
        distance_rate = c.center_rate_mps - center_rate
        parallel = np.diag([distance**2, distance**2, 0.0])
        parallel_rate = np.diag([2*distance*distance_rate]*2 + [0.0])
        inertia += c.inertia_at_center_kgm2 + c.mass_kg*parallel
        inertia_rate += (c.inertia_rate_at_center_kgm2ps
                         + c.mass_rate_kgps*parallel + c.mass_kg*parallel_rate)
    return MassProperties(mass, np.array([0.0, 0.0, center]), inertia, inertia_rate,
                          mass_rate, np.array([0.0, 0.0, center_rate]), components)
