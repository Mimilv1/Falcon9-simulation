"""Axisymmetric component mass properties in the team's Body +x convention.

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
    center_axial_m: float
    center_rate_mps: float
    inertia_at_center_kgm2: np.ndarray
    inertia_rate_at_center_kgm2ps: np.ndarray
    length_m: float
    radius_m: float
    mass_acceleration_kgps2: float = 0.0
    center_acceleration_mps2: float = 0.0


@dataclass(frozen=True)
class MassProperties:
    mass_kg: float
    cg_station_body_m: np.ndarray
    inertia_body_kgm2: np.ndarray
    inertia_rate_body_kgm2ps: np.ndarray
    mass_rate_kgps: float = 0.0
    cg_rate_body_mps: np.ndarray | None = None
    components: tuple[MassComponent, ...] = ()
    mass_acceleration_kgps2: float = 0.0
    cg_acceleration_body_mps2: np.ndarray | None = None

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
    return np.diag([radius**2/2] + [(3*radius**2 + length**2)/12]*2)


def fixed_cylinder(name, mass_kg, mass_rate_kgps, center_axial_m, length_m, radius_m,
                   *, mass_acceleration_kgps2=0.0):
    values = [mass_kg, mass_rate_kgps, center_axial_m, length_m, radius_m,mass_acceleration_kgps2]
    if not np.all(np.isfinite(values)) or mass_kg < 0:
        raise ValueError("Component properties must be finite, with nonnegative mass.")
    shape = cylinder_inertia_per_kg(length_m, radius_m)
    return MassComponent(name, float(mass_kg), float(mass_rate_kgps),
                         float(center_axial_m), 0.0, mass_kg*shape,
                         mass_rate_kgps*shape, float(length_m), float(radius_m),float(mass_acceleration_kgps2))


def liquid_column(name, mass_kg, mass_rate_kgps, bottom_axial_m, radius_m,
                  density_kgm3, max_height_m, *, mass_acceleration_kgps2=0.0):
    """A uniform liquid cylinder draining from a fixed bottom along Body +x.

    Rates are derivatives of the geometric mass-property model, not additional
    forces or torques. Their use in open-system equations needs a separate
    momentum-flux derivation. At zero mass, an imposed outflow is a left limit.
    """
    values = [mass_kg, mass_rate_kgps, bottom_axial_m, radius_m, density_kgm3, max_height_m,mass_acceleration_kgps2]
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
    local_rate += np.diag([0.0] + [mass_kg*height*height_rate/6]*2)
    return MassComponent(name, float(mass_kg), float(mass_rate_kgps),
                         bottom_axial_m + height/2, height_rate/2,
                         mass_kg*shape, local_rate, height, radius_m,
                         float(mass_acceleration_kgps2),mass_acceleration_kgps2/(2*mass_per_height))


def combine_components(components):
    """Combine inertias about the instantaneous CG using the parallel-axis theorem."""
    components = tuple(components)
    if not components:
        raise ValueError("At least one mass component is required.")
    mass = sum(c.mass_kg for c in components)
    mass_rate = sum(c.mass_rate_kgps for c in components)
    mass_acceleration = sum(c.mass_acceleration_kgps2 for c in components)
    if mass <= 0:
        raise ValueError("Attached mass must be positive.")
    center = sum(c.mass_kg*c.center_axial_m for c in components)/mass
    center_rate = (sum(c.mass_rate_kgps*c.center_axial_m
                       + c.mass_kg*c.center_rate_mps for c in components)
                   - center*mass_rate)/mass
    center_acceleration = (sum(c.mass_acceleration_kgps2*c.center_axial_m
                               + 2*c.mass_rate_kgps*c.center_rate_mps
                               + c.mass_kg*c.center_acceleration_mps2 for c in components)
                           - center*mass_acceleration - 2*mass_rate*center_rate)/mass
    inertia = np.zeros((3, 3))
    inertia_rate = np.zeros((3, 3))
    for c in components:
        distance = c.center_axial_m - center
        distance_rate = c.center_rate_mps - center_rate
        parallel = np.diag([0.0, distance**2, distance**2])
        parallel_rate = np.diag([0.0] + [2*distance*distance_rate]*2)
        inertia += c.inertia_at_center_kgm2 + c.mass_kg*parallel
        inertia_rate += (c.inertia_rate_at_center_kgm2ps
                         + c.mass_rate_kgps*parallel + c.mass_kg*parallel_rate)
    return MassProperties(mass, np.array([center, 0.0, 0.0]), inertia, inertia_rate,
                          mass_rate, np.array([center_rate, 0.0, 0.0]), components,
                          mass_acceleration,np.array([center_acceleration,0.0,0.0]))
