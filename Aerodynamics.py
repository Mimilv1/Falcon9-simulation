from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np


class MachDragCurve:
    def __init__(self, mach, coefficients):
        self.mach = np.asarray(mach, dtype=float)
        self.coefficients = np.asarray(coefficients, dtype=float)
        if (
            self.mach.ndim != 1
            or len(self.mach) < 2
            or self.mach.shape != self.coefficients.shape
            or not np.all(np.isfinite(self.mach))
            or not np.all(np.isfinite(self.coefficients))
            or np.any(np.diff(self.mach) <= 0)
            or self.mach[0] < 0
            or np.any(self.coefficients < 0)
        ):
            raise ValueError("Invalid Mach/drag-coefficient table.")

    @classmethod
    def from_json(cls, path, reference_area_m2):
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        if not np.isclose(
            record["reference_area_m2"], reference_area_m2, rtol=0, atol=1e-9
        ):
            raise ValueError("Cd curve and reference-area conventions must match.")
        return cls(record["mach"], record["drag_coefficient"])

    def __call__(self, mach):
        if not np.isfinite(mach) or mach < self.mach[0] or mach > self.mach[-1]:
            raise ValueError("Mach number is outside the adopted drag curve.")
        return float(np.interp(mach, self.mach, self.coefficients))


@dataclass(frozen=True)
class AeroOutput:
    force_body_n: np.ndarray
    relative_speed_mps: float
    mach: float
    dynamic_pressure_pa: float
    angle_of_attack_rad: float
    drag_coefficient: float
    normal_force_coefficient: float = 0.0
    free_moment_body_nm: np.ndarray | None = None
    grid_fin_control_force_n: float = 0.0


class DragOnlyAerodynamics:
    def __init__(self, reference_area_m2, drag_coefficient):
        if not np.isfinite(reference_area_m2) or reference_area_m2 <= 0:
            raise ValueError("Invalid reference area or drag coefficient.")
        if not callable(drag_coefficient) and (
            not np.isfinite(drag_coefficient) or drag_coefficient < 0
        ):
            raise ValueError("Invalid reference area or drag coefficient.")
        self.reference_area_m2 = float(reference_area_m2)
        self.drag_coefficient = drag_coefficient

    def evaluate(self, state, properties, atmosphere, wind_body_mps):

        migration = properties.cg_rate_body_mps
        relative = state.velocity_body_mps - migration - wind_body_mps
        speed = float(np.linalg.norm(relative))
        mach = speed / atmosphere.sound_speed_mps
        cd = float(
            self.drag_coefficient(mach)
            if callable(self.drag_coefficient)
            else self.drag_coefficient
        )
        if not np.isfinite(cd) or cd < 0:
            raise ValueError("Drag coefficient must be finite and nonnegative.")
        dynamic_pressure = 0.5 * atmosphere.density_kgm3 * speed * speed
        drag = (
            -dynamic_pressure * self.reference_area_m2 * cd * relative / speed
            if speed > 1e-10
            else np.zeros(3)
        )
        aoa = (
            float(np.arctan2(np.linalg.norm(relative[1:]), relative[0]))
            if speed > 1e-10
            else 0.0
        )
        return AeroOutput(drag, speed, mach, dynamic_pressure, aoa, cd)


@dataclass(frozen=True)
class SideSegment:
    bottom_m: float
    top_m: float
    width_m: float
    crossflow_drag_coefficient: float

    @property
    def side_area_m2(self):
        return float(self.width_m * (self.top_m - self.bottom_m))

    @property
    def centroid_m(self):
        return float(0.5 * (self.bottom_m + self.top_m))


@dataclass(frozen=True)
class GridFinControl:
    fin_count: int
    azimuths_rad: tuple
    fin_area_m2: float
    lift_slope_per_rad: float
    pivot_radius_m: float
    station_m: float
    maximum_deflection_rad: float

    def __post_init__(self):
        azimuths = tuple(self.azimuths_rad)
        if (
            not isinstance(self.fin_count, int)
            or self.fin_count < 1
            or len(azimuths) != self.fin_count
            or not np.all(np.isfinite(azimuths))
            or not np.all(
                np.isfinite(
                    [
                        self.fin_area_m2,
                        self.lift_slope_per_rad,
                        self.pivot_radius_m,
                        self.station_m,
                        self.maximum_deflection_rad,
                    ]
                )
            )
            or not self.fin_area_m2 > 0
            or not self.lift_slope_per_rad > 0
            or not self.pivot_radius_m >= 0
            or not self.maximum_deflection_rad > 0
        ):
            raise ValueError("Invalid grid fin control declaration.")
        object.__setattr__(self, "azimuths_rad", azimuths)


@dataclass(frozen=True)
class AeroConfiguration:
    name: str
    reference_area_m2: float
    reference_length_m: float
    axial_drag: MachDragCurve
    potential_normal_slope_per_rad: float
    potential_center_station_m: float
    side_segments: tuple[SideSegment, ...]
    pitch_damping_coefficient: float
    roll_damping_coefficient: float = 0.0
    grid_fin_control: "GridFinControl | None" = None

    def __post_init__(self):
        segments = tuple(self.side_segments)
        if (
            not isinstance(self.name, str)
            or not self.name
            or not np.all(
                np.isfinite(
                    [
                        self.reference_area_m2,
                        self.reference_length_m,
                        self.potential_normal_slope_per_rad,
                        self.potential_center_station_m,
                        self.pitch_damping_coefficient,
                        self.roll_damping_coefficient,
                    ]
                )
            )
            or self.reference_area_m2 <= 0
            or self.reference_length_m <= 0
            or not 0 <= self.potential_normal_slope_per_rad
            or not 0 <= self.potential_center_station_m <= self.reference_length_m
            or self.pitch_damping_coefficient > 0
            or self.roll_damping_coefficient > 0
        ):
            raise ValueError("Invalid aerodynamic configuration parameters.")
        if not isinstance(self.axial_drag, MachDragCurve):
            raise ValueError("The axial drag law must be a MachDragCurve.")
        if self.axial_drag.mach[0] > 0:
            raise ValueError(
                "The axial drag curve must cover Mach 0 for the rest case."
            )
        if not segments:
            raise ValueError("A configuration needs at least one side segment.")
        for segment in segments:
            if (
                not isinstance(segment, SideSegment)
                or not np.all(
                    np.isfinite(
                        [
                            segment.bottom_m,
                            segment.top_m,
                            segment.width_m,
                            segment.crossflow_drag_coefficient,
                        ]
                    )
                )
                or not 0 <= segment.bottom_m < segment.top_m <= self.reference_length_m
                or segment.width_m <= 0
                or not 0 <= segment.crossflow_drag_coefficient
            ):
                raise ValueError("Invalid configuration side segment.")
        fin = self.grid_fin_control
        if fin is not None and (
            not isinstance(fin, GridFinControl)
            or not 0 <= fin.station_m <= self.reference_length_m
        ):
            raise ValueError(
                "The grid fin control station must be a GridFinControl "
                "station within the configuration length."
            )
        object.__setattr__(self, "side_segments", segments)


class Aerodynamics:

    def __init__(self, configurations):
        adopted = dict(configurations)
        if not adopted:
            raise ValueError("At least one aerodynamic configuration is required.")
        for name, configuration in adopted.items():
            if (
                not isinstance(configuration, AeroConfiguration)
                or configuration.name != name
            ):
                raise ValueError(
                    "Configurations must be AeroConfiguration objects keyed by name."
                )
        self.configurations = adopted

    @classmethod
    def from_json(cls, path):
        source = Path(path)
        record = json.loads(source.read_text(encoding="utf-8"))
        base_directory = source.resolve().parent
        configurations = {}
        for name, spec in record["configurations"].items():
            curve_record = spec["axial_drag_curve"]
            if "file" in curve_record:
                axial = MachDragCurve.from_json(
                    base_directory / curve_record["file"], spec["reference_area_m2"]
                )
            else:
                axial = MachDragCurve(
                    curve_record["mach"], curve_record["drag_coefficient"]
                )
            segments = tuple(
                SideSegment(
                    segment["bottom_m"],
                    segment["top_m"],
                    segment["width_m"],
                    segment["crossflow_drag_coefficient"],
                )
                for segment in spec["side_segments"]
            )
            fin_spec = spec.get("grid_fin_control")
            grid_fin = (
                None
                if fin_spec is None
                else GridFinControl(
                    fin_count=fin_spec["fin_count"],
                    azimuths_rad=tuple(fin_spec["azimuths_rad"]),
                    fin_area_m2=fin_spec["fin_area_m2"],
                    lift_slope_per_rad=fin_spec["lift_slope_per_rad"],
                    pivot_radius_m=fin_spec["pivot_radius_m"],
                    station_m=fin_spec["station_m"],
                    maximum_deflection_rad=fin_spec["maximum_deflection_rad"],
                )
            )
            configurations[name] = AeroConfiguration(
                name=name,
                reference_area_m2=spec["reference_area_m2"],
                reference_length_m=spec["reference_length_m"],
                axial_drag=axial,
                potential_normal_slope_per_rad=spec["potential_normal_slope_per_rad"],
                potential_center_station_m=spec["potential_center_station_m"],
                side_segments=segments,
                pitch_damping_coefficient=spec["pitch_damping_coefficient"],
                roll_damping_coefficient=spec.get("roll_damping_coefficient", 0.0),
                grid_fin_control=grid_fin,
            )
        return cls(configurations)

    @staticmethod
    def _validated_fin_commands(configuration, grid_fin_deflections_rad):
        control = configuration.grid_fin_control
        if grid_fin_deflections_rad is None:
            return np.zeros(0 if control is None else control.fin_count)
        values = np.asarray(grid_fin_deflections_rad, dtype=float)
        if control is None:
            if values.size and np.any(values):
                raise ValueError("This configuration declares no grid fin control.")
            return np.zeros(0)
        if values.shape != (control.fin_count,):
            raise ValueError("Grid fin commands need one deflection per declared fin.")
        if not np.all(np.isfinite(values)):
            raise ValueError("Grid fin commands must be finite.")
        if np.any(np.abs(values) > control.maximum_deflection_rad):
            raise ValueError("Grid fin commands exceed the declared deflection travel.")
        return values

    def evaluate(
        self,
        state,
        properties,
        atmosphere,
        wind_body_mps,
        configuration_name,
        grid_fin_deflections_rad=None,
    ):
        configuration = self.configurations.get(configuration_name)
        if configuration is None:
            raise ValueError(
                "Unknown aerodynamic configuration: " + str(configuration_name)
            )
        fin_commands = self._validated_fin_commands(
            configuration, grid_fin_deflections_rad
        )

        migration = properties.cg_rate_body_mps
        relative = state.velocity_body_mps - migration - wind_body_mps
        speed = float(np.linalg.norm(relative))
        reference_area = configuration.reference_area_m2
        if speed <= 1e-10:
            return AeroOutput(
                np.zeros(3),
                0.0,
                0.0,
                0.0,
                0.0,
                float(configuration.axial_drag(0.0)),
                0.0,
                np.zeros(3),
                0.0,
            )
        direction = relative / speed
        cosine = float(direction[0])
        sine = float(np.linalg.norm(direction[1:]))
        total_angle = float(np.arctan2(sine, cosine))
        mach = speed / atmosphere.sound_speed_mps
        axial_coefficient = float(configuration.axial_drag(mach))
        dynamic_pressure = 0.5 * atmosphere.density_kgm3 * speed * speed
        force = np.zeros(3)
        moment = np.zeros(3)

        force[0] -= (
            dynamic_pressure * reference_area * axial_coefficient * cosine * abs(cosine)
        )
        normal_force = 0.0
        cg = properties.cg_station_body_m
        if sine > 0.0:
            side = np.zeros(3)
            side[1:] = -direction[1:] / sine

            def arm(station_m):
                return np.array([station_m, 0.0, 0.0]) - cg

            potential = (
                dynamic_pressure
                * reference_area
                * configuration.potential_normal_slope_per_rad
                * sine
                * abs(cosine)
            )
            force += potential * side

            station = (
                configuration.potential_center_station_m
                if cosine >= 0.0
                else 2.0 * cg[0] - configuration.potential_center_station_m
            )
            station = min(max(station, 0.0), configuration.reference_length_m)
            moment += np.cross(arm(station), potential * side)
            normal_force += potential
            for segment in configuration.side_segments:
                crossflow = (
                    dynamic_pressure
                    * segment.crossflow_drag_coefficient
                    * segment.side_area_m2
                    * sine
                    * sine
                )
                force += crossflow * side
                moment += np.cross(arm(segment.centroid_m), crossflow * side)
                normal_force += crossflow

        damping_scale = (
            dynamic_pressure
            * reference_area
            * configuration.reference_length_m**2
            / (2.0 * speed)
        )
        omega = state.omega_body_radps
        moment += damping_scale * np.array(
            [
                configuration.roll_damping_coefficient * omega[0],
                configuration.pitch_damping_coefficient * omega[1],
                configuration.pitch_damping_coefficient * omega[2],
            ]
        )

        fin_control_force = 0.0
        control = configuration.grid_fin_control
        if control is not None:
            for azimuth, deflection in zip(control.azimuths_rad, fin_commands):
                if deflection == 0.0:
                    continue
                magnitude = (
                    dynamic_pressure
                    * control.fin_area_m2
                    * control.lift_slope_per_rad
                    * deflection
                    * cosine
                    * abs(cosine)
                )
                fin_force = (
                    np.array([0.0, -np.sin(azimuth), np.cos(azimuth)]) * magnitude
                )
                fin_station = np.array(
                    [
                        control.station_m,
                        control.pivot_radius_m * np.cos(azimuth),
                        control.pivot_radius_m * np.sin(azimuth),
                    ]
                )
                force += fin_force
                moment += np.cross(fin_station - cg, fin_force)
                fin_control_force += abs(magnitude)
        return AeroOutput(
            force,
            speed,
            mach,
            dynamic_pressure,
            total_angle,
            axial_coefficient,
            normal_force / (dynamic_pressure * reference_area),
            moment,
            fin_control_force,
        )
