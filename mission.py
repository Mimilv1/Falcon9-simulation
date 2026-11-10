from contextlib import ExitStack
from copy import deepcopy
from dataclasses import asdict, is_dataclass
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import numpy as np

from Aerodynamics import Aerodynamics
from Environement import Environment
from Events import (
    EventConditionError,
    MissionSequence,
    landing_verdict,
    load_mission,
    restart_conditions,
)
from VehicleData import build_falcon9
from ascent import SteadyPoweredAscent, ControlledAscent, PitchProgram
from dynamics import rk4_step
from engine_schedule import SmoothShutdown
from frames import body_to_ecef_matrix
from kalman import NavigationKF, misalignment_terms
from mass_flow import PropellantLimitError
from sensors import (
    GNSS,
    GNSSSpec,
    GyroAttitude,
    GyroAttitudeSpec,
    attitude_error_rotation_vector,
    IMU,
    IMUSpec,
    Lidar,
    LidarSpec,
    StarTracker,
    StarTrackerSpec,
)
from staging import separate_coasting_vehicle

ROOT = Path(__file__).parent


def json_value(value):
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_value(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if callable(value):
        return "constant_at_temp_ref_c"
    return value


def save_json(path, value):
    path.write_text(
        json.dumps(json_value(value), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def vector_row(prefix, vector, axes="xyz"):
    return {f"{prefix}_{axis}": float(value) for axis, value in zip(axes, vector)}


def structural_state(state, properties, station):
    rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
    arm = np.asarray(station) - properties.cg_station_body_m
    position = rotation @ (state.position_body_m + arm)
    velocity = rotation @ (
        state.velocity_body_mps
        + np.cross(state.omega_body_radps, arm)
        - properties.cg_rate_body_mps
    )
    return position, velocity


def prepare_launch(plan, vehicle=None):
    vehicle = build_falcon9() if vehicle is None else vehicle
    if plan.get("vehicle_configuration") == "booster":
        vehicle, _ = vehicle.separate()
        vehicle.stages[0].propellant_mass = plan["terminal_test"][
            "stage1_propellant_kg"
        ]
    initial = plan["initial_state"]
    vehicle.position = np.asarray(initial["structural_base_position_body_m"], float)
    vehicle.velocity = np.asarray(initial["structural_base_velocity_body_mps"], float)
    vehicle.quaternion = np.asarray(initial["quaternion_body_to_ecef_wxyz"], float)
    vehicle.omega = np.asarray(initial["omega_body_radps"], float)
    for stage in vehicle.stages:
        for engine in stage.engines:
            engine.running = stage.name == "stage1" and (
                plan.get("vehicle_configuration") != "booster"
                or engine is stage.engines[0]
            )
            engine.throttle = initial["stage1_throttle"] if engine.running else 0.0
    return vehicle


def default_provider(vehicle, phase, plan):
    if phase not in {
        "powered_ascent",
        "engine_shutdown",
        "full_stack_coast",
        "booster_coast",
        "terminal_braking",
    }:
        raise EventConditionError("Current B models do not implement " + phase + ".")
    environment = Environment(**plan["environment"])
    aerodynamics = Aerodynamics.from_json(ROOT / "data" / "aero_configurations.json")
    envelope = None
    if phase != "terminal_braking":
        end = next(e["time_s"] for e in plan["transitions"] if e["name"] == "meco")
        envelope = SmoothShutdown(
            end - plan["shutdown"]["duration_s"], plan["shutdown"]["duration_s"]
        )
    base = SteadyPoweredAscent(
        vehicle,
        environment,
        aerodynamics,
        command_envelope=envelope,
        aero_configuration=(
            "booster_landing"
            if phase == "terminal_braking"
            else "booster_alone" if phase == "booster_coast" else "full_stack"
        ),
    )
    if phase in {"powered_ascent", "engine_shutdown"}:
        p = plan["attitude_program"]
        return ControlledAscent(
            base,
            PitchProgram(p["start_s"], p["tau_s"], np.deg2rad(p["tilt_max_deg"])),
            p["kp"],
            p["kd"],
            p.get("heading_deg", 90.0),
        )
    return base


def physical_conditions(vehicle, state, loads):
    lower = vehicle.mass_model().propellant_bounds(None)[0]
    return {
        "fuel_available": bool(state.propellant_kg[0] > lower[0]),
        "zero_flow": bool(
            np.all(loads.consumption_kgps == 0)
            and np.all(loads.consumption_rate_kgps2 == 0)
        ),
        "engines_off": all(
            not e.running and e.throttle == 0
            for stage in vehicle.stages
            for e in stage.engines
        ),
        "full_stack_attached": [s.name for s in vehicle.stages] == ["stage1", "stage2"],
        "detached": [s.name for s in vehicle.stages] == ["stage1"],
    }


class NavigationRun:
    def __init__(self, plan, state, properties, seed):
        settings = plan["integration"]
        self.environment = Environment(**plan["environment"])
        self.imu = IMU(
            IMUSpec(
                settings["imu_station_body_m"],
                rate_hz=settings["imu_rate_hz"],
                **settings["imu_errors"],
            ),
            seed=seed,
        )
        self.gnss = GNSS(GNSSSpec(settings["gnss_station_body_m"]), seed=seed + 1)
        self.lidar = Lidar(
            LidarSpec(
                settings["lidar_station_body_m"],
                planet_radius_m=self.environment.earth_radius_m,
            ),
            seed=seed + 2,
        )
        p, v = structural_state(state, properties, self.imu.spec.station_body_m)
        # Open-loop gyro propagation: kept as a dead-reckoning baseline to compare
        # against the MEKF. It starts from the same (possibly misaligned) attitude.
        self.attitude = GyroAttitude(
            GyroAttitudeSpec(
                initial_error_rotvec_rad=settings.get(
                    "initial_attitude_error_rad", (0.0, 0.0, 0.0)
                )
            ),
            state.quaternion_body_to_ecef_wxyz,
        )
        # True : the navigation filter is a MEKF; it estimates attitude and gyro
        #        bias from the noisy gyros + GNSS/lidar and uses its own attitude.
        # False: the filter uses the true attitude (original behaviour).
        self.use_estimated_attitude = settings.get("use_estimated_attitude", True)
        attitude_kwargs = {}
        if self.use_estimated_attitude:
            attitude_kwargs = dict(
                q0=self.attitude.q.copy(),
                attitude_sigma0=settings.get("initial_attitude_sigma_rad", 0.01),
                gyro_bias_sigma0=settings.get("initial_gyro_bias_sigma_radps", 1e-3),
                # Gyro scale factor + misalignment estimated by the MEKF (the truth
                # values live in plan["integration"]["imu_errors"]; filter starts at 0)
                estimate_gyro_calibration=settings.get("estimate_gyro_calibration", True),
                gyro_sf_sigma0=settings.get("initial_gyro_scale_factor_sigma", 1e-3),
                gyro_mis_sigma0=settings.get(
                    "initial_gyro_misalignment_sigma_rad", 1e-3
                ),
                gyro_sf_rw=settings.get("gyro_scale_factor_rw", 0.0),
                gyro_mis_rw=settings.get("gyro_misalignment_rw", 0.0),
            )
        self.filter = NavigationKF.from_imu_spec(
            self.imu.spec,
            p + settings["initial_position_error_m"],
            v,
            settings["initial_sigma"],
            environment=self.environment,
            **attitude_kwargs,
        )
        # Star tracker (absolute attitude aiding of the MEKF). Optional keys in
        # plan["integration"]["star_tracker"]: "enabled" (default True), "gate_nis"
        # (default 16.27 = 99.9 % for 3 dof) and any StarTrackerSpec field.
        st_cfg = dict(settings.get("star_tracker", {}))
        st_enabled = st_cfg.pop("enabled", True) and self.use_estimated_attitude
        self.startracker_gate_nis = st_cfg.pop("gate_nis", 16.27)
        st_cfg.setdefault("planet_radius_m", self.environment.earth_radius_m)
        self.startracker = (
            StarTracker(
                StarTrackerSpec(**st_cfg),
                earth_rotation_ecef_radps=self.environment.earth_rotation_ecef_radps,
                seed=seed + 3,
            )
            if st_enabled
            else None
        )
        self.last_time = None
        self.last_omega = np.zeros(3)
        self.lidar_updates = 0
        self.startracker_updates = 0
        self.startracker_nis = []
        self.nis = []

    @property
    def sensors(self):
        """Sensors that schedule their own sample times."""
        sensors = [self.imu, self.gnss, self.lidar]
        return sensors + ([self.startracker] if self.startracker is not None else [])

    def observe(self, time_s, state, loads, properties):
        truth_rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        imu = self.imu.update(time_s, state, loads, properties)
        if imu is not None:
            self.attitude.update(time_s, imu, state.quaternion_body_to_ecef_wxyz)
        # With the MEKF the filter ignores `rotation` and uses its own attitude
        rotation = (
            self.filter.rotation if self.use_estimated_attitude else truth_rotation
        )
        if imu is not None:
            if self.last_time is not None:
                self.filter.predict(
                    imu.accel_meas_mps2,
                    rotation,
                    time_s - self.last_time,
                    gyro_meas=imu.gyro_meas_radps,
                )
            self.last_time, self.last_omega = time_s, imu.gyro_meas_radps
        gnss = self.gnss.update(time_s, state, properties)
        if gnss is not None and gnss.valid:
            arm = self.gnss.spec.antenna_station_body_m - self.imu.spec.station_body_m
            self.filter.update_gnss(
                gnss.pos_meas_m,
                gnss.vel_meas_mps,
                rotation,
                self.last_omega,
                arm,
                gnss.pos_sigma_m,
                gnss.vel_sigma_mps,
                self.gnss.spec.latency_s,
            )
            self.nis.append(self.filter.last_nis)
        lidar = self.lidar.update(time_s, state, properties)
        if lidar is not None and lidar.valid:
            arm = self.lidar.spec.station_body_m - self.imu.spec.station_body_m
            update = self.filter.update_lidar(
                lidar.range_meas_m,
                rotation,
                arm,
                self.lidar.spec.boresight_body,
                self.lidar.spec.planet_center_ecef_m,
                self.lidar.spec.planet_radius_m,
                lidar.sigma_m,
            )
            self.lidar_updates += int(update is not None)
        star = None
        if self.startracker is not None:
            star = self.startracker.update(time_s, state)
            if star is not None and star.valid:
                residual = self.filter.update_startracker(
                    star.q_meas_sensor_to_eci_wxyz,
                    star.t_meas_s,
                    self.startracker.spec.q_body_from_sensor,
                    star.sigma_sensor_rad,
                    self.environment.earth_rotation_ecef_radps,
                    star.latency_s,
                    gate_nis=self.startracker_gate_nis,
                )
                if residual is not None:
                    self.startracker_updates += 1
                    self.startracker_nis.append(self.filter.last_st_nis)
        return imu, gnss, lidar, star

    def row(self, time_s, state, properties):
        p, v = structural_state(state, properties, self.imu.spec.station_body_m)
        x = self.filter.x
        q_true = state.quaternion_body_to_ecef_wxyz
        gyro_only_error = attitude_error_rotation_vector(self.attitude.q, q_true)
        if self.use_estimated_attitude:
            q_est = self.filter.q
            attitude_sigma = self.filter.attitude_sigma_rad
        else:
            q_est = q_true
            attitude_sigma = np.full(3, np.nan)
        attitude_error = attitude_error_rotation_vector(q_est, q_true)
        return {
            "time_s": time_s,
            # Attitude used by the navigation filter (MEKF estimate, or truth)
            **vector_row("quaternion_est_body_to_ecef", q_est, "wxyz"),
            **vector_row("attitude_error_body_rad", attitude_error),
            "attitude_error_deg": float(np.rad2deg(np.linalg.norm(attitude_error))),
            **vector_row("attitude_sigma_body_rad", attitude_sigma),
            **vector_row("gyro_bias_est_body_radps", self.filter.bg),
            **vector_row("gyro_scale_factor_est", self.filter.gyro_sf),
            **vector_row(
                "gyro_misalignment_est_rad",
                self.filter.gyro_mis,
                ("xy", "xz", "yx", "yz", "zx", "zy"),
            ),
            # Open-loop gyro propagation baseline (no aiding)
            "gyro_only_attitude_error_deg": float(
                np.rad2deg(np.linalg.norm(gyro_only_error))
            ),
            **vector_row("position_ecef_m", x[:3]),
            **vector_row("velocity_ecef_mps", x[3:6]),
            **vector_row("accel_bias_body_mps2", x[6:9]),
            **vector_row("position_error_ecef_m", x[:3] - p),
            **vector_row("velocity_error_ecef_mps", x[3:6] - v),
            **{
                f"sigma_{i}": float(s)
                for i, s in enumerate(np.sqrt(np.diag(self.filter.P)))
            },
        }


def measurement_rows(samples):
    imu, gnss, lidar, star = samples
    rows = []
    if imu is not None:
        rows.append(
            {
                "time_s": imu.t,
                "sensor": "imu",
                "valid": True,
                **vector_row("accel_body_mps2", imu.accel_meas_mps2),
                **vector_row("gyro_body_radps", imu.gyro_meas_radps),
            }
        )
    if gnss is not None:
        rows.append(
            {
                "time_s": gnss.t,
                "sensor": "gnss",
                "valid": gnss.valid,
                **vector_row("position_ecef_m", gnss.pos_meas_m),
                **vector_row("velocity_ecef_mps", gnss.vel_meas_mps),
            }
        )
    if lidar is not None:
        rows.append(
            {
                "time_s": lidar.t,
                "sensor": "lidar",
                "valid": lidar.valid,
                "range_m": lidar.range_meas_m,
                "range_sigma_m": lidar.sigma_m,
            }
        )
    if star is not None:
        rows.append(
            {
                "time_s": star.t,
                "sensor": "startracker",
                "valid": star.valid,
                "startracker_reason": star.reason,
                **vector_row(
                    "quaternion_meas_sensor_to_eci",
                    star.q_meas_sensor_to_eci_wxyz,
                    "wxyz",
                ),
            }
        )
    return rows


def truth_row(time_s, phase, state, evaluation):
    properties = evaluation.properties
    rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
    base_position, base_velocity = structural_state(state, properties, [0, 0, 0])
    return {
        "time_s": time_s,
        "phase": phase,
        **vector_row("cg_position_body_m", state.position_body_m),
        **vector_row("cg_velocity_body_mps", state.velocity_body_mps),
        **vector_row(
            "quaternion_body_to_ecef",
            state.quaternion_body_to_ecef_wxyz,
            "wxyz",
        ),
        **vector_row("omega_body_radps", state.omega_body_radps),
        "stage1_propellant_kg": float(state.propellant_kg[0]),
        "stage2_propellant_kg": float(state.propellant_kg[1]),
        **vector_row("cg_position_ecef_m", rotation @ state.position_body_m),
        **vector_row("cg_velocity_ecef_mps", rotation @ state.velocity_body_mps),
        **vector_row("force_body_n", evaluation.loads.force_body_n),
        **vector_row("moment_body_nm", evaluation.loads.moment_body_nm),
        **vector_row("gravity_body_mps2", evaluation.loads.gravity_body_mps2),
        **vector_row(
            "frame_acceleration_body_mps2",
            evaluation.loads.frame_acceleration_body_mps2,
        ),
        "base_vertical_velocity_mps": float(
            np.dot(base_position, base_velocity) / np.linalg.norm(base_position)
        ),
        "speed_mps": float(np.linalg.norm(state.velocity_body_mps)),
        "mach": evaluation.aero.mach,
        "mass_kg": properties.mass_kg,
        "cg_station_m": float(properties.cg_station_body_m[0]),
        **vector_row("inertia_kgm2", np.diag(properties.inertia_body_kgm2)),
        "base_altitude_m": evaluation.base_altitude_m,
        "cg_altitude_m": evaluation.cg_altitude_m,
        "thrust_n": evaluation.thrust_n,
        "dynamic_pressure_pa": evaluation.aero.dynamic_pressure_pa,
        "angle_of_attack_rad": evaluation.aero.angle_of_attack_rad,
    }


def run_mission(
    output,
    *,
    until_s=160.0,
    step_s=None,
    seed=None,
    navigation=True,
    plan=None,
    provider_factory=default_provider,
    conditions_provider=None,
    recovery=False,
):
    plan = load_mission() if plan is None else plan
    controller = None
    if recovery:
        from recovery import RecoveryController, recovery_plan

        plan = recovery_plan(plan)
        controller = RecoveryController(plan)
        provider_factory = controller.provider
        conditions_provider = controller.conditions
    settings = plan["integration"]
    step_s = settings["step_s"] if step_s is None else step_s
    seed = settings["seed"] if seed is None else seed
    if (
        not np.isfinite(until_s)
        or until_s <= 0
        or not np.isfinite(step_s)
        or step_s <= 0
    ):
        raise ValueError("Require positive finite duration and step.")
    if not np.isfinite(settings["imu_rate_hz"]) or settings["imu_rate_hz"] <= 0:
        raise ValueError("IMU rate must be positive and finite.")
    if navigation and step_s > min(1 / settings["imu_rate_hz"], 0.1) + 1e-12:
        raise ValueError("Step must not exceed the configured IMU period.")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    names = (
        "truth.csv",
        "measurements.csv",
        "navigation.csv",
        "events.json",
        "summary.json",
    )
    if any((output / name).exists() for name in names):
        raise FileExistsError("Choose an output directory without existing run files.")
    input_hashes = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (
            *ROOT.glob("*.py"),
            *ROOT.glob("data/*.json"),
            *ROOT.glob("data/*.csv"),
        )
    }
    vehicle = prepare_launch(plan)
    upper = None
    sequence = MissionSequence(plan)
    provider = provider_factory(vehicle, sequence.phase, plan)
    flow, derivative = provider.flow_at(0.0)
    state = vehicle.to_state(
        consumption_by_stage={"stage1": flow[0]},
        consumption_derivative_by_stage={"stage1": derivative[0]},
    )
    evaluation = provider.evaluate(0.0, state)
    nav = (
        NavigationRun(plan, state, evaluation.properties, seed) if navigation else None
    )
    time_s, next_output, next_progress = 0.0, 0.0, 10.0
    status, detail = "requested_interval_complete", ""
    truth_rows, nav_rows, sensor_counts = (
        [],
        [],
        {"imu": 0, "gnss": 0, "lidar": 0, "startracker": 0},
    )
    with ExitStack() as stack:

        def writer(name, fields):
            stream = stack.enter_context(
                (output / name).open("w", newline="", encoding="utf-8")
            )
            result = csv.DictWriter(stream, fieldnames=fields)
            result.writeheader()
            return result

        tw = writer("truth.csv", list(truth_row(0, sequence.phase, state, evaluation)))
        fields = ["time_s", "sensor", "valid", "range_m", "range_sigma_m"]
        for prefix in (
            "accel_body_mps2",
            "gyro_body_radps",
            "position_ecef_m",
            "velocity_ecef_mps",
        ):
            fields.extend(f"{prefix}_{axis}" for axis in "xyz")
        fields.append("startracker_reason")
        fields.extend(f"quaternion_meas_sensor_to_eci_{a}" for a in "wxyz")
        sw = writer("measurements.csv", fields)
        nw = (
            writer("navigation.csv", list(nav.row(0, state, evaluation.properties)))
            if nav
            else None
        )
        while True:
            conditions = (
                {}
                if conditions_provider is None
                else dict(
                    conditions_provider(time_s, vehicle, state, evaluation, sequence)
                )
            )
            conditions.update(physical_conditions(vehicle, state, evaluation.loads))
            event = sequence.due(time_s, conditions)
            clock = sequence.next_clock_time
            if (
                sequence.next_event is not None
                and sequence.next_event.time_s is None
                and clock is not None
            ):
                if time_s >= clock - 1e-9 and event is None:
                    status, detail = (
                        "unresolved_event",
                        sequence.next_event.name + " was not completed.",
                    )
                    break
            if event is not None:
                try:
                    if event.name in {
                        "boostback_start",
                        "entry_start",
                        "landing_start",
                    }:
                        position, velocity = structural_state(
                            state, evaluation.properties, [0, 0, 0]
                        )
                        descending = np.dot(position, velocity) < 0
                        conditions.update(
                            restart_conditions(
                                vehicle,
                                state,
                                evaluation.loads,
                                required_propellant_kg=conditions.get(
                                    "required_propellant_kg"
                                ),
                                attitude_ready=conditions.get("attitude_ready", False),
                                propulsion_ready=conditions.get(
                                    "propulsion_ready", False
                                ),
                                descending=descending,
                            )
                        )
                    if event.name == "touchdown":
                        landing = plan["landing_criteria"]
                        if (
                            landing["contact_height_m"] is None
                            or landing["landing_target_ecef_m"] is None
                        ):
                            raise EventConditionError(
                                "Define the landing target and leg-contact height first."
                            )
                        metrics = dict(conditions.get("landing_metrics", {}))
                        metrics["ground_contact"] = conditions.get(
                            "ground_contact", False
                        )
                        conditions.update(landing_verdict(metrics, landing["limits"]))
                    record = sequence.accept(event.name, time_s, conditions)
                    record["published_reference_time_s"] = next(
                        e.get("reference_time_s", e.get("time_s"))
                        for e in plan["transitions"]
                        if e["name"] == event.name
                    )
                    if event.name == "meco":
                        for stage in vehicle.stages:
                            for engine in stage.engines:
                                engine.running, engine.throttle = False, 0.0
                        vehicle.apply_state(state, time_s)
                    elif event.name == "stage_separation":
                        record["parent_state"] = state.to_vector().tolist()
                        split = separate_coasting_vehicle(
                            vehicle, state, time_s, evaluation.loads
                        )
                        vehicle, state = split.lower, split.lower_state
                        upper = (split.upper, split.upper_state, time_s)
                        record["conservation_audit"] = split.audit
                        record["booster_state"] = state.to_vector().tolist()
                        record["upper_state"] = split.upper_state.to_vector().tolist()
                    if event.name == "touchdown":
                        record["landing_metrics"] = json_value(metrics)
                        status = (
                            "landing_success"
                            if conditions["landing_success"]
                            else "landing_failed"
                        )
                    else:
                        provider = provider_factory(vehicle, sequence.phase, plan)
                        evaluation = provider.evaluate(time_s, state)
                except EventConditionError as error:
                    status, detail = "event_conditions_unmet", str(error)
                    break
            samples = (
                nav.observe(time_s, state, evaluation.loads, evaluation.properties)
                if nav
                else None
            )
            if samples:
                for row in measurement_rows(samples):
                    sw.writerow(row)
                    sensor_counts[row["sensor"]] += 1
            if time_s >= next_output - 1e-9 or time_s >= until_s - 1e-9:
                row = truth_row(time_s, sequence.phase, state, evaluation)
                tw.writerow(row)
                truth_rows.append(row)
                if nav:
                    row = nav.row(time_s, state, evaluation.properties)
                    nw.writerow(row)
                    nav_rows.append(row)
                next_output += settings["output_period_s"]
            if time_s >= until_s - 1e-9 or sequence.phase in {
                "landed",
                "landing_failed",
            }:
                break
            clock = sequence.next_clock_time
            boundary = until_s if clock is None else min(until_s, clock)
            h = min(step_s, boundary - time_s)
            if controller is not None and controller.phase is not None:
                for knot in (
                    controller.next_update,
                    controller.ramp[0] + controller.ramp[1],
                ):
                    if knot > time_s + 1e-9:
                        h = min(h, knot - time_s)
            if next_output > time_s + 1e-9:
                h = min(h, next_output - time_s)
            if controller is not None and sequence.phase == "landing_burn":
                clearance = controller.contact_altitude(state, evaluation.properties)
                _, velocity = controller.point_state(state, evaluation.properties)
                h = min(
                    h,
                    max(
                        1e-8, 0.8 * max(clearance, 0) / max(np.linalg.norm(velocity), 1)
                    ),
                )
            if nav:
                future = [
                    sensor.next_t
                    for sensor in nav.sensors
                    if sensor.next_t > time_s + 1e-9
                ]
                h = min(h, min(future) - time_s)
            if h <= 1e-12:
                raise RuntimeError("Mission did not advance past an event boundary.")
            try:
                candidate = rk4_step(
                    time_s, state, h, vehicle.mass_model(), None, provider
                )
                candidate_time = min(time_s + h, until_s)
                if abs(candidate_time - boundary) < 1e-10:
                    candidate_time = boundary
                candidate_evaluation = provider.evaluate(candidate_time, candidate)
                if (
                    controller is not None
                    and sequence.phase == "landing_burn"
                    and controller.contact_altitude(
                        candidate, candidate_evaluation.properties
                    )
                    <= 0
                ):
                    lo, hi = 0.0, h
                    for _ in range(35):
                        mid = (lo + hi) / 2
                        trial = rk4_step(
                            time_s, state, mid, vehicle.mass_model(), None, provider
                        )
                        trial_ev = provider.evaluate(time_s + mid, trial)
                        if controller.contact_altitude(trial, trial_ev.properties) > 0:
                            lo = mid
                        else:
                            hi = mid
                        if hi - lo < 1e-9:
                            break
                    candidate_time = time_s + hi
                    candidate = rk4_step(
                        time_s, state, hi, vehicle.mass_model(), None, provider
                    )
                    candidate_evaluation = provider.evaluate(candidate_time, candidate)
            except PropellantLimitError as error:
                status, detail = "propellant_limit", str(error)
                break
            except ValueError as error:
                if not any(
                    text in str(error)
                    for text in (
                        "supports altitude",
                        "outside the adopted drag curve",
                        "outside the configuration",
                    )
                ):
                    raise
                status, detail = "model_domain_limit", str(error)
                break
            state, time_s, evaluation = candidate, candidate_time, candidate_evaluation
            if time_s >= next_progress - 1e-9:
                print(
                    f"  {time_s:7.2f} s | {sequence.phase:18s} | h={evaluation.cg_altitude_m/1000:6.2f} km",
                    flush=True,
                )
                next_progress += 10.0
            vehicle.apply_state(
                state,
                time_s,
                consumption_by_stage=dict(
                    zip(("stage1", "stage2"), evaluation.loads.consumption_kgps)
                ),
                consumption_derivative_by_stage=dict(
                    zip(("stage1", "stage2"), evaluation.loads.consumption_rate_kgps2)
                ),
            )
        last = truth_row(time_s, sequence.phase, state, evaluation)
        if not truth_rows or abs(truth_rows[-1]["time_s"] - time_s) > 1e-9:
            tw.writerow(last)
            truth_rows.append(last)
            if nav:
                row = nav.row(time_s, state, evaluation.properties)
                nw.writerow(row)
                nav_rows.append(row)
    if (
        recovery
        and status == "requested_interval_complete"
        and sequence.phase != "landed"
    ):
        status, detail = (
            "mission_time_limit",
            "Return not completed by requested simulation endpoint.",
        )
    summary = {
        "mission": plan["mission"],
        "status": status,
        "detail": detail,
        "requested_end_s": until_s,
        "actual_end_s": time_s,
        "phase": sequence.phase,
        "final_state_vector": state.to_vector().tolist(),
        "full_mission_completed": sequence.phase == "landed",
        "step_s": step_s,
        "seed": seed,
        "sensor_counts": sensor_counts,
        "initial_truth": truth_rows[0],
        "final_truth": last,
        "plan": plan,
        "limitations": [
            "B multilayer atmosphere, inverse-square gravity, Earth-rotation accelerations and multi-configuration aero are active.",
            "Generic fairing-based full-stack Cd used for Dragon as A; high-Mach extension and booster coefficients estimated. Mach >10 rejected.",
            "Ideal pitch PD adds an unlimited pure couple; no physical TVC or booster flip. Return burns and touchdown are not implemented.",
            "Liquid geometry remains a settled-fluid approximation during coast.",
            "Navigation uses truth attitude; sensor specifications remain placeholders.",
            "Published planned times are not actual observed event times.",
        ],
        "input_sha256": input_hashes,
    }
    if controller is not None:
        if sequence.phase in {"landed", "landing_failed"}:
            controller.burns.append(
                {
                    "phase": "landing_burn",
                    "start_s": controller.started,
                    "end_s": time_s,
                    "consumed_kg": float(
                        controller.phase_fuel - state.propellant_kg[0]
                    ),
                    "reason": "contact",
                    "contact_passed": sequence.phase == "landed",
                }
            )
        summary["recovery_burns"] = controller.burns
        summary["limitations"][
            2
        ] = "Project-designed recovery using true-state feedback and ideal bounded moments; not physical TVC/RCS or SpaceX flight software. Outcomes are measured, not forced."
    summary["physical_state_finite"] = bool(np.all(np.isfinite(state.to_vector())))
    summary["max_dynamic_pressure_pa"] = max(
        r["dynamic_pressure_pa"] for r in truth_rows
    )
    summary["max_dynamic_pressure_time_s"] = max(
        truth_rows, key=lambda r: r["dynamic_pressure_pa"]
    )["time_s"]
    summary["next_event"] = sequence.next_event.name if sequence.next_event else None
    if nav:
        summary["lidar_updates"] = nav.lidar_updates
        summary["gnss_nis_mean"] = float(np.mean(nav.nis)) if nav.nis else None
        summary["sensor_specs"] = {
            "imu": nav.imu.spec,
            "gnss": nav.gnss.spec,
            "lidar": nav.lidar.spec,
            **(
                {"startracker": nav.startracker.spec}
                if nav.startracker is not None
                else {}
            ),
        }
        summary["startracker_updates"] = nav.startracker_updates
        summary["startracker_rejected"] = (
            nav.filter.startracker_rejected if nav.startracker is not None else 0
        )
        summary["startracker_nis_mean"] = (
            float(np.mean(nav.startracker_nis)) if nav.startracker_nis else None
        )
        errors = np.array(
            [[r[f"position_error_ecef_m_{a}"] for a in "xyz"] for r in nav_rows]
        )
        summary["navigation_position_rms_m"] = float(
            np.sqrt(np.mean(np.sum(errors**2, axis=1)))
        )
        summary["attitude_source"] = (
            "mekf_gyro_gnss_lidar" if nav.use_estimated_attitude else "truth"
        )
        summary["limitations"][4] = (
            "Navigation attitude and gyro bias are estimated by a MEKF from noisy gyros "
            "with GNSS/lidar aiding (roll about the thrust axis is weakly observable); "
            "sensor specifications remain placeholders."
            if nav.use_estimated_attitude
            else "Navigation uses truth attitude; sensor specifications remain placeholders."
        )
        summary["attitude_error_final_deg"] = nav_rows[-1]["attitude_error_deg"]
        summary["attitude_error_rms_deg"] = float(
            np.sqrt(np.mean([r["attitude_error_deg"] ** 2 for r in nav_rows]))
        )
        summary["gyro_only_attitude_error_final_deg"] = nav_rows[-1][
            "gyro_only_attitude_error_deg"
        ]
        summary["gyro_bias_estimate_final_radps"] = [
            nav_rows[-1][f"gyro_bias_est_body_radps_{a}"] for a in "xyz"
        ]
        if nav.filter.calibration_enabled:
            pairs = ("xy", "xz", "yx", "yz", "zx", "zy")
            summary["gyro_scale_factor_estimate_final"] = [
                nav_rows[-1][f"gyro_scale_factor_est_{a}"] for a in "xyz"
            ]
            summary["gyro_scale_factor_true"] = nav.imu.spec.gyro_scale_factor
            summary["gyro_scale_factor_sigma_final"] = nav.filter.gyro_sf_sigma
            summary["gyro_misalignment_estimate_final_rad"] = [
                nav_rows[-1][f"gyro_misalignment_est_rad_{p}"] for p in pairs
            ]
            summary["gyro_misalignment_true_rad"] = misalignment_terms(
                nav.imu.spec.gyro_misalignment
            )
            summary["gyro_misalignment_sigma_final_rad"] = nav.filter.gyro_mis_sigma_rad
    summary["physics_case"] = (
        "independent_terminal_test"
        if plan.get("vehicle_configuration") == "booster"
        else "continuous_project_recovery" if recovery else "ascent_to_booster_coast"
    )
    if upper is not None:
        from upper_stage import propagate_upper_stage

        print("Propagating separated stage 2...", flush=True)
        summary["stage2"] = propagate_upper_stage(*upper, until_s, plan, output)
        save_json(output / "stage2_events.json", summary["stage2"]["events"])
    summary["reference_comparison"] = (
        compare_reference(output, truth_rows)
        if plan.get("vehicle_configuration") != "booster"
        else {
            "count": 0,
            "status": "Not applicable: independent assumed terminal-descent case.",
        }
    )
    try:
        summary["git_commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        summary["git_commit"] = None
    save_json(output / "events.json", sequence.log)
    save_json(output / "summary.json", summary)
    return json_value(summary)


def compare_reference(output, truth_rows):
    path = ROOT / "data" / "crs10_ascent_reference.csv"
    with path.open(encoding="utf-8-sig") as stream:
        reference = list(csv.DictReader(stream))
    times = np.array([r["time_s"] for r in truth_rows])
    rows = []
    for r in reference:
        t = float(r["time_s"])
        if times[0] <= t <= times[-1]:
            item = {"time_s": t}
            for name, field in (
                ("altitude_m", "cg_altitude_m"),
                ("speed_mps", "speed_mps"),
            ):
                value = float(r[name])
                simulated = float(np.interp(t, times, [x[field] for x in truth_rows]))
                item.update(
                    {
                        "reference_" + name: value,
                        "simulated_" + name: simulated,
                        "error_" + name: simulated - value,
                    }
                )
            rows.append(item)
    if not rows:
        return {"count": 0, "status": "no overlapping reference samples"}
    with (output / "validation.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "count": len(rows),
        "interval_s": [rows[0]["time_s"], rows[-1]["time_s"]],
        "basis": "Video-digitized ascent telemetry, E. Speed reference frame not independently verified. Not first-stage return data.",
        "final_sample": rows[-1],
    }
    for name in ("altitude_m", "speed_mps"):
        result["rmse_" + name] = float(
            np.sqrt(np.mean([r["error_" + name] ** 2 for r in rows]))
        )
    return result


def plot_results(output, *, show=False, save=None):
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)

    def read(name):
        with (output / name).open(encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    truth = read("truth.csv")
    nav = read("navigation.csv") if (output / "navigation.csv").exists() else []
    events = json.loads((output / "events.json").read_text())
    validation = read("validation.csv") if (output / "validation.csv").exists() else []
    stage2 = read("stage2_truth.csv") if (output / "stage2_truth.csv").exists() else []
    t = np.array([float(r["time_s"]) for r in truth])
    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True)
    for ax, field, label, scale in (
        (axes[0, 0], "cg_altitude_m", "CG altitude [km]", 0.001),
        (axes[0, 1], "speed_mps", "Ground-relative CG speed [m/s]", 1),
        (axes[1, 0], "stage1_propellant_kg", "Stage-1 propellant [t]", 0.001),
        (axes[1, 1], "cg_station_m", "CG above structural base [m]", 1),
    ):
        ax.plot(t, [float(r[field]) * scale for r in truth], label="Stack / stage 1")
        ax.set_ylabel(label)
    for ax, field, scale in (
        (axes[0, 0], "altitude_m", 0.001),
        (axes[0, 1], "speed_mps", 1),
    ):
        if stage2:
            ax.plot(
                [float(r["time_s"]) for r in stage2],
                [
                    float(r["cg_" + field] if field == "altitude_m" else r[field])
                    * scale
                    for r in stage2
                ],
                color="tab:purple",
                label="Stage 2 + Dragon (A)",
            )
        if validation:
            ax.plot(
                [float(r["time_s"]) for r in validation],
                [float(r["reference_" + field]) * scale for r in validation],
                "--",
                label="Video readout (E)",
            )
        ax.legend()
    if nav:
        nt = [float(r["time_s"]) for r in nav]
        for ax, prefix, label in (
            (axes[2, 0], "position_error_ecef_m", "Navigation position error [m]"),
            (axes[2, 1], "velocity_error_ecef_mps", "Navigation velocity error [m/s]"),
        ):
            error = np.array([[float(r[prefix + "_" + a]) for a in "xyz"] for r in nav])
            ax.plot(nt, np.linalg.norm(error, axis=1))
            ax.set_ylabel(label)
    else:
        for ax in axes[2]:
            ax.text(
                0.5, 0.5, "Navigation disabled", ha="center", transform=ax.transAxes
            )
    for ax in axes.flat:
        ax.grid(alpha=0.3)
        for event in events[1:]:
            ax.axvline(event["time_s"], color=".6", lw=0.7, ls=":")
        if stage2:
            for event in json.loads((output / "stage2_events.json").read_text())[1:]:
                ax.axvline(
                    event["time_s"], color="tab:purple", lw=0.7, ls="--", alpha=0.5
                )
    for ax in axes[2]:
        ax.set_xlabel("Simulation time [s]")
    fig.suptitle(
        "Team integration: "
        + truth[0]["phase"]
        + " to "
        + truth[-1]["phase"]
        + "\nAssumed commands; navigation attitude source: see summary.json (attitude_source)"
    )
    fig.tight_layout()
    destination = Path(save) if save else output / "overview.png"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    if stage2:
        plot_trajectories(output, truth, stage2, events)
    return destination


def plot_trajectories(output, truth, stage2, events):
    import matplotlib.pyplot as plt

    radius = 6371000.0
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for rows, label in ((truth, "Stack / stage 1"), (stage2, "Stage 2 + Dragon (A)")):
        p = np.array(
            [[float(r["cg_position_ecef_m_" + a]) for a in "xyz"] for r in rows]
        )
        r = np.linalg.norm(p, axis=1)
        downrange = radius * np.arccos(np.clip(p[:, 2] / r, -1, 1)) / 1000
        axes[0].plot(downrange, (r - radius) / 1000, label=label)
        axes[1].plot(-p[:, 1] / 1000, p[:, 0] / 1000, label=label)
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend()
    axes[0].set(
        xlabel="Surface angular distance from launch [km]", ylabel="CG altitude [km]"
    )
    axes[1].set(
        xlabel="External east component [km]", ylabel="External north component [km]"
    )
    fig.suptitle("Separated trajectories; assumed stage-2 commands")
    fig.tight_layout()
    fig.savefig(output / "trajectories.png", dpi=150)
    plt.close(fig)
    upper_events = json.loads((output / "stage2_events.json").read_text())
    names = {
        "liftoff": "Liftoff",
        "shutdown_start": "Shutdown",
        "meco": "MECO",
        "stage_separation": "Separation",
        "boostback_start": "Boostback ignition",
        "boostback_end": "Boostback end",
        "entry_start": "Entry ignition",
        "entry_end": "Entry end",
        "landing_start": "Landing ignition",
        "touchdown": "Touchdown",
        "stage2_ignition": "Stage-2 ignition",
        "seco": "SECO",
    }
    timeline = [
        {
            "vehicle": "Stack / stage 1",
            "time_s": e["time_s"],
            "event": names.get(e["event"], e["event"]),
        }
        for e in events
    ]
    timeline.extend(
        {"vehicle": "Stage 2", "time_s": e["time_s"], "event": names[e["event"]]}
        for e in upper_events[1:]
    )
    timeline.sort(key=lambda e: e["time_s"])
    with (output / "event_timeline.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=["vehicle", "time_s", "event"])
        writer.writeheader()
        writer.writerows(timeline)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis("off")
    table = ax.table(
        cellText=[[f'{e["time_s"]:.1f}', e["vehicle"], e["event"]] for e in timeline],
        colLabels=["Time [s]", "Vehicle", "Event"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.65)
    ax.set_title("Event times for overview.png; simulated events / adopted commands")
    fig.tight_layout()
    fig.savefig(output / "flight_events.png", dpi=150)
    plt.close(fig)


def terminal_plan(plan=None):
    plan = deepcopy(load_mission() if plan is None else plan)
    case = plan["terminal_test"]
    plan["initial_phase"] = "terminal_braking"
    plan["vehicle_configuration"] = "booster"
    plan["scope"] = case["basis"]
    plan["transitions"] = []
    initial = plan["initial_state"]
    initial["structural_base_position_body_m"][0] += case["base_altitude_m"]
    initial["structural_base_velocity_body_mps"] = [
        case["base_vertical_velocity_mps"],
        0,
        0,
    ]
    initial["stage1_throttle"] = case["central_engine_throttle"]
    return plan