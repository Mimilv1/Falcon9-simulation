import csv
import numpy as np

from Environement import Environment
from dynamics import rk4_step
from frames import body_to_ecef_matrix
from load_interface import AppliedForce, assemble_loads
from propulsion import evaluate_steady_engine, scale_shutdown_output


class UpperStageProvider:
    def __init__(self, vehicle, state, time_s, plan):
        self.vehicle = vehicle
        self.model = vehicle.mass_model()
        self.cfg = plan["upper_stage"]
        self.environment = Environment(**plan["environment"])
        self.stage = vehicle.stages[0]
        self.engine = self.stage.engines[0]
        self.engine.running = True
        self.engine.throttle = self.cfg["throttle"]
        self.start = time_s
        rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        up = rotation @ state.position_body_m
        up /= np.linalg.norm(up)
        self.heading = rotation[:, 0] - np.dot(rotation[:, 0], up) * up
        if np.linalg.norm(self.heading) < 1e-10:
            self.heading = np.array([0.0, -1.0, 0.0])
            self.heading -= np.dot(self.heading, up) * up
        self.heading /= np.linalg.norm(self.heading)
        self.initial_pitch = np.arcsin(np.clip(np.dot(rotation[:, 0], up), -1, 1))

    def envelope(self, time_s):
        start, end, ramp = (self.cfg[k] for k in ("ignition_s", "cutoff_s", "ramp_s"))
        if time_s <= start or time_s >= end:
            return 0.0, 0.0
        u = min((time_s - start) / ramp, (end - time_s) / ramp, 1.0)
        value = 10 * u**3 - 15 * u**4 + 6 * u**5
        rate = 30 * u**2 * (1 - u) ** 2 / ramp
        return value, (
            rate if time_s < start + ramp else -rate if time_s > end - ramp else 0.0
        )

    def evaluate(self, time_s, state):
        air = self.environment.atmosphere(self.environment.altitude(state, None))
        output = evaluate_steady_engine(
            self.engine,
            air.pressure_pa,
            stage_index=1,
            stage_offset_m=0.0,
            nozzle_exit_axial_m=0.0,
            g0=self.vehicle.parameters.G0,
        )
        output = scale_shutdown_output(output, *self.envelope(time_s))
        flow = np.array([0.0, output.mass_flow_kgps])
        derivative = np.array([0.0, output.exhaust_port.mass_flow_rate_kgps2])
        props = self.model.mass_properties(
            None, state.propellant_kg, flow, consumption_rate_kgps2=derivative
        )
        env = self.environment
        rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        position = rotation @ state.position_body_m
        up = position / np.linalg.norm(position)
        height = env.altitude(state, props)
        air = env.atmosphere(height)
        relative = (
            state.velocity_body_mps - props.cg_rate_body_mps - env.wind_body(state)
        )
        speed = np.linalg.norm(relative)
        area = np.pi * (self.stage.diameter / 2) ** 2
        drag = (
            -0.5
            * air.density_kgm3
            * area
            * self.cfg["drag_coefficient"]
            * speed
            * relative
        )
        heading = self.heading - np.dot(self.heading, up) * up
        heading /= np.linalg.norm(heading)
        velocity = rotation @ state.velocity_body_mps
        radial_speed = np.dot(velocity, up)
        horizontal_speed2 = np.dot(velocity, velocity) - radial_speed**2
        wanted_speed = np.clip(
            (self.cfg["target_altitude_m"] - height) / self.cfg["height_time_s"],
            -100,
            600,
        )
        radial_acceleration = (wanted_speed - radial_speed) / self.cfg[
            "velocity_time_s"
        ]
        background = env.gravity_ecef(position) + env.frame_acceleration_ecef(
            position, velocity
        )
        required = radial_acceleration - horizontal_speed2 / np.linalg.norm(position)
        required -= np.dot(background + rotation @ drag / props.mass_kg, up)
        pitch = (
            np.arcsin(
                np.clip(required / max(output.thrust_n / props.mass_kg, 1e-6), 0, 1)
            )
            if output.thrust_n > 0
            else self.initial_pitch
        )
        wanted = np.cos(pitch) * heading + np.sin(pitch) * up
        error = rotation.T @ np.cross(rotation[:, 0], wanted)
        alpha = (
            self.cfg["attitude_kp"] * error
            - self.cfg["attitude_kd"] * state.omega_body_radps
        )
        moment = props.inertia_body_kgm2 @ alpha
        loads = assemble_loads(
            props,
            [output.applied_force, AppliedForce(drag, props.cg_station_body_m)],
            gravity_body_mps2=env.gravity_body(state),
            free_moment_body_nm=moment,
            consumption_kgps=flow,
            exhaust_ports=(output.exhaust_port,),
            frame_acceleration_body_mps2=env.frame_acceleration_body(state),
        )
        return loads, props, output.thrust_n, height

    def __call__(self, time_s, state):
        return self.evaluate(time_s, state)[0]


def propagate_upper_stage(vehicle, state, start_s, end_s, plan, output):
    provider = UpperStageProvider(vehicle, state, start_s, plan)
    cfg = provider.cfg
    rows = []
    events = [{"event": "stage_separation", "time_s": start_s}]
    time_s = start_s
    knots = [
        cfg["ignition_s"],
        cfg["ignition_s"] + cfg["ramp_s"],
        cfg["cutoff_s"] - cfg["ramp_s"],
        cfg["cutoff_s"],
        end_s,
    ]
    next_output = start_s
    status = "requested_interval_complete"
    while True:
        loads, props, thrust, height = provider.evaluate(time_s, state)
        if time_s >= next_output - 1e-9 or time_s >= end_s - 1e-9:
            rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
            p = rotation @ state.position_body_m
            row = {
                "time_s": time_s,
                "phase": (
                    "stage2_coast"
                    if time_s < cfg["ignition_s"]
                    else (
                        "stage2_burn"
                        if time_s < cfg["cutoff_s"]
                        else "stage2_post_seco_coast"
                    )
                ),
                "cg_altitude_m": height,
                "speed_mps": float(np.linalg.norm(state.velocity_body_mps)),
                "stage2_propellant_kg": state.propellant_kg[1],
                "mass_kg": props.mass_kg,
                "thrust_n": thrust,
                "cg_station_m": props.cg_station_body_m[0],
            }
            row.update({"cg_position_ecef_m_" + a: p[i] for i, a in enumerate("xyz")})
            rows.append(row)
            next_output = time_s + plan["integration"]["output_period_s"]
        if time_s >= end_s - 1e-9:
            break
        future = [k for k in knots + [next_output] if k > time_s + 1e-9]
        step = min(cfg["step_s"], min(future) - time_s)
        state = rk4_step(time_s, state, step, provider.model, None, provider)
        time_s += step
    for key, name in (("ignition_s", "stage2_ignition"), ("cutoff_s", "seco")):
        if start_s <= cfg[key] <= time_s:
            events.append({"event": name, "time_s": cfg[key]})
    with (output / "stage2_truth.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return {
        "status": status,
        "actual_end_s": time_s,
        "initial": rows[0],
        "final": rows[-1],
        "events": events,
        "settings": cfg,
        "finite_state": bool(np.all(np.isfinite(state.to_vector()))),
    }
