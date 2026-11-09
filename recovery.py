from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import numpy as np

from Aerodynamics import Aerodynamics
from Environement import Environment
from ascent import AscentEvaluation
from frames import body_to_ecef_matrix
from load_interface import AppliedForce, assemble_loads
from mass_flow import ExhaustPort
from propulsion import evaluate_steady_engine, steady_mass_flow


def smooth_command(t, start, duration, y0, d0, y1):
    u = np.clip((t - start) / duration, 0, 1)
    if t >= start + duration:
        return float(y1), 0.0
    if t <= start:
        return float(y0), float(d0)
    a0 = y0
    a1 = d0 * duration
    a3 = 10 * (y1 - y0) - 6 * a1
    a4 = -15 * (y1 - y0) + 8 * a1
    a5 = 6 * (y1 - y0) - 3 * a1
    y = a0 + a1 * u + a3 * u**3 + a4 * u**4 + a5 * u**5
    dy = (a1 + 3 * a3 * u * u + 4 * a4 * u**3 + 5 * a5 * u**4) / duration
    return float(y), float(dy)


def recovery_plan(plan):
    plan = deepcopy(plan)
    plan["scope"] = (
        "Continuous project-designed return; source times are reference markers, not reconstructed commands."
    )
    for event in plan["transitions"]:
        event["reference_time_s"] = event["time_s"]
        if event["name"] in {"entry_start", "landing_start"}:
            event["time_s"] = None
        if event["name"] == "stage_separation":
            event["phase"] = "booster_coast"
    r = plan["recovery"]
    R = plan["initial_state"]["structural_base_position_body_m"][0]
    x = r["landing_target_north_m"]
    y = -r["landing_target_east_m"]
    plan["landing_criteria"]["landing_target_ecef_m"] = [
        x,
        y,
        float(np.sqrt(R * R - x * x - y * y)),
    ]
    plan["landing_criteria"]["contact_height_m"] = r["contact_below_base_m"]
    return plan


class RecoveryController:
    def __init__(self, plan):
        self.plan = plan
        self.cfg = plan["recovery"]
        self.environment = Environment(**plan["environment"])
        self.aero = Aerodynamics.from_json(
            Path(__file__).parent / "data/aero_configurations.json"
        )
        self.target = np.array(plan["landing_criteria"]["landing_target_ecef_m"])
        self.phase = None
        self.time = 0.0
        self.started = 0.0
        self.next_update = 0.0
        self.ramp = (0.0, 1.0, 0.0, 0.0, 0.0)
        self.active = ()
        self.ending = False
        self.reason = ""
        self.phase_fuel = 0.0
        self.burns = []
        self.landing_end = None
        self.legs = False
        self.landing_ready = False
        self.predicted_miss = None
        self.prediction_time = -1e9
        self.steering = None
        self.control_period = self.cfg["guidance_period_s"]

    def command(self, t):
        value, derivative = smooth_command(t, *self.ramp)
        if value < 0 or value > 1:
            return float(np.clip(value, 0, 1)), 0.0
        return value, derivative

    def set_throttle(self, t, value, duration=None):
        if abs(value - self.ramp[-1]) < 1e-12:
            return
        current, derivative = self.command(t)
        self.ramp = (
            t,
            duration or self.control_period,
            current,
            derivative,
            value,
        )

    def flow_at(self, t):
        value, derivative = self.command(t)
        q = sum(
            steady_mass_flow(
                replace(self.vehicle.stages[0].engines[i], running=True, throttle=1.0)
            )
            for i in self.active
        )
        return np.array([q * value, 0.0]), np.array([q * derivative, 0.0])

    def provider(self, vehicle, phase, plan):
        from mission import default_provider

        if phase in {"powered_ascent", "engine_shutdown", "full_stack_coast"}:
            return default_provider(vehicle, phase, plan)
        self.vehicle = vehicle
        self.model = vehicle.mass_model()
        if phase != self.phase:
            self.phase = phase
            self.steering = None
            self.started = self.time
            self.ending = False
            self.reason = ""
            self.phase_fuel = vehicle.stages[0].propellant_mass
            self.active = (
                (0, 1, 5)
                if phase in {"boostback_burn", "entry_burn"}
                else (0,) if phase == "landing_burn" else ()
            )
            self.ramp = (self.time, self.cfg["ramp_s"], 0.0, 0.0, 0.0)
            for i, engine in enumerate(vehicle.stages[0].engines):
                engine.running = i in self.active
                engine.throttle = 1.0 if engine.running else 0.0
            if phase == "boostback_burn":
                self.set_throttle(
                    self.time, self.cfg["boostback_throttle"], self.cfg["ramp_s"]
                )
            if phase == "entry_burn":
                self.set_throttle(
                    self.time, self.cfg["entry_throttle"], self.cfg["ramp_s"]
                )
            if phase == "landing_burn":
                self.legs = True
                p, v = self.point_state(self.state, self.properties)
                up = p / np.linalg.norm(p)
                h = max(self.contact_altitude(self.state, self.properties), 1)
                self.landing_end = self.time + 2 * h / max(-np.dot(v, up) + 1.0, 1.0)
                self.set_throttle(self.time, 1.0, self.cfg["ramp_s"])
            self.next_update = self.time
        return self

    def point_state(self, state, properties, station=(0, 0, 0)):
        rot = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        arm = np.asarray(station) - properties.cg_station_body_m
        return rot @ (state.position_body_m + arm), rot @ (
            state.velocity_body_mps
            + np.cross(state.omega_body_radps, arm)
            - properties.cg_rate_body_mps
        )

    def contact_altitude(self, state, properties):
        c = self.cfg
        foot = np.array(
            [
                [-c["contact_below_base_m"], y, z]
                for y, z in [
                    (c["contact_radius_m"], 0),
                    (-c["contact_radius_m"], 0),
                    (0, c["contact_radius_m"]),
                    (0, -c["contact_radius_m"]),
                ]
            ]
        )
        points = (
            body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
            @ (state.position_body_m - properties.cg_station_body_m + foot).T
        ).T
        return float(
            np.min(np.linalg.norm(points, axis=1)) - self.environment.earth_radius_m
        )

    def geometry(self, state, props):
        p, v = self.point_state(state, props)
        up = p / np.linalg.norm(p)
        h = np.linalg.norm(p) - self.environment.earth_radius_m
        vr = np.dot(v, up)
        vh = v - vr * up
        g = np.linalg.norm(self.environment.gravity_ecef(p))
        flight = max((vr + np.sqrt(max(0, vr * vr + 2 * g * max(h, 0)))) / g, 1)
        error = p - self.target
        horizontal = error - np.dot(error, up) * up
        miss = horizontal + vh * flight
        return p, v, up, h, vr, vh, g, flight, horizontal, miss

    def direction(self, state, props, t, aero_force=None):
        if self.steering is not None and self.phase in {
            "entry_burn",
            "aerodynamic_descent",
            "landing_burn",
        }:
            return self.steering.copy()
        p, v, up, h, vr, vh, g, flight, error, miss = self.geometry(state, props)
        aero_accel = (
            np.zeros(3)
            if aero_force is None
            else body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
            @ aero_force
            / props.mass_kg
        )
        if self.phase in {"booster_coast", "boostback_burn"}:
            miss = miss if self.predicted_miss is None else self.predicted_miss
            acceleration = -miss / max(flight * 8.0, 1) + 0.8 * up
        elif self.phase == "landing_burn":
            rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
            goal = self.target + (
                props.cg_station_body_m[0] + self.cfg["contact_below_base_m"]
            ) * self.target / np.linalg.norm(self.target)
            pcg = rotation @ state.position_body_m
            vcg = rotation @ state.velocity_body_mps
            tgo = max(
                2
                * max(self.contact_altitude(state, props), 0.05)
                / max(-vr - self.cfg["landing_target_vertical_velocity_mps"], 1),
                0.3,
            )
            gravity = self.environment.gravity_ecef(pcg)
            acceleration = (
                6 * (goal - pcg) / tgo**2
                - 4 * vcg / tgo
                - 2 * self.cfg["landing_target_vertical_velocity_mps"] * up / tgo
                - gravity
                - aero_accel
            )
            vertical = max(np.dot(acceleration, up), 0.1)
            lateral = acceleration - np.dot(acceleration, up) * up
            max_lateral = vertical * np.tan(
                np.deg2rad(self.cfg["landing_tilt_limit_deg"])
            )
            if np.linalg.norm(lateral) > max_lateral:
                lateral *= max_lateral / np.linalg.norm(lateral)
            acceleration = vertical * up + lateral
        elif self.phase == "entry_burn":
            desired_horizontal = -error / max(flight + 25.0, 30.0)
            lateral = (desired_horizontal - vh) / 8.0
            acceleration = 30.0 * up + lateral - aero_accel
        elif h < 8000:
            acceleration = up
        else:
            acceleration = -v / max(np.linalg.norm(v), 1)
        return acceleration

    def predict_impact_miss(self, state, props):
        p, v = self.point_state(state, props)
        mass = props.mass_kg
        radius = self.environment.earth_radius_m
        entry_spent = 0.0
        entry_done = False
        model = self.aero.configurations["booster_return"]
        area = model.reference_area_m2
        side_area = sum(
            x.side_area_m2 * x.crossflow_drag_coefficient for x in model.side_segments
        )
        engine = self.vehicle.stages[0].engines[0]
        full_flow = steady_mass_flow(replace(engine, running=True, throttle=1.0))
        for _ in range(1400):
            up = p / np.linalg.norm(p)
            height = np.linalg.norm(p) - radius
            if height <= 0:
                break
            speed = max(np.linalg.norm(v), 1)
            vr = v @ up
            vh = v - vr * up
            air = self.environment.atmosphere(max(height, 0))
            q = 0.5 * air.density_kgm3 * speed**2
            nose = -v / speed if height > 8000 else up
            cosine = v @ nose / speed
            lateral = v - (v @ nose) * nose
            sine = np.linalg.norm(lateral) / speed
            coefficient = model.axial_drag(min(speed / air.sound_speed_mps, 10))
            force = -q * area * coefficient * cosine * abs(cosine) * nose
            if sine > 1e-9:
                force -= (
                    (
                        q
                        * area
                        * model.potential_normal_slope_per_rad
                        * sine
                        * abs(cosine)
                        + q * side_area * sine * sine
                    )
                    * lateral
                    / np.linalg.norm(lateral)
                )
            gravity = self.environment.gravity_ecef(p)
            thrust = np.zeros(3)
            if height < self.cfg["entry_altitude_m"] and vr < 0 and not entry_done:
                if (
                    speed <= self.cfg["entry_end_speed_mps"]
                    or entry_spent >= self.cfg["entry_budget_kg"]
                ):
                    entry_done = True
                else:
                    g = np.linalg.norm(gravity)
                    flight = max((vr + np.sqrt(vr * vr + 2 * g * height)) / g, 1)
                    error = p - self.target
                    error -= error @ up * up
                    desired_horizontal = -error / max(flight + 25, 30)
                    acceleration = (
                        30 * up + (desired_horizontal - vh) / 8 - force / mass
                    )
                    thrust = (
                        3
                        * (
                            self.cfg["entry_throttle"] * engine.spec.thrust_vac
                            - engine.spec.exit_area * air.pressure_pa
                        )
                        * acceleration
                        / np.linalg.norm(acceleration)
                    )
                    burned = 3 * self.cfg["entry_throttle"] * full_flow * 0.5
                    mass -= burned
                    entry_spent += burned
            acceleration = (
                gravity
                + self.environment.frame_acceleration_ecef(p, v)
                + (force + thrust) / mass
            )
            p = p + v * 0.5 + 0.125 * acceleration
            v = v + 0.5 * acceleration
        up = p / np.linalg.norm(p)
        error = p - self.target
        return error - (error @ up) * up

    def solve_force(self, state, props, desired_acceleration, throttle=None):
        p, v, up, h, vr, vh, g, flight, error, miss = self.geometry(state, props)
        first = np.array([1.0, 0, 0])
        first -= first @ up * up
        first /= np.linalg.norm(first)
        second = np.cross(up, first)
        model = self.aero.configurations[
            "booster_landing" if self.legs else "booster_return"
        ]
        air = self.environment.atmosphere(max(h, 0))
        speed = max(np.linalg.norm(v), 0.01)
        q = 0.5 * air.density_kgm3 * speed**2
        area = model.reference_area_m2
        side = sum(
            x.side_area_m2 * x.crossflow_drag_coefficient for x in model.side_segments
        )
        cd = model.axial_drag(min(speed / air.sound_speed_mps, 10))
        engine = self.vehicle.stages[0].engines[0]
        count = len(self.active)
        goal = props.mass_kg * (
            desired_acceleration
            - self.environment.gravity_ecef(p)
            - self.environment.frame_acceleration_ecef(p, v)
        )

        def force(x):
            nose = up + x[0] * first + x[1] * second
            nose /= np.linalg.norm(nose)
            cosine = v @ nose / speed
            lateral = v - (v @ nose) * nose
            sine = np.linalg.norm(lateral) / speed
            result = -q * area * cd * cosine * abs(cosine) * nose
            if sine > 1e-9:
                result -= (
                    (
                        q
                        * area
                        * model.potential_normal_slope_per_rad
                        * sine
                        * abs(cosine)
                        + q * side * sine * sine
                    )
                    * lateral
                    / np.linalg.norm(lateral)
                )
            value = x[2] if throttle is None else throttle
            if count:
                result += (
                    count
                    * (
                        value * engine.spec.thrust_vac
                        - engine.spec.exit_area * air.pressure_pa
                    )
                    * nose
                )
            return result, nose

        limit = np.tan(
            np.deg2rad(
                min(self.cfg["landing_tilt_limit_deg"], max(2.0, 0.02 * max(h, 0)))
                if self.phase == "landing_burn"
                else self.cfg["landing_tilt_limit_deg"]
            )
        )
        best = None
        center = np.zeros(2)
        step = limit / 4
        for iteration in range(3):
            values = range(-4, 5) if iteration == 0 else range(-2, 3)
            for i in values:
                for j in values:
                    angles = center + step * np.array([i, j])
                    if np.linalg.norm(angles) > limit + 1e-12:
                        continue
                    x = np.r_[angles, 0.0]
                    value, nose = force(x)
                    if throttle is None:
                        pressure_force = (
                            count * engine.spec.exit_area * air.pressure_pa * nose
                        )
                        passive = value + pressure_force
                        required = (goal - passive) @ up / max(nose @ up, 0.1)
                        x[2] = np.clip(
                            (required / count + engine.spec.exit_area * air.pressure_pa)
                            / engine.spec.thrust_vac,
                            engine.spec.throttle_min,
                            1,
                        )
                        value, nose = force(x)
                    error_force = value - goal
                    horizontal = error_force - (error_force @ up) * up
                    score = (
                        horizontal @ horizontal
                        + (100 * (error_force @ up) ** 2 if throttle is None else 0)
                        + 1e-4 * (angles @ angles)
                    )
                    if best is None or score < best[0]:
                        best = (score, x.copy(), nose.copy())
            center = best[1][:2]
            step *= 0.25
        return best[2], float(best[1][2] if throttle is None else throttle)

    def stop_burn(self, t, reason):
        if not self.ending:
            self.ending = True
            self.reason = reason
            self.set_throttle(t, 0.0, self.cfg["ramp_s"])

    def conditions(self, t, vehicle, state, ev, sequence):
        self.time = t
        self.state = state
        self.properties = ev.properties
        if sequence.phase not in {
            "booster_coast",
            "boostback_burn",
            "return_coast",
            "entry_burn",
            "aerodynamic_descent",
            "landing_burn",
        }:
            return {}
        self.vehicle = vehicle
        p, v, up, h, vr, vh, g, flight, error, miss = self.geometry(
            state, ev.properties
        )
        if (
            self.phase in {"booster_coast", "boostback_burn"}
            and t - self.prediction_time >= 1 - 1e-9
        ):
            self.predicted_miss = self.predict_impact_miss(state, ev.properties)
            self.prediction_time = t
        if (
            self.phase in {"booster_coast", "boostback_burn"}
            and self.predicted_miss is not None
        ):
            miss = self.predicted_miss
        rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        want = self.direction(state, ev.properties, t, ev.aero.force_body_n)
        want /= np.linalg.norm(want)
        angle = np.arccos(np.clip(rotation[:, 0] @ want, -1, 1))
        spent = self.phase_fuel - state.propellant_kg[0]
        guidance_due = t >= self.next_update - 1e-9
        if guidance_due:
            self.control_period = (
                min(
                    self.cfg["guidance_period_s"],
                    max(0.01, max(h, 0) / max(-vr, 1) / 10),
                )
                if self.phase == "landing_burn"
                else self.cfg["guidance_period_s"]
            )
            self.next_update = t + self.control_period
            if self.phase == "boostback_burn":
                if (
                    t - self.started >= self.cfg["boostback_min_s"]
                    and np.linalg.norm(miss) < self.cfg["boostback_miss_limit_m"]
                ):
                    self.stop_burn(t, "predicted ballistic miss tolerance")
                elif (
                    spent >= self.cfg["boostback_budget_kg"] - 500
                    or t - self.started >= self.cfg["boostback_max_s"]
                ):
                    self.stop_burn(t, "fuel or duration cap")
            elif self.phase == "entry_burn":
                tgo = max(2 * h / max(-vr, 1), 20)
                desired_horizontal = -6 * error / tgo**2 - 4 * vh / tgo
                self.steering, _ = self.solve_force(
                    state,
                    ev.properties,
                    20 * up + desired_horizontal,
                    self.cfg["entry_throttle"],
                )
                if np.linalg.norm(v) <= self.cfg["entry_end_speed_mps"]:
                    self.stop_burn(t, "entry speed target")
                elif spent >= self.cfg["entry_budget_kg"] - 500:
                    self.stop_burn(t, "entry budget cap")
            elif self.phase == "landing_burn" and not self.ending:
                rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
                goal = self.target + (
                    ev.properties.cg_station_body_m[0]
                    + self.cfg["contact_below_base_m"]
                ) * self.target / np.linalg.norm(self.target)
                pcg = rotation @ state.position_body_m
                vcg = rotation @ state.velocity_body_mps
                tgo = max(
                    2
                    * max(self.contact_altitude(state, ev.properties), 0.05)
                    / max(-vr - self.cfg["landing_target_vertical_velocity_mps"], 1),
                    0.3,
                )
                net = (
                    6 * (goal - pcg) / tgo**2
                    - 4 * vcg / tgo
                    - 2 * self.cfg["landing_target_vertical_velocity_mps"] * up / tgo
                )
                horizontal_time = max(
                    2.0, tgo * self.cfg["landing_horizontal_horizon_fraction"]
                )
                horizontal = -6 * error / horizontal_time**2 - 4 * vh / horizontal_time
                net = (net @ up) * up + horizontal
                self.steering, throttle = self.solve_force(state, ev.properties, net)
                self.set_throttle(t, throttle)
                if (
                    spent
                    >= min(
                        self.phase_fuel
                        - self.model.propellant_bounds(None)[0][0]
                        - self.cfg["reserve_kg"],
                        self.cfg["landing_budget_kg"]
                        + sum(
                            max(
                                0,
                                self.cfg[
                                    (
                                        "boostback_budget_kg"
                                        if b["phase"] == "boostback_burn"
                                        else "entry_budget_kg"
                                    )
                                ]
                                - b["consumed_kg"],
                            )
                            for b in self.burns
                        ),
                    )
                    - 200
                ):
                    self.stop_burn(t, "landing budget cap")
            if self.phase == "aerodynamic_descent":
                tgo = max(2 * h / max(-vr, 1), 5)
                desired_horizontal = -6 * error / tgo**2 - 4 * vh / tgo
                self.steering, _ = self.solve_force(
                    state, ev.properties, -g * up + desired_horizontal, 0.0
                )
        amplitude, _ = self.command(t)
        complete = self.ending and amplitude == 0
        if complete:
            for e in vehicle.stages[0].engines:
                e.running = False
                e.throttle = 0.0
            if not self.burns or self.burns[-1]["phase"] != self.phase:
                self.burns.append(
                    {
                        "phase": self.phase,
                        "start_s": self.started,
                        "end_s": t,
                        "consumed_kg": float(spent),
                        "reason": self.reason,
                        "ballistic_miss_m": float(np.linalg.norm(miss)),
                    }
                )
        engine = vehicle.stages[0].engines[0]
        max_thrust = (
            engine.spec.thrust_vac
            - engine.spec.exit_area * self.environment.atmosphere(max(h, 0)).pressure_pa
        )
        max_net = max(max_thrust / ev.properties.mass_kg - g, 0.1)
        landing_ready = self.landing_ready
        if guidance_due and self.phase == "aerodynamic_descent" and vr < 0:
            height = max(h, 0)
            speed = vr
            mass = ev.properties.mass_kg
            engine = vehicle.stages[0].engines[0]
            q = steady_mass_flow(replace(engine, running=True, throttle=1.0))
            curve = self.aero.configurations["booster_alone"].axial_drag
            area = self.aero.configurations["booster_alone"].reference_area_m2
            for k in range(1000):
                air = self.environment.atmosphere(max(height, 0))
                mach = min(abs(speed) / air.sound_speed_mps, 10)
                drag = 0.5 * air.density_kgm3 * speed * speed * area * curve(mach)
                amplitude, _ = smooth_command(k * 0.1, 0, self.cfg["ramp_s"], 0, 0, 1)
                thrust = (
                    (engine.spec.thrust_vac - engine.spec.exit_area * air.pressure_pa)
                    * amplitude
                    * np.cos(np.deg2rad(self.cfg["landing_tilt_limit_deg"]))
                )
                acceleration = (thrust + drag) / mass - g
                height += speed * 0.1 + 0.005 * acceleration
                speed += acceleration * 0.1
                mass -= q * amplitude * 0.1
                if speed >= -1 or height <= 0:
                    break
            landing_ready = height <= self.cfg["landing_margin_m"]
            self.landing_ready = landing_ready

        budgets = {
            "boostback_start": sum(
                self.cfg[x]
                for x in [
                    "boostback_budget_kg",
                    "entry_budget_kg",
                    "landing_budget_kg",
                    "reserve_kg",
                ]
            ),
            "entry_start": sum(
                self.cfg[x]
                for x in ["entry_budget_kg", "landing_budget_kg", "reserve_kg"]
            ),
            "landing_start": self.cfg["landing_budget_kg"] + self.cfg["reserve_kg"],
        }
        next_name = sequence.next_event.name if sequence.next_event else ""
        tilt = np.degrees(np.arccos(np.clip(rotation[:, 0] @ up, -1, 1)))
        metrics = {
            "ground_contact": self.contact_altitude(state, ev.properties) <= 1e-5,
            "legs_deployed": self.legs,
            "vertical_speed_mps": abs(float(vr)),
            "horizontal_speed_mps": float(np.linalg.norm(vh)),
            "tilt_deg": float(tilt),
            "angular_rate_degps": float(
                np.degrees(np.linalg.norm(state.omega_body_radps))
            ),
            "target_error_m": float(np.linalg.norm(error)),
        }
        return {
            "required_propellant_kg": budgets.get(next_name),
            "attitude_ready": angle < np.deg2rad(25),
            "propulsion_ready": self.cfg["assume_restart_ready"],
            "burn_complete": complete,
            "entry_start": vr < 0 and h <= self.cfg["entry_altitude_m"],
            "landing_start": landing_ready,
            "ground_contact": metrics["ground_contact"],
            "landing_metrics": metrics,
        }

    def evaluate(self, t, state):
        value, derivative = self.command(t)
        flow, flowdot = self.flow_at(t)
        props = self.model.mass_properties(
            None, state.propellant_kg, flow, consumption_rate_kgps2=flowdot
        )
        altitude = self.environment.altitude(state, props)
        air = self.environment.atmosphere(altitude)
        forces = []
        ports = []
        thrust = 0.0
        for i in self.active:
            engine = self.vehicle.stages[0].engines[i]
            minimum = engine.spec.throttle_min
            actual = max(value, minimum)
            factor = value / actual
            height = self.environment.altitude(
                state, props, [0, engine.position[1], engine.position[2]]
            )
            result = evaluate_steady_engine(
                replace(engine, running=True, throttle=actual),
                self.environment.atmosphere(height).pressure_pa,
                stage_index=0,
                stage_offset_m=0,
                nozzle_exit_axial_m=0,
            )
            forces.append(
                AppliedForce(
                    result.applied_force.force_body_n * factor,
                    result.applied_force.station_body_m,
                )
            )
            q = steady_mass_flow(replace(engine, running=True, throttle=1.0))
            ports.append(
                ExhaustPort(
                    0,
                    result.exhaust_port.station_body_m,
                    q * value,
                    q * derivative,
                    result.exhaust_port.gyration_body_m2,
                )
            )
            thrust += result.thrust_n * factor
        configuration = (
            "booster_landing"
            if self.legs
            else (
                "booster_return"
                if self.phase in {"return_coast", "entry_burn", "aerodynamic_descent"}
                else "booster_alone"
            )
        )
        aero = self.aero.evaluate(
            state, props, air, self.environment.wind_body(state), configuration
        )
        rotation = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        desired = self.direction(state, props, t, aero.force_body_n)
        desired /= np.linalg.norm(desired)
        error = rotation.T @ np.cross(rotation[:, 0], desired)
        desired_rate = np.clip(
            self.cfg["attitude_kp"] * error,
            -self.cfg["max_angular_rate_radps"],
            self.cfg["max_angular_rate_radps"],
        )
        alpha = self.cfg["attitude_kd"] * (desired_rate - state.omega_body_radps)
        length = np.linalg.norm(alpha)
        if length > self.cfg["max_angular_accel_radps2"]:
            alpha *= self.cfg["max_angular_accel_radps2"] / length
        moment = props.inertia_body_kgm2 @ alpha + np.cross(
            state.omega_body_radps, props.inertia_body_kgm2 @ state.omega_body_radps
        )
        moment -= aero.free_moment_body_nm
        fin = self.aero.configurations[configuration].grid_fin_control
        if fin is not None and aero.dynamic_pressure_pa > 1:
            basis = []
            cosine = np.cos(aero.angle_of_attack_rad)
            scale = (
                aero.dynamic_pressure_pa
                * fin.fin_area_m2
                * fin.lift_slope_per_rad
                * fin.maximum_deflection_rad
                * cosine
                * abs(cosine)
            )
            for azimuth in fin.azimuths_rad:
                force = scale * np.array([0, -np.sin(azimuth), np.cos(azimuth)])
                station = np.array(
                    [
                        fin.station_m,
                        fin.pivot_radius_m * np.cos(azimuth),
                        fin.pivot_radius_m * np.sin(azimuth),
                    ]
                )
                basis.append(np.cross(station - props.cg_station_body_m, force))
            fraction = np.clip(
                np.linalg.lstsq(np.array(basis).T, moment, rcond=None)[0], -1, 1
            )
            controlled = self.aero.evaluate(
                state,
                props,
                air,
                self.environment.wind_body(state),
                configuration,
                fraction * fin.maximum_deflection_rad,
            )
            moment -= controlled.free_moment_body_nm - aero.free_moment_body_nm
            aero = controlled
        forces.append(AppliedForce(aero.force_body_n, props.cg_station_body_m))
        limit = (
            self.cfg["max_control_moment_nm"]
            if self.active
            else self.cfg["coast_control_moment_nm"]
        )
        if np.linalg.norm(moment) > limit:
            moment *= limit / np.linalg.norm(moment)
        loads = assemble_loads(
            props,
            forces,
            gravity_body_mps2=self.environment.gravity_body(state),
            frame_acceleration_body_mps2=self.environment.frame_acceleration_body(
                state
            ),
            free_moment_body_nm=aero.free_moment_body_nm + moment,
            consumption_kgps=flow,
            exhaust_ports=ports,
        )
        return AscentEvaluation(
            loads,
            props,
            air,
            aero,
            thrust,
            altitude,
            self.environment.altitude(state, props, [0, 0, 0]),
        )

    def __call__(self, t, state):
        return self.evaluate(t, state).loads
