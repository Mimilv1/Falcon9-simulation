from dataclasses import dataclass, replace
import numpy as np
from propulsion import steady_mass_flow, evaluate_steady_engine, scale_shutdown_output
from load_interface import AppliedForce, assemble_loads
from Aerodynamics import Aerodynamics
from frames import body_to_ecef_matrix


@dataclass(frozen=True)
class AscentEvaluation:
    loads: object
    properties: object
    atmosphere: object
    aero: object
    thrust_n: float
    cg_altitude_m: float
    base_altitude_m: float


class SteadyPoweredAscent:

    def __init__(
        self,
        vehicle,
        environment,
        aerodynamics,
        *,
        nozzle_exit_axial_m=0.0,
        max_angle_of_attack_rad=None,
        command_envelope=None,
        aero_configuration="full_stack"
    ):
        self.vehicle = vehicle
        self.model = vehicle.mass_model()
        self.environment = environment
        self.aerodynamics = aerodynamics
        if isinstance(aerodynamics, Aerodynamics):
            self.aero_configuration = aero_configuration
            if aero_configuration not in aerodynamics.configurations:
                raise ValueError(
                    "Unknown aero configuration for this provider: "
                    + str(aero_configuration)
                )
        else:
            self.aero_configuration = None
        self.nozzle_exit_axial_m = float(nozzle_exit_axial_m)
        self.command_envelope = command_envelope
        self.max_angle_of_attack_rad = max_angle_of_attack_rad
        if max_angle_of_attack_rad is not None and (
            not np.isfinite(max_angle_of_attack_rad) or max_angle_of_attack_rad <= 0
        ):
            raise ValueError(
                "A configured angle-of-attack check limit must be positive."
            )
        self.stage1 = next((s for s in vehicle.stages if s.name == "stage1"))
        self._commands = tuple(
            (
                (e.running, e.throttle, e.gimbal_pitch, e.gimbal_yaw)
                for e in self.stage1.engines
            )
        )
        self.consumption_kgps = np.array(
            [
                sum(
                    (
                        steady_mass_flow(e, vehicle.parameters.G0)
                        for e in self.stage1.engines
                    )
                ),
                0.0,
            ]
        )

    def flow_at(self, time_s):
        amplitude, derivative = (
            (1.0, 0.0)
            if self.command_envelope is None
            else self.command_envelope(time_s)
        )
        return (self.consumption_kgps * amplitude, self.consumption_kgps * derivative)

    def evaluate(self, time_s, state):
        if self._commands != tuple(
            (
                (e.running, e.throttle, e.gimbal_pitch, e.gimbal_yaw)
                for e in self.stage1.engines
            )
        ):
            raise ValueError(
                "Steady commands changed: construct a new properly bounded phase provider."
            )
        consumption, consumption_derivative = self.flow_at(time_s)
        properties = self.model.mass_properties(
            None,
            state.propellant_kg,
            consumption,
            consumption_rate_kgps2=consumption_derivative,
        )
        altitude = self.environment.altitude(state, properties)
        atmosphere = self.environment.atmosphere(altitude)
        forces = []
        ports = []
        thrust = 0.0
        for engine in self.stage1.engines:
            exit_station = np.array(
                [
                    self.stage1.offset + self.nozzle_exit_axial_m,
                    engine.position[1],
                    engine.position[2],
                ]
            )
            exit_height = self.environment.altitude(state, properties, exit_station)
            ambient = self.environment.atmosphere(exit_height).pressure_pa
            result = evaluate_steady_engine(
                engine,
                ambient,
                stage_index=0,
                stage_offset_m=self.stage1.offset,
                nozzle_exit_axial_m=self.nozzle_exit_axial_m,
                g0=self.vehicle.parameters.G0,
            )
            if self.command_envelope is not None:
                result = scale_shutdown_output(result, *self.command_envelope(time_s))
            forces.append(result.applied_force)
            ports.append(result.exhaust_port)
            thrust += result.thrust_n
        wind = self.environment.wind_body(state)
        if self.aero_configuration is not None:
            aero = self.aerodynamics.evaluate(
                state, properties, atmosphere, wind, self.aero_configuration
            )
        else:
            aero = self.aerodynamics.evaluate(state, properties, atmosphere, wind)
        if (
            self.max_angle_of_attack_rad is not None
            and aero.relative_speed_mps > 1e-06
            and (aero.angle_of_attack_rad > self.max_angle_of_attack_rad)
        ):
            raise ValueError(
                "Attitude check exceeded its declared angle-of-attack limit."
            )
        forces.append(AppliedForce(aero.force_body_n, properties.cg_station_body_m))
        aero_moment = aero.free_moment_body_nm
        if aero_moment is None:
            aero_moment = np.zeros(3)
        loads = assemble_loads(
            properties,
            forces,
            gravity_body_mps2=self.environment.gravity_body(state),
            free_moment_body_nm=aero_moment,
            consumption_kgps=consumption,
            exhaust_ports=ports,
            frame_acceleration_body_mps2=self.environment.frame_acceleration_body(
                state
            ),
        )
        return AscentEvaluation(
            loads,
            properties,
            atmosphere,
            aero,
            thrust,
            altitude,
            self.environment.altitude(state, properties, [0, 0, 0]),
        )

    def __call__(self, time_s, state):
        return self.evaluate(time_s, state).loads


@dataclass(frozen=True)
class PitchProgram:
    start_s: float = 8.0
    tau_s: float = 60.0
    tilt_max_rad: float = np.deg2rad(80.0)

    def tilt(self, time_s):
        if time_s <= self.start_s:
            return 0.0
        x = (time_s - self.start_s) / self.tau_s
        return self.tilt_max_rad * x * x / (1.0 + x * x)


class ControlledAscent:
    EAST_ECEF = np.array([0.0, -1.0, 0.0])

    def __init__(self, ascent, program, kp=2.0, kd=2.5):
        self.ascent, self.program, self.kp, self.kd = (ascent, program, kp, kd)
        self.last = None

    def flow_at(self, time_s):
        return self.ascent.flow_at(time_s)

    @property
    def environment(self):
        return self.ascent.environment

    def evaluate(self, time_s, state):
        ev = self.ascent.evaluate(time_s, state)
        self.last = ev
        rot = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        up = rot @ state.position_body_m
        up = up / np.linalg.norm(up)
        east = self.EAST_ECEF - self.EAST_ECEF @ up * up
        east /= np.linalg.norm(east)
        theta = self.program.tilt(time_s)
        wanted = np.cos(theta) * up + np.sin(theta) * east
        error_body = rot.T @ np.cross(rot[:, 0], wanted)
        alpha = self.kp * error_body - self.kd * state.omega_body_radps
        torque = ev.properties.inertia_body_kgm2 @ alpha
        loads = replace(ev.loads, moment_body_nm=ev.loads.moment_body_nm + torque)
        return replace(ev, loads=loads)

    def __call__(self, time_s, state):
        return self.evaluate(time_s, state).loads
