"""First B propulsion model: steady throttle, fixed nozzles, axial thrust."""
from dataclasses import dataclass
import numpy as np
from load_interface import AppliedForce
from mass_flow import ExhaustPort,uniform_disk_gyration


@dataclass(frozen=True)
class EngineOutput:
    applied_force: AppliedForce
    exhaust_port: ExhaustPort
    thrust_n: float
    mass_flow_kgps: float


def steady_mass_flow(engine,g0=9.80665):
    if not engine.running or engine.throttle == 0:
        return 0.0
    throttle=float(engine.throttle)
    minimum=engine.spec.throttle_min or 0.0
    if not np.isfinite(throttle) or not minimum <= throttle <= 1:
        raise ValueError('Steady engine command outside adopted throttle limits.')
    if engine.spec.isp_vac <= 0 or engine.spec.thrust_vac <= 0 or g0 <= 0:
        raise ValueError('Positive vacuum thrust, Isp and standard gravity are required.')
    return throttle*engine.spec.thrust_vac/(engine.spec.isp_vac*g0)


def evaluate_steady_engine(engine,ambient_pressure_pa, *, stage_index,stage_offset_m,
                           nozzle_exit_axial_m,g0=9.80665):
    """T=throttle*Tvac-Ae*p, q=throttle*Tvac/(Isp_vac*g0).

    Constant vacuum Isp, fixed exit area, uniform exit mass-flux disk and a
    coplanar exit station are declared approximations. No transient/gimbal law.
    """
    if not np.isfinite(ambient_pressure_pa) or ambient_pressure_pa < 0:
        raise ValueError('Ambient pressure must be finite and nonnegative.')
    if engine.gimbal_pitch != 0 or engine.gimbal_yaw != 0:
        raise ValueError('This first propulsion model supports zero fixed gimbal only.')
    area=engine.spec.exit_area
    if area is None or not np.isfinite(area) or area <= 0:
        raise ValueError('A positive adopted nozzle exit area is required.')
    flow=steady_mass_flow(engine,g0)
    thrust=engine.throttle*engine.spec.thrust_vac-area*ambient_pressure_pa if flow > 0 else 0.0
    if thrust < 0:
        raise ValueError('Engine pressure/throttle combination is outside this model.')
    station=np.asarray(engine.position,dtype=float)+np.array([stage_offset_m,0,0])
    exit_station=np.array([stage_offset_m+nozzle_exit_axial_m,engine.position[1],engine.position[2]])
    port=ExhaustPort(stage_index,exit_station,flow,0.0,uniform_disk_gyration(np.sqrt(area/np.pi)))
    return EngineOutput(AppliedForce([thrust,0,0],station),port,float(thrust),float(flow))


def scale_shutdown_output(steady_output,amplitude,amplitude_rate_per_s):
    """Verification envelope scales BOTH full thrust and flow to zero.

    This is not steady operation below the minimum throttle, or a measured
    shutdown law. Scaling the entire ambient-corrected thrust prevents a
    residual -Ae*p force when an engine is off. Derivative propagates into CG.
    """
    if not np.all(np.isfinite([amplitude,amplitude_rate_per_s])) or not 0<=amplitude<=1:
        raise ValueError('Invalid shutdown amplitude or derivative.')
    port=steady_output.exhaust_port
    return EngineOutput(
        AppliedForce(amplitude*steady_output.applied_force.force_body_n,
                     steady_output.applied_force.station_body_m),
        ExhaustPort(port.stage_index,port.station_body_m,amplitude*port.mass_flow_kgps,
                    amplitude_rate_per_s*port.mass_flow_kgps,port.gyration_body_m2),
        amplitude*steady_output.thrust_n,amplitude*steady_output.mass_flow_kgps)
