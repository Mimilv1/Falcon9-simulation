"""Team integration entry point: Stage 1 powered ascent to MECO.

Pipeline: VehicleData -> ascent loads (thrust + drag + gravity + Earth rotation)
-> six-DOF RK4 -> IMU / GNSS / LIDAR truth + error models -> navigation Kalman
filter -> 3D matplotlib comparison of the true and the estimated trajectory.

Flight events are NOT evidence-based yet; the explicit assumptions are below and
every one of them can be overridden from the command line (see --help):

* all 9 Merlin engines at full throttle, fixed nozzles (no gimbal law in the
  propulsion model);
* attitude: an ideal pitch program (vertical, then smooth tilt toward the East)
  tracked by a PD law that applies a PURE COUPLE about the CG. It stands in for
  the missing thrust-vector control; no force is added and no authority limit
  is enforced;
* MECO: commanded at --t-meco with the smooth quintic shutdown envelope of
  engine_schedule.SmoothShutdown (verification envelope, not a measured law);
* no landing reserve is modeled: MECO time is not the real CRS-10 one.
"""
import argparse
import json
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from Aerodynamics import DragOnlyAerodynamics, MachDragCurve
from ascent import SteadyPoweredAscent
from dynamics import rk4_step
from engine_schedule import SmoothShutdown
from Environement import Environment
from frames import body_to_ecef_matrix
from kalman import NavigationKF
from mass_flow import PropellantLimitError
from sensors import (GNSS, IMU, GNSSSpec, IMUSpec, Lidar, LidarSpec,
                     quat_to_matrix, truth_inputs)
from VehicleData import build_falcon9

DRAG_CURVE_PATH = Path(__file__).parent / "data" / "drag_curve_pezent.json"
FALLBACK_CD = 0.3  # used only when the sourced drag table is missing

# Fixed structural stations (Body +x from the stage-1 base). PLACEHOLDERS.
IMU_STATION_M = [38.0, 0.0, 0.0]
GNSS_ANTENNA_STATION_M = [41.0, 0.0, 0.0]
LIDAR_STATION_M = [0.5, 0.0, 0.0]


# --------------------------------------------------------------------------- #
# Aerodynamics
# --------------------------------------------------------------------------- #
class ClampedDrag:
    """Mach readout clamped to the table range (the raw curve rejects outside)."""

    def __init__(self, curve):
        self.curve = curve

    def __call__(self, mach):
        return self.curve(float(np.clip(mach, self.curve.mach[0], self.curve.mach[-1])))


def build_aerodynamics(vehicle):
    if DRAG_CURVE_PATH.exists():
        area = json.loads(DRAG_CURVE_PATH.read_text(encoding="utf-8"))["reference_area_m2"]
        curve = MachDragCurve.from_json(DRAG_CURVE_PATH, area)
        note = (f"Cd(M) from {DRAG_CURVE_PATH.name} (Mach {curve.mach[0]:.2f}-{curve.mach[-1]:.2f}, "
                "clamped outside)")
        return DragOnlyAerodynamics(area, ClampedDrag(curve)), note
    area = vehicle.parameters.body_area
    return (DragOnlyAerodynamics(area, FALLBACK_CD),
            f"drag table not found: constant Cd = {FALLBACK_CD}")


# --------------------------------------------------------------------------- #
# Attitude program (stands in for thrust-vector control)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PitchProgram:
    """Tilt from local vertical toward the East: 0 until start_s, then smooth."""
    start_s: float = 8.0
    tau_s: float = 60.0
    tilt_max_rad: float = np.deg2rad(80.0)

    def tilt(self, time_s):
        if time_s <= self.start_s:
            return 0.0
        x = (time_s - self.start_s) / self.tau_s
        return self.tilt_max_rad * x * x / (1.0 + x * x)


class ControlledAscent:
    """Wraps the B load provider and adds an ideal attitude-control couple."""

    EAST_ECEF = np.array([0.0, -1.0, 0.0])  # see Environement frame convention

    def __init__(self, ascent, program, kp=2.0, kd=2.5):
        self.ascent, self.program, self.kp, self.kd = ascent, program, kp, kd
        self.last = None  # last AscentEvaluation (for logging)

    def __call__(self, time_s, state):
        ev = self.ascent.evaluate(time_s, state)
        self.last = ev
        rot = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz)
        up = rot @ state.position_body_m
        up = up / np.linalg.norm(up)
        east = self.EAST_ECEF - (self.EAST_ECEF @ up) * up
        east /= np.linalg.norm(east)
        theta = self.program.tilt(time_s)
        wanted = np.cos(theta) * up + np.sin(theta) * east
        error_body = rot.T @ np.cross(rot[:, 0], wanted)  # nose -> wanted axis
        alpha = self.kp * error_body - self.kd * state.omega_body_radps
        torque = ev.properties.inertia_body_kgm2 @ alpha
        return replace(ev.loads, moment_body_nm=ev.loads.moment_body_nm + torque)


# --------------------------------------------------------------------------- #
# Simulation
# --------------------------------------------------------------------------- #
def structural_truth(state, props, station):
    """Truth position/velocity (ECEF) of a fixed structural point."""
    rot = quat_to_matrix(state.quaternion_body_to_ecef_wxyz)
    arm = np.asarray(station) - props.cg_station_body_m
    return (rot @ (state.position_body_m + arm),
            rot @ (state.velocity_body_mps + np.cross(state.omega_body_radps, arm)
                   - props.cg_rate_body_mps))


def simulate(args):
    dt = args.dt
    if not 0 < dt <= 0.1:
        raise ValueError("--dt must be in (0, 0.1] s (sensors run at up to 10 Hz).")
    for name in ("t_meco", "shutdown"):
        n = getattr(args, name) / dt
        if abs(n - round(n)) > 1e-9:
            raise ValueError(f"--{name.replace('_', '-')} must be a multiple of --dt.")

    vehicle = build_falcon9()
    for engine in vehicle.stages[0].engines:  # all stage-1 engines, full throttle
        engine.running, engine.throttle = True, 1.0

    env = Environment(launch_latitude_deg=args.latitude)
    aero, drag_note = build_aerodynamics(vehicle)
    shutdown = SmoothShutdown(args.t_meco, args.shutdown)
    ascent = SteadyPoweredAscent(vehicle, env, aero, command_envelope=shutdown)
    program = PitchProgram(args.pitch_start, args.pitch_tau, np.deg2rad(args.tilt_max))
    provider = ControlledAscent(ascent, program)
    model = ascent.model

    # Inventory check: MECO must come before the limiting tank is empty.
    flow = float(ascent.consumption_kgps[0])
    lower = float(model.propellant_bounds(None)[0][0])
    prop0 = float(vehicle.stages[0].propellant_mass)
    t_limit = (prop0 - lower) / flow - 0.5 * args.shutdown
    if args.t_meco > t_limit - 0.5:
        raise ValueError(f"--t-meco too late: the limiting tank empties; use < {t_limit - 0.5:.1f} s.")

    state = vehicle.to_state()
    t_end = args.t_meco + args.shutdown
    n_steps = int(round(t_end / dt))

    rate_hz = 1.0 / dt
    imu_spec = IMUSpec(IMU_STATION_M, rate_hz=rate_hz, accel_bias_mps2=[0.05, -0.03, 0.02],
                       accel_scale_factor=1e-3, accel_bias_rw=1e-4,
                       gyro_bias_radps=[0.0, 0.0, 1e-4])  # PLACEHOLDER errors
    gnss_spec = GNSSSpec(GNSS_ANTENNA_STATION_M)
    lidar_spec = LidarSpec(LIDAR_STATION_M, planet_radius_m=env.earth_radius_m)
    imu = IMU(imu_spec, seed=args.seed)
    gnss = GNSS(gnss_spec, seed=args.seed + 1)
    lidar = Lidar(lidar_spec, seed=args.seed + 2)
    gnss_arm = np.asarray(GNSS_ANTENNA_STATION_M) - np.asarray(IMU_STATION_M)
    lidar_arm = np.asarray(LIDAR_STATION_M) - np.asarray(IMU_STATION_M)

    t = 0.0
    loads, props = truth_inputs(t, state, model, None, provider)
    p0, v0 = structural_truth(state, props, IMU_STATION_M)
    kf = NavigationKF.from_imu_spec(  # initial error: pad survey offset, at rest
        imu_spec, p0 + [1.0, -1.0, 0.5], v0, [2, 2, 2, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1],
        environment=env)

    log = {k: [] for k in ("t", "truth_p", "truth_v", "cg", "alt", "speed", "mach",
                           "q", "thrust", "mass", "tilt", "gnss")}
    memo = {"t": None, "w": np.zeros(3)}
    nis, lidar_updates = [], 0

    def process(t, state, loads, props):
        nonlocal lidar_updates
        rot = quat_to_matrix(state.quaternion_body_to_ecef_wxyz)  # attitude truth (V1)
        s = imu.update(t, state, loads, props)
        if s is not None:
            if memo["t"] is not None:
                kf.predict(s.accel_meas_mps2, rot, s.t - memo["t"])
            memo["t"], memo["w"] = s.t, s.gyro_meas_radps
        g = gnss.update(t, state, props)
        if g is not None and g.valid:
            kf.update_gnss(g.pos_meas_m, g.vel_meas_mps, rot, memo["w"], gnss_arm,
                           g.pos_sigma_m, g.vel_sigma_mps, gnss_spec.latency_s)
            nis.append(kf.last_nis)
            log["gnss"].append(g.pos_meas_m + g.vel_meas_mps * gnss_spec.latency_s - rot @ gnss_arm)
        lid = lidar.update(t, state, props)
        if lid is not None and lid.valid:
            if kf.update_lidar(lid.range_meas_m, rot, lidar_arm, lidar_spec.boresight_body,
                               lidar_spec.planet_center_ecef_m, lidar_spec.planet_radius_m,
                               lid.sigma_m) is not None:
                lidar_updates += 1
        if s is not None:
            kf.record(t)
            tp, tv = structural_truth(state, props, IMU_STATION_M)
            ev = provider.last
            cg_ecef = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz) @ state.position_body_m
            log["t"].append(t)
            log["truth_p"].append(tp)
            log["truth_v"].append(tv)
            log["cg"].append(cg_ecef)
            log["alt"].append(ev.cg_altitude_m)
            log["speed"].append(float(np.linalg.norm(state.velocity_body_mps)))
            log["mach"].append(ev.aero.mach)
            log["q"].append(ev.aero.dynamic_pressure_pa)
            log["thrust"].append(ev.thrust_n)
            log["mass"].append(props.mass_kg)
            nose = rot[:, 0]
            up = cg_ecef / np.linalg.norm(cg_ecef)
            log["tilt"].append(float(np.degrees(np.arccos(np.clip(nose @ up, -1, 1)))))

    process(t, state, loads, props)
    reason = "MECO (shutdown complete)"
    for k in range(1, n_steps + 1):
        try:
            state = rk4_step((k - 1) * dt, state, dt, model, None, provider)
        except PropellantLimitError:
            reason = "propellant limit reached before the commanded MECO"
            break
        t = k * dt
        loads, props = truth_inputs(t, state, model, None, provider)
        process(t, state, loads, props)
        if k % int(round(10.0 / dt)) == 0:
            print(f"  t = {t:6.1f} s | alt {log['alt'][-1]/1e3:6.1f} km | "
                  f"v {log['speed'][-1]:6.0f} m/s | Mach {log['mach'][-1]:4.1f} | "
                  f"tilt {log['tilt'][-1]:4.1f} deg")

    ts, X, SIG = kf.to_arrays()
    out = {key: np.asarray(val) for key, val in log.items()}
    out.update(est_p=X[:, 0:3], est_v=X[:, 3:6], est_sigma=SIG, nis=np.asarray(nis),
               lidar_updates=lidar_updates, reason=reason, drag_note=drag_note,
               earth_radius_m=env.earth_radius_m,
               propellant_left_kg=float(state.propellant_kg[0]))
    return out


# --------------------------------------------------------------------------- #
# Report and plots
# --------------------------------------------------------------------------- #
def to_local(p_ecef, radius_m):
    """ECEF -> launch tangent frame (East, North, Up) [m]; East = -y, North = +x, Up = z - R."""
    p = np.atleast_2d(p_ecef)
    return np.column_stack([-p[:, 1], p[:, 0], p[:, 2] - radius_m])


def report(r):
    err_p = np.linalg.norm(r["est_p"] - r["truth_p"], axis=1)
    err_v = np.linalg.norm(r["est_v"] - r["truth_v"], axis=1)
    local = to_local(r["truth_p"], r["earth_radius_m"])
    print(f"\nStop: {r['reason']}")
    print(f"MECO at t = {r['t'][-1]:.1f} s | altitude {r['alt'][-1]/1e3:.1f} km | "
          f"downrange {np.hypot(local[-1, 0], local[-1, 1])/1e3:.1f} km | "
          f"speed {r['speed'][-1]:.0f} m/s | Mach {r['mach'][-1]:.1f}")
    print(f"Mass at MECO {r['mass'][-1]:,.0f} kg | stage-1 propellant left "
          f"{r['propellant_left_kg']:,.0f} kg | max q {r['q'].max()/1e3:.1f} kPa "
          f"at t = {r['t'][r['q'].argmax()]:.0f} s")
    print(f"Navigation error: final position {err_p[-1]:.2f} m / velocity {err_v[-1]:.3f} m/s | "
          f"position rms {np.sqrt(np.mean(err_p**2)):.2f} m")
    print(f"GNSS mean NIS {r['nis'].mean():.2f} (6 dof -> ~6) | LIDAR updates {r['lidar_updates']}")
    print(f"Drag model: {r['drag_note']}")
    return err_p, err_v


def plot(r, err_p, err_v, equal_axes, save, show):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radius = r["earth_radius_m"]
    truth = to_local(r["truth_p"], radius) / 1e3
    est = to_local(r["est_p"], radius) / 1e3
    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    if len(r["gnss"]):  # drawn first so the two trajectories stay visible on top
        gnss = to_local(r["gnss"], radius)[::5] / 1e3
        ax.scatter(*gnss.T, s=2, color="0.65", alpha=0.35, label="Raw GNSS measurements")
    ax.plot(*truth.T, color="tab:blue", lw=2.2, label="True trajectory (IMU station)")
    ax.plot(*est.T, color="tab:red", lw=1.2, ls="--", label="Estimated trajectory (Kalman IMU+GNSS+LIDAR)")
    ax.scatter(*truth[0], marker="^", s=80, color="green", label="Launch")
    ax.scatter(*truth[-1], marker="*", s=160, color="black", label=f"MECO (t = {r['t'][-1]:.0f} s)")
    ax.set_xlabel("East [km]")
    ax.set_ylabel("North [km]")
    ax.set_zlabel("Altitude (tangent plane) [km]")
    ax.set_title("Stage 1 ascent to MECO")
    allp = np.vstack([truth, est])
    if equal_axes:
        center = (allp.max(axis=0) + allp.min(axis=0)) / 2
        half = max(np.ptp(allp, axis=0).max() / 2, 1e-3)
        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(max(center[2] - half, -0.5), center[2] + half)
        ax.set_box_aspect((1, 1, 1))
    else:
        mid_n = (allp[:, 1].max() + allp[:, 1].min()) / 2
        half_n = max(np.ptp(allp[:, 1]) / 2, 5.0)  # keep >= 10 km of North so ticks stay readable
        ax.set_ylim(mid_n - half_n, mid_n + half_n)
        spans = np.array([np.ptp(allp[:, 0]), 2 * half_n, np.ptp(allp[:, 2])])
        ax.set_box_aspect(np.maximum(spans, 1e-3) ** 0.5)
    ax.legend(loc="upper left", fontsize=8)

    ax2 = fig.add_subplot(2, 2, 2)
    sigma3 = 3 * np.linalg.norm(r["est_sigma"][:, 0:3], axis=1)
    ax2.plot(r["t"], err_p, color="tab:red", label="Position error (norm)")
    ax2.plot(r["t"], sigma3, color="0.5", ls=":", label="Filter 3σ (norm)")
    ax2.set_ylabel("Position [m]")
    ax2.set_yscale("log")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8)
    ax2.set_title("Estimation error")

    ax3 = fig.add_subplot(2, 2, 4, sharex=ax2)
    ax3.plot(r["t"], err_v, color="tab:orange", label="Velocity error (norm)")
    ax3.set_ylabel("Velocity [m/s]")
    ax3.set_xlabel("Time [s]")
    ax3.grid(alpha=0.3)
    ax3.legend(fontsize=8)

    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=150)
        print(f"Figure saved: {save}")
    if show:
        plt.show()


def parse_args():
    p = argparse.ArgumentParser(description="Stage 1 ascent to MECO with sensors and 3D plot.")
    p.add_argument("--t-meco", type=float, default=150.0, help="MECO command time [s] (default 150)")
    p.add_argument("--shutdown", type=float, default=1.0, help="smooth shutdown duration [s]")
    p.add_argument("--dt", type=float, default=0.01, help="integration/IMU step [s] (default 0.01)")
    p.add_argument("--pitch-start", type=float, default=8.0, help="start of the pitch-over [s]")
    p.add_argument("--pitch-tau", type=float, default=60.0, help="pitch-over time constant [s]")
    p.add_argument("--tilt-max", type=float, default=80.0, help="asymptotic tilt from vertical [deg]")
    p.add_argument("--latitude", type=float, default=28.6, help="launch latitude [deg]")
    p.add_argument("--seed", type=int, default=1, help="sensor noise seed")
    p.add_argument("--equal-axes", action="store_true", help="same scale on the 3 axes of the 3D plot")
    p.add_argument("--save", type=str, default=None, help="save the figure to this path")
    p.add_argument("--no-show", action="store_true", help="do not open the matplotlib window")
    return p.parse_args()


def main():
    args = parse_args()
    vehicle = build_falcon9()
    print(f"{vehicle.parameters.mission} model loaded: {vehicle.getmass():,.0f} kg. "
          f"Simulating to MECO (t = {args.t_meco:.0f} s)...")
    result = simulate(args)
    err_p, err_v = report(result)
    plot(result, err_p, err_v, args.equal_axes, args.save, show=not args.no_show)


if __name__ == "__main__":
    main()