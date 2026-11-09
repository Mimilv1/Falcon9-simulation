# CRS-10 mission configuration

The simulator follows the full vehicle until first/second-stage separation,
then follows the first-stage booster. Project 2 requires simulation and
navigation for selected representative phases; the full-return controller
is an additional demonstration.

## Flight phases

The full case preserves the propagated state through every transition.

| Phase / event | Transition condition | Nominal simulated time [s] |
| --- | --- | --- |
| Powered ascent | Liftoff initial state | 0 |
| Engine shutdown | Start the 1 s shutdown ramp | 140 |
| Full-stack coast | Zero thrust and propellant flow | 141 |
| Stage separation | 144 s, zero flow, engines off | 144 |
| Boostback ignition | 161 s, separated, fuel and attitude checks | 161 |
| Boostback end | Predicted miss goal, fuel cap or duration cap; smooth shutdown | 208.810 |
| Return coast | Continue first-stage propagation | 208.810–370.140 |
| Entry burn | Descending through 45 km, fuel and attitude checks | 370.140 |
| Entry burn end | Speed reaches 500 m/s or fuel cap; smooth shutdown | 392.950 |
| Aerodynamic descent | Continue to predicted landing-ignition condition | 392.950–426.160 |
| Landing burn | Predicted stopping clearance plus 1.5 km allowance | 426.160 |
| Touchdown | Lowest deployed leg tip reaches the spherical ground | 465.083 |

The listed times after boostback ignition are simulated outcomes, not observed
CRS-10 event times. The nominal boostback ends on its planning cap. Per-burn
caps are checked periodically and followed by a smooth shutdown. The published
preflight times adopted for MECO, separation and boostback ignition are
approximate references, used as A-class simulation commands.

## Configuration

Values are in `data/crs10_mission.json`; the full case builds its target and
state-triggered return events in `recovery.py`.

| Setting | Value / meaning |
| --- | --- |
| Body convention | +x toward nose; structural base is the station datum |
| External frame | Earth-centred, launch-aligned; initial base [0,0,6371000] m |
| Initial base velocity and angular rate | Zero relative to rotating Earth |
| Initial quaternion | Body-to-ECEF Hamilton [√0.5,0,−√0.5,0] |
| Propellant planning allocations | Boostback 38 t, entry 16 t, landing 10 t, reserve 2 t |
| Return target | Assumed offset 14 km south / 6 km east of launch, on the sphere |
| Leg geometry | Four assumed tips, 1 m below base and 3.5 m radial offset |
| Maximum integration step | 0.01 s |
| Random seed | 1 |

Unused boostback/entry allocation can transfer to landing. The reserve and
limiting-tank residual are protected. Restart health and propellant settling
are assumed; the 1 s flow ramps are numerical commands. Separation uses a
zero-impulse map with structural-point continuity and a conservation audit.

## Sensors and navigation

| Sensor | Sampling rate | Fixed structural station, Body [m] |
| --- | --- | --- |
| IMU | 100 Hz | [38,0,0] |
| GNSS | 10 Hz | [41,0,0] |
| LiDAR | 10 Hz | [0.5,0,0] |

The same sensor stations and filter instance continue through separation.
Sampling occurs at accepted integration times. Noise and bias settings are
recorded in each run's `summary.json` and remain provisional sensor assumptions.

The nine-state filter estimates IMU-station position, velocity and accelerometer
bias. It receives truth attitude. It does not currently estimate attitude or
gyro bias. Return guidance uses the true vehicle state rather than the filter
output. Navigation error and disagreement with public flight data are separate
comparisons.

## Verified nominal result

Computational source: `54d6312`; default 100 Hz sensor integration, seed 1.
The run starts at liftoff and ends at actual leg contact.

| Contact metric | Result | Acceptance limit |
| --- | --- | --- |
| Vertical speed [m/s] | 0.497 | ≤2 |
| Horizontal speed [m/s] | 0.902 | ≤1 |
| Target error [m] | 3.306 | ≤10 |
| Tilt [deg] | 0.535 | ≤5 |
| Angular rate [deg/s] | 0.013 | ≤1 |

Result: `landing_success`, contact at 465.083 s, propellant remaining 8.931 t.
Saved-sample navigation position RMS: 0.918 m with truth attitude.
These are nominal project-model results, not actual Falcon 9 landing accuracy.

The compact verification record is in `data/ascent_model_audit.json`.
Generated trajectories remain local under `results/`.

## Model scope and evidence

- Vehicle masses, geometry, CG and inertia: [parameter evidence](parameter_evidence.md).
- Thrust, generic full-stack drag, estimated booster coefficients and Mach limits:
  [aerodynamic evidence](aero_model.md) and the parameter record.
- Moving-CG dynamics, exhaust assumptions and ideal attitude moments:
  [dynamics](powered_dynamics.md) and [A/B/C interface](core_api.md).
- Landing target, legs, return commands, restart readiness and controller gains
  are declared project assumptions. Wind/parameter robustness and physical
  TVC/RCS/contact dynamics are not established by one successful run.

Mission references:

- [SpaceX CRS-10 press kit, mission timeline, PDF p.2](https://www.nasa.gov/wp-content/uploads/2018/07/crs10presskitfinal.pdf): approximate preflight schedule.
- [SpaceX CRS-10 technical webcast](https://www.youtube.com/watch?v=rUDLxFUMC9c): primary flight record.
- `data/crs10_ascent_reference.csv` and `data/crs10_ascent_reference_source.json`:
  video-derived ascent data and provenance, classified E rather than official
  onboard telemetry. The ascent operating point and pitch program were calibrated
  against this reference, so the comparison is not independent validation.
