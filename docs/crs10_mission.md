# CRS-10 execution and team handoff

## What is connected

`main.py` → `mission.run_mission` → vehicle/mass properties → B loads → A RK4
→ accepted state → existing C sensor measurements and filter → CSV/JSON/plot.
Every RK4 trial rebuilds mass, CG, inertia, flow, gravity, rotation effects and
aerodynamics from the trial state. Random sensor samples are generated only at
accepted sampling epochs. The filter uses its own estimated position/velocity
for environment propagation, with truth attitude supplied as in C's current code.

There is one CRS-10 vehicle parameter set. The mission plan is
`data/crs10_mission.json`; it contains commands and explicit assumptions, not a
second copy of vehicle masses. No new comments/docstrings were added to the code.

## Source-timed events and implemented scope

**M1, K:** [SpaceX CRS-10 press kit](https://www.nasa.gov/wp-content/uploads/2018/07/crs10presskitfinal.pdf),
February 2017, PDF p.2, “Mission Timeline (all times approximate)”. MECO 141 s,
stage separation 144 s, boostback start 161 s, entry start 392 s, landing start
453 s are official approximate preflight milestones. Adopting them as exact
simulation commands is **A**, not knowledge of actual flight-computer thresholds.

**M2:** [SpaceX technical webcast](https://www.youtube.com/watch?v=rUDLxFUMC9c).
No new frame-by-frame event timing was extracted in this revision.

**M3, E:** Existing `data/crs10_ascent_reference.csv` is published webcast
readout, pinned author revision and processing in
`data/crs10_ascent_reference_source.json`. It contains 140 samples in
0–139.973 s. It is independent ascent reference data, not a source for a
post-separation booster state. The speed reference frame is not independently
verified. An author's estimated 143 s MECO is not mixed with M1's planned 141 s.

| Phase / event | Initial state | Transition / success check | Current implementation |
| --- | --- | --- | --- |
| Powered ascent | Full stack, structural base at external [0,0,R], base velocity zero relative Earth, upright; initial CG velocity includes fuel-driven CG migration | Reach 140 s with valid fuel/model state | 9 engines, full throttle; ideal eastward pitch program; changing mass/CG/inertia; B full_stack aero |
| Engine shutdown | Inherit the accepted ascent state at 140 s | Smooth quintic envelope reaches zero flow and derivative at 141 s | 1 s duration A; thrust and flow scaled consistently; ideal pitch control retained |
| Full-stack coast | Inherit shutdown state; all engine commands off | Reach 144 s with zero flow/derivative and full stack attached | Gravity, rotation and aero continue; no zero-velocity reset; no control couple after MECO |
| Separation | Inherit accepted coasting state at 144 s | No-impulse/no-discard map; mass/momentum/energy audit and structural continuity pass | Child CG states calculated, only lower stage propagated; first-stage sensors/filter preserved; aero switches to booster_alone |
| Booster coast | Inherit actual lower-stage state | Reach requested 160 s endpoint without model/fuel failure | Free rotation under B aero; no flip command is claimed |
| Boostback | Must inherit preceding state, detached, engines off, settled/ready and sufficient declared fuel budget | Scheduled reference 161 s plus physical restart guards; burn end unresolved | Not implemented; a timer alone cannot restart a burn |
| Return coast / entry | Must inherit a completed boostback and a physically descending state for entry | Entry reference 392 s plus guards; entry-burn end unresolved | Not implemented |
| Aero descent / landing burn | Must inherit completed entry, then descent | Landing-start reference 453 s plus guards; touchdown from leg-contact geometry | Not implemented as continuous return |
| Touchdown | Pre-contact state at actual descending contact crossing | All declared speed/tilt/target/leg limits satisfied | Criteria helper only; no contact detector/controller in default provider |

Initial attitude is Hamilton [√0.5,0,−√0.5,0] Body → launch-aligned ECEF.
The launch latitude is 28.6° A, Earth rotation 7.292115×10⁻⁵ rad/s. The external
origin is Earth centre and the initial base lies on external +z. Body +x is the
nose direction. Structural stations remain distances from the vehicle base.

The ideal pitch program is retained from team main: begins at 8 s,
τ = 60 s, asymptotic tilt 80°, PD gains 2 and 2.5, unlimited CG pure couple.
It is **A**, not measured CRS-10 guidance or a physical TVC actuator. No
full-mission match can be inferred merely from a successful numerical run.
B's full-stack Cd originated with a generic fairing; high-Mach and booster
coefficients, side-force terms and rotary damping contain declared estimates.

## Independent terminal-descent fixture

`python main.py --case terminal --no-show`

This uses the same first-stage parameters, B environment/aero, A dynamics and C
sensor/filter implementation. It starts a separate simulation; it does not
reset or splice an ascent state to claim a continuous CRS-10 trajectory.

| Setting | Value | Class and rationale |
| --- | --- | --- |
| Structural-base altitude | 1500 m | A; within the current C LiDAR 1800 m envelope |
| Initial base vertical velocity | −100 m/s | A; descending, for braking/sensor integration |
| Propellant | 15000 kg | A; above limiting-tank residual, enough for this short check; not a measured landing reserve |
| Initial attitude / angular rate | Upright / zero | A; tail-down powered descent fixture |
| Active engine | Central engine only, throttle 0.65 | A; symmetric axial thrust with no invented gimbal command |
| Aerodynamics | booster_landing | B model with fins/legs assumed deployed from start; no active fin control |
| Interval | 10 s | A; exercises descent and deceleration without contact or restart transients |

Its success criterion is integration/interface operation: finite state,
positive fuel margin, decreasing downward speed, valid LiDAR updates, and saved
truth/measurements/estimates. It is not a touchdown-success test. A complete
landing needs contact geometry, target and a control law.

## Run and outputs

```text
python main.py --output results/team-ascent --no-show
python main.py --case terminal --output results/team-terminal --no-show
python main.py --plot-only results/team-ascent
python main.py --until 20 --no-show
```

| Output | Meaning |
| --- | --- |
| truth.csv | 0.1 s accepted-state samples plus endpoint: phase, Body/ECEF CG position/velocity, quaternion, angular rate, fuel, mass, CG, inertia, forces/moments, gravity/frame acceleration, altitude, Mach, thrust, q̄ |
| measurements.csv | 100 Hz IMU and 10 Hz GNSS/LiDAR, validity flags; invalid numeric entries remain NaN |
| navigation.csv | IMU-station ECEF estimates, accelerometer bias, errors relative to that same structural point, state sigmas |
| events.json | Accepted events/guards; at separation, parent/child states and conservation audit |
| validation.csv | Time-aligned simulation versus ascent reference, only overlapping times; no extrapolation; absent for terminal case |
| summary.json | Scope/status, seed, actual sensor specs, plan/sources, input hashes, final state, internal navigation metrics and independent reference-error metrics |
| overview.png | Altitude/speed/reference, fuel, CG station and navigation errors |

Sensors remain at current C-main first-stage stations: IMU 38 m, GNSS 41 m,
LiDAR 0.5 m from base. IMU rate is fixed at 100 Hz instead of changing with
integration step. Existing C-main biases/scale/random-walk settings are
retained explicitly in the plan. IMU errors and sensor stations are assumptions,
not Falcon 9 certified data. Navigation truth attitude remains an explicit limit.

One `NavigationRun` instance survives separation: do not initialize the filter
from booster true CG. The estimated point remains the fixed IMU station.
Sampling, outputs and events are aligned exactly; dt >0.01 s is rejected when
navigation is enabled. Dynamics-only convergence runs can use a larger dt.

Results are ignored by Git. Use an empty output directory. Noncompletion saves
partial accepted results, marks `full_mission_completed=false` and exits with
code 2. Requesting >161 s cannot silently skip an unsupported restart.

## Checks and what remains

The regression suite covers moving-CG launch initialization, event ordering,
zero-flow guards, fuel budget checks, no skipped unknown burn endpoints,
separation conservation/structural continuity, filter preservation, sensor
sampling, model-domain stops, and independent terminal setup. B's ideal
vertical tests explicitly disable Earth rotation; the production run retains it.
Altitude test tolerance is 1 μm where Earth-radius subtraction causes nm
roundoff. Neither adjustment changes production physics.

For the final team run, first pull main, install requirements, and use the two
commands above. A/B must still address the independent ascent-reference
mismatch, justify a physical pitch/throttle history, and decide which flight
phases the report claims. Continuous boostback/entry/landing additionally needs
burn budgets and goals, attitude/flip commands, restart envelopes and contact
geometry. Do not claim those phases from the independent terminal fixture.
C's algorithm design remains C's work; this merge only adapts its existing API
and writes the measurements/truth needed to verify the final integration.


## Verified run, 9 November 2026

Execution source commit: `986fc50` (later edits only remove stale frame comments
and update documentation). All 100 core/mission tests and 36 environment tests
passed. Default 0.01 s integration with fixed 100 Hz IMU, seed 1:

| Check | Result |
| --- | --- |
| Ascent → shutdown → separation → booster coast | Completed to 160 s, no model-domain or fuel stop |
| Event times | Shutdown start 140 s, zero-thrust MECO 141 s, separation 144 s |
| Final CG altitude / ground-relative speed | 74.03 km / 2271.9 m/s |
| Final booster mass / propellant | 54148.2 kg / 31948.2 kg |
| Coast | Zero thrust and constant propellant after 141 s |
| Separation | All mass/momentum/energy relative residuals below 1×10⁻¹⁰; structural position continuity error below 1×10⁻⁸ m |
| Samples | IMU 16001; GNSS 1601; LiDAR 1601 generated samples, including invalid out-of-range returns |
| Internal navigation position RMS | 0.690 m over saved 0.1 s samples; truth attitude and provisional sensor errors |
| Independent terminal check | Base altitude 1500 → 777.6 m; vertical velocity −100 → −44.69 m/s over 10 s; LiDAR updates 101 |
| Numerical step comparison | 0.05 versus 0.01 s, same event times; final altitude difference <0.001 m and speed difference <0.00001 m/s |

Independent reference comparison exposes a substantial modeling mismatch:
at 139.973 s the reference altitude/speed are 59.4 km and 1609.444 m/s;
the simulation gives about 57.53 km and 2317.03 m/s. Speed is about 44% high.
Across the 140 reference samples, altitude RMSE is about 2.402 km and speed
RMSE 297.1 m/s. Simulated maximum q̄ is about 35.7 kPa at 63.4 s, compared with
the approximate 75 s preflight marker. Numerical consistency is therefore
verified, but CRS-10 trajectory fidelity is not. Review the assumed pitch and
full-throttle history, reference-frame interpretation and B-model applicability
before any physical tuning or claim of validation.

`verification.json` in the ascent result directory records runtime assertions
and the step comparison. Generated CSV/JSON/PNG files remain local under
`results/`; team members regenerate them with the documented commands.
