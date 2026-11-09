# CRS-10 model audit and recovery design

## Status and execution

`python main.py --no-show` runs ascent, shutdown, separation and booster coast
to 160 s. The booster is still climbing at this point. The terminal case is an
independent fixture; it is not spliced into the ascent trajectory.

`python main.py --case full --no-show` attempts the continuous return with
project-designed guidance. It uses the inherited state at every phase boundary;
no position, speed, attitude or fuel reset creates a landing. A failed contact
is recorded as `landing_failed`, not relabelled as success. Current development
trials fail landing criteria. This is an experimental recovery controller.

The existing C sensor/filter algorithms are unchanged. A revised 100 Hz IMU /
10 Hz GNSS/LiDAR run through 160 s completed, seed 1, with 0.618 m saved-sample
position RMS at the IMU station. This is an internal metric under assumed noise
and truth attitude, not real Falcon navigation accuracy. The fetched Emilien-Branch has no unmerged commits relative to main as of
9 November2026. Further C branch changes are not presumed. Closed-loop
navigation feedback remains a distinct model-design choice.

## Why the previous speed was 44% high

The comparison uses geometric-CG ground-relative speed, which is invariant to
rotation of the Body basis. It is simulation truth, not a Kalman estimate.

| Audit | Evidence and conclusion |
| --- | --- |
| Input data and units | Original JSON was downloaded again: SHA256 matches `data/crs10_ascent_reference_source.json`. Author specifies velocity m/s and altitude km. Our export only converts height ×1000; it does not convert velocity a second time. Raw-file stage mixing remains a source risk, but the comparison ends before stage separation. No new frame-by-frame video verification is claimed. |
| Effective propulsion level | Old per-engine vacuum thrust 914 kN implied 845.099 kN at sea level; 9 engines gave 7.606 MN. S10's historical sea-level value is 756 kN, or 6.804 MN total. The old full-throttle effective operating level was 11.8% above that historical total. |
| Early divergence | At 10 s, old speed was 43.2 versus reference 28.3 m/s. At 60 s, 407.2 versus 281.9 m/s. The discrepancy was present early, not introduced by separation or sensors. |
| Fuel coupling | Old full flow 2697.17 kg/s left only 31948 kg at MECO. Revised effective flow 2434.24 kg/s leaves 68889 kg. Greater old consumption reduced mass faster and amplified later acceleration. |
| Attitude program | Old arbitrary pitch τ=60 s produced too much horizontal speed. Using τ=100 s shifts more impulse vertically. This value is calibrated against the same ascent data, not measured guidance. |
| Gravity/rotation/drag | In the old 1–139.9 s trace, projected velocity contributions were thrust +3418.0, gravity −1010.7, aero −100.7 and frame acceleration +2.7 m/s. These are contributions along the actual trajectory, not counterfactual runs. Rotation cannot explain a 708 m/s excess; increasing drag purely to force a match is not justified. |
| Integrator | Prior 0.05 versus 0.01 s check gave <0.001 m final altitude difference. It does not explain the large physical mismatch. |

Source-pinned numerical audit: `data/ascent_model_audit.json`.

| Controlled sensitivity at 139.973 s | Speed | Error versus 1609.444 m/s | Altitude |
| --- | --- | --- | --- |
| Previous 914 kN / pitch τ=60 s | 2317.0 m/s | +44.0% | 57.53 km |
| Historical effective thrust only / τ=60 s | 1731.9 m/s | +7.6% | 38.58 km |
| Historical effective thrust / τ=100 s | 1574.2 m/s | −2.2% | 54.63 km |
| Reference webcast readout | 1609.444 m/s | — | 59.40 km |

The isolated sensitivity cases retain the previous eastward heading to avoid
mixing heading with thrust/pitch changes. The revised team default uses a 45°
northeast heading as an explicit project assumption. Small rotation differences
do not change the main conclusion. The revised ascent is still about 8% low in
height near 140 s. Better agreement is calibration, not independent validation.

**Identifiability:** this does not prove a particular engine hardware rating was
wrong on CRS-10. A higher-rated engine at lower throttle can give the same
operating point. Actual throttle, loaded mass and engine condition are not
uniquely recoverable from this speed curve. The defensible finding is that the
old effective thrust/consumption plus assumed pitch history were inconsistent
with this reference trajectory.

## Updated operating-point evidence

**S10:** [SpaceX Falcon 9 Payload User Guide Rev.2, 21 October 2015, p.11](https://device.report/m/08e51aa8dc266d701619331c9fdf0c23ce8456f61c08bd14084c00cc7c4d36ae.pdf),
section2.2, historical sea-level thrust 756 kN per engine. SpaceX authorship;
mirror hosting. [NASA SSP51071 section6.1](https://www.nasa.gov/wp-content/uploads/2017/06/ssp_51071-baseline_0.pdf)
independently identifies the original 2015 guide URL.

The nominal model adopts 824901 N vacuum equivalent = 756000 + 0.68×101325 N.
This is **E/A for CRS-10**, not a certified mission rating. Isp311 s, nozzle
area0.68 m² and minimum throttle0.60 remain mixed-source estimates from the
previous literature set. Original masses are unchanged. See `parameter_evidence.md`.

**M1:** [Official CRS-10 press kit](https://www.nasa.gov/wp-content/uploads/2018/07/crs10presskitfinal.pdf),
p.2, explicitly approximate preflight schedule: MECO141, separation144,
boostback start161, entry start392 and landing start453 s. Adoption of exact
ascent/separation clock commands is A. Entry/landing in the recovery design use
state conditions instead; their published times remain comparison markers.

**M2:** [SpaceX technical webcast](https://www.youtube.com/watch?v=rUDLxFUMC9c).
**M3:** [Pinned telemetry author documentation](https://github.com/shahar603/Telemetry-Data/blob/b245d3b81aa36b7941ec10f3f4b508999d106a6d/README.md).
The author's estimated MECO143 s is not silently mixed with the official planned
141 s. No post-separation state is copied from upper-stage telemetry.

## Continuous recovery design — all project choices are A

| Phase | Commands | End condition and protection |
| --- | --- | --- |
| Ascent | 9 engines, historical effective level; northeast ideal pitch τ=100 s | Smooth shutdown starts140 s, ends141 s with zero flow/derivative |
| Coast / separation | Engines off; inherit state | Separation144 s; zero-impulse/no-discard map audits mass, momentum, energy and structural continuity |
| Flip / booster coast | B booster aero and bounded ideal attitude moments; point for reverse-horizontal impulse | Boostback161 s only if detached, engines off/zero flow, attitude ready and complete return budget fits |
| Boostback | Centre plus diametrically opposite engines [0,1,5], throttle1; predictor includes spherical gravity, rotation, B force coefficients and prospective entry | Predicted miss ≤500 m after ≥10 s, or 38 t budget /60 s cap; 1 s smooth shutdown. Reaching a budget cap is recorded separately from reaching the target. |
| Return coast | Grid fins deployed, legs stowed via `booster_return`; attitude control uses B fin effectiveness | Entry on descending crossing of45 km and restart guards, not a forced392 s state |
| Entry burn | Same 3 engines, throttle0.8; bounded force/attitude allocation includes aerodynamic forces | Speed ≤500 m/s or16 t budget; smooth shutdown. Goal/cap reason recorded. |
| Aero descent | Engine off; fin/attitude commands correct predicted lateral motion | Predict full-thrust stopping trajectory including altitude-dependent drag, mass decrease and startup ramp; ignite when predicted clearance reaches1500 m allowance |
| Landing burn | Centre engine, target throttle0.60–1; legs deployed; net-force allocation prioritizes vertical braking | Guidance updates0.2 s, accelerating to0.01 s near contact; thrust plus B aero determines tilt/throttle. Stop free-flight integration at first leg contact. 10 t planning allocation may use unused earlier budget, while final reserve is protected. |
| Contact assessment | No velocity/attitude reset | Pass all contact speed, tilt, rate and position limits; otherwise `landing_failed` |

Initial return budget: boostback38 t + entry16 t + landing10 t + reserve2 t =66 t.
Revised MECO inventory is68.889 t; the model's unburnable limiting-tank residual
is about1.202 t, leaving67.687 t available. The 2 t reserve is additional to that
residual. Budgets are project allocations, not SpaceX flight records. Unused
boostback/entry allocation can transfer to landing; the reserve cannot.

Restart health/settling is explicitly assumed ready after zero-flow and attitude
checks. The 1 s ignition/shutdown ramps are numerical transients, not measured
engine data. Engine groups are symmetric so the existing balanced exhaust
contract remains valid; no unsupported asymmetric-port gimbal force is inserted.

## Guidance, actuators and contact

The landing command uses a ZEM/ZEV-type kinematic terminal law with target
position rₜ and target velocity vₜ (vertical component −1 m/s). The lateral terminal horizon is half the vertical horizon so position correction
occurs earlier. Required net acceleration is
6(rₜ−r)/τ² −4v/τ −2vₜ/τ. Gravity/frame effects are removed before requesting
nongravitational force. A bounded search allocates **thrust plus aerodynamic
force**, not thrust alone, to this request. Vertical-force error receives100:1
weight relative to horizontal error. This is our project controller, not G-FOLD
or SpaceX software.

Attitude control allocates B's linearized four-fin moments with ±0.35 rad
travel; remaining ideal moment is capped at2 MN·m with powered engines and
0.2 MN·m in coast. Maximum desired angular rate is0.22 rad/s and angular
acceleration0.12 rad/s². These are A-class virtual actuator limits. Landing tilt authority tapers from20°
to2° in the last1000 m. Physical
TVC/RCS force coupling, valve dynamics and fin slew-rate limits are not modeled.
Fins deploy after boostback; legs deploy at landing ignition. Deployed-leg mass
redistribution and contact reaction forces are not resolved.

The target is an A-class rounded offset14 km south /6 km east of the launch
site, on the same spherical Earth. It is not a surveyed LZ-1 coordinate. The
contact geometry is four assumed leg tips1 m below the structural base, with
3.5 m radius. Contact uses the lowest tip, not CG altitude zero. Trials are
localized by step reduction and a contact tolerance; no penetration is reset
into an artificial successful landing.

Landing criteria (A): vertical speed ≤2 m/s, horizontal speed ≤1 m/s, tilt ≤5°,
angular rate ≤1°/s, target error ≤10 m, legs deployed and actual contact.
Current prototype trials do not satisfy these criteria. Do not loosen the
criteria after a failed run to manufacture success.

## Output and handoff

Each run saves truth/measurement/navigation CSVs, accepted events and separation
audit, summary with plan/sources/input hashes, ascent comparison CSV and
`overview.png`. Results remain ignored by Git. `--case full` returns a nonzero
exit code unless landing succeeds. Default ascent can complete its shorter
requested interval without claiming a complete return.

The runner uses the same fixed first-stage IMU/GNSS/LiDAR stations across
separation. Random sensor generation occurs only at accepted sampling epochs.
The recovery controller currently uses truth-state feedback; connecting C's
estimates to guidance is a distinct final integration step. Its current filter
also receives truth attitude. Neither fact should be hidden in the report.

Before submission: finish recovery tuning/robustness or explicitly limit claims
to demonstrated phases; check any newly published C changes; run complete sensor/filter
integration and error analysis; confirm clean installation/reproducibility;
write the manual, assumptions/sources, results and dangerous-error discussion.


## Numerical and runtime verification

The Earth-centred RK4 safeguard is described in `powered_dynamics.md`.
Body state storage is retained; intermediate position/velocity integration is
performed externally so that rotating the basis cannot invent a below-ground
trial position. The high-radius/free-spin regression checks every trial.

Default sensor scheduling remains100 Hz IMU /10 Hz GNSS/LiDAR. For a faster
engineering integration check use:

```text
python main.py --case full --imu-rate 20 --dt 0.05 --no-show
```

This explicitly records20 Hz IMU in the run plan/specs. It does not change C's
sensor/filter implementation or establish a final hardware rate. Final team
validation should also run the agreed100 Hz setting. A successful unit suite
or a complete saved trajectory does not imply a successful landing.


## Latest continuous integration result

The development run used0.05 s maximum dynamics step,20 Hz IMU and10 Hz
GNSS/LiDAR, seed1, with the Earth-centred RK4 safeguard. It visited all phases
in one trajectory and preserved the same navigation instance through separation.
No state was reset to create a nominal return.

| Event | Simulated time [s] |
| --- | --- |
| MECO | 141.000 |
| Separation | 144.000 |
| Boostback start / end | 161.000 /208.850 |
| Entry start / end | 370.200 /393.050 |
| Landing burn start | 426.300 |
| Leg contact | 465.069 |

These are simulated outcomes of the project controller, not recovered actual
CRS-10 event states. Boostback ended on its approximate planning cap rather
than the500 m miss goal. Per-burn caps are sampled decisions followed by a
smooth shutdown, so a small overrun is possible; observed boostback consumption
was38.015 t. A strict per-burn cutoff allowance remains a refinement.

| Contact metric | Actual | Criterion | Result |
| --- | --- | --- | --- |
| Vertical speed | 2.625 m/s | ≤2 m/s | Fail |
| Horizontal speed | 1.0015 m/s | ≤1 m/s | Fail, including the small excess |
| Tilt | 1.946° | ≤5° | Pass |
| Angular rate | 0.289°/s | ≤1°/s | Pass |
| Target error | 17.943 m | ≤10 m | Fail |

Status is `landing_failed`; `full_mission_completed` is false. Remaining fuel
was8.979 t, above the declared reserve plus limiting-tank residual. Consumption:
boostback38.015 t, entry14.184 t and landing7.712 t. Fuel shortage did not cause
this last failed contact; terminal guidance accuracy remains the principal
unfinished control task.

The run generated9302 IMU,4651 GNSS and4651 LiDAR samples;492 LiDAR updates
were accepted. Navigation position RMS over saved samples was1.154 m with
truth attitude and provisional noise. This checks integration, not mission
accuracy or estimator robustness. The earlier count is superseded by
**144 passing regression tests** after the geocentric-trial regression was added.

A final short checkpoint-only tuning experiment was not promoted: it did not
meet all acceptance criteria. The shared parameters remain those of the full
recorded integration run. All generated test trajectories stay under ignored
`results/`, not the shared repository history.


## Landing-controller revision under verification

The acceptance limits and landing target remain unchanged. The vertical loop
now uses the actual lowest leg-tip clearance hₗ and radial velocity wₗ, with
wₜ=−0.5 m/s and descent braking acceleration (wₗ²−wₜ²)/(2hₗ). Rotational,
CG-migration and radial-curvature terms convert the contact-point request to a
CG acceleration request. This removes the former CG-target/contact-height mismatch.

Horizontal control is a critically damped fixed-gain loop:
a⊥=−0.0576e⊥−0.48v⊥. Its natural rate is0.24 rad/s; it does not increase gains
without bound as touchdown approaches. Landing attitude gains are1.5 and6,
with the original0.22 rad/s desired-rate,0.12 rad/s² acceleration and2 MN·m
virtual moment limits retained. Finite-travel grid fins supply residual moment
after virtual TVC authority during landing; coast allocation is unchanged.

Landing tilt tapers to4°, within the existing5° acceptance limit. Below20 m
base altitude the last steering reference is held to let angular motion settle;
actual attitude, velocity and position continue to evolve through dynamics.
There is no state reset or replacement of the touchdown metrics.

The fetched remote C branch `Emilien-Branch` has zero commits absent from main.
No additional C merge is currently required. This statement refers to published
GitHub code, not any unpublished work on a teammate's computer.
