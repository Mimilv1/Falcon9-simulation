# Smooth shutdown and zero-impulse separation check

This connects a sequence of accepted flight states through engine shutdown,
attached coasting and stage separation, then propagates **both** resulting
vehicles. It uses the existing Body +x convention and powered dynamics. The
schedule is deliberately an **A-class event verification**, not Falcon 9 mission
timing, an engine-transient reconstruction or a complete coast model.

```text
python main.py --config config/shutdown_separation_check.json
```

Default output is `results/shutdown_separation/`. Existing vertical and attitude
configurations keep their separate output folders.

## Schedule and assumptions

| Interval / event | Behavior | Class / basis |
| --- | --- | --- |
| 0–30 s | Nine first-stage engines at the existing steady nominal command; the 10–14 s attitude-check doublet remains enabled | Engine inputs E from the existing source profile; commands and diagnostic doublet A |
| 30–32 s | Smooth amplitude reduces total thrust and propellant flow to zero | A: 2 s analytic verification envelope; no measured engine response supplied |
| 32–34 s | Attached-stack coast, with the existing full-stack drag surrogate and zero propulsion | A: frozen liquid geometry; not a fluid settling model |
| 34 s | Zero-impulse stage separation; no discarded mass | A: selected test time, no pneumatic/pyrotechnic/spring impulse or interstage shedding |
| 34–36 s | Booster and upper stack propagated separately under gravity only | A: no separated-vehicle aerodynamic data adopted; liquid distributions frozen, no slosh or re-settling |
| Split station | Upper-stage base at original axial station 42.6 m | Existing E/A contiguous-stage geometry; no new measured separation-plane dimension |
| Sensors | Original base sensor stays on booster; upper stack has an independent virtual sensor at its own base | A: arbitrary fixed verification stations |

The 30/32/34/36 s times are configurable and deliberately early. They are not
assertions about real Falcon 9 MECO, separation or second-stage ignition. Much
first-stage fuel remains at this test separation. The upper-stage engine stays
off throughout. Both pieces retain the same fixed initial origin and initial
attitude reference; only their own structural station datums differ.

The literature data, K/E/A classifications and links remain in
[the parameter evidence](parameter_evidence.md) and
[the first-ascent model guide](first_ascent_segment.md). This step introduces
verification assumptions rather than new claimed vehicle measurements.

## Shutdown envelope and mass derivatives

Define u = (t − tₛ)/Δt, where tₛ = 30 s and Δt = 2 s. During 0 < u < 1:

- a = 1 − 10u³ + 15u⁴ − 6u⁵.
- da/dt = −30u²(1 − u)²/Δt.

Before the window, a = 1 and da/dt = 0. Afterward, both are zero. Amplitude and
its first two derivatives join continuously at both endpoints. The implementation
uses an algebraically equivalent factored form for numerical stability.

For nominal engine consumption Q and nominal ambient-corrected thrust Tₙ:

- q = aQ and dq/dt = Q·da/dt.
- T = aTₙ, with Tₙ = τTᵥ − Aₑp.

Here q is positive outward mass flow in kg/s, Q is the existing steady flow,
τ is the fixed nominal throttle, Tᵥ is vacuum thrust, Aₑ is nozzle exit area and
p is ambient pressure. The entire corrected thrust is scaled, so an off engine
does not leave an artificial negative Aₑp force. This envelope is not a claim of
steady operation below the adopted minimum throttle.

The analytic dq/dt enters the exhaust-port contract, mass second derivative and
CG migration acceleration. It is not silently set to zero during shutdown.
The total effective full-flow duration is 30 + 2/2 = **31 s**. This independently
predicts the remaining propellant and is checked against integration. Breakpoints
split the integrator exactly at scheduled transitions. Trial evaluations never
mutate engine commands or accepted Stage inventories.

## Separation state map

The map is allowed only after flow and its derivative are zero and engine
commands are off. It synchronizes a copy of the accepted parent state and uses
the existing `Vehicle.separate` geometry. No velocity kick is invented.

For a rigid material point at displacement r from a reference point, its velocity
is v + ω × r. Therefore upper- and lower-stage base velocities can differ when
the stack rotates. Copying the same CG velocity into both vehicles would lose
this kinematic relationship.

Each child receives its own CG, mass and inertia. The sum of child CG momenta
must equal parent momentum. Angular momentum includes **both** each child's spin
and the orbital contribution from its CG motion about the parent CG. Checking
only the spin terms would incorrectly report a loss. The audit also checks:

- Total mass and recombined CG position.
- Positions and velocities of corresponding structural base points.
- Attitude and angular-rate continuity.
- Kinetic energy, for this zero-impulse/no-discard map.

The map rejects a failed audit. `events.json` contains absolute and relative
residuals and the event assumptions. Relative conservation tolerance defaults
to 10⁻¹⁰; structural continuity is checked within 10⁻⁸ m and m/s. These are
numerical verification tolerances, not physical sensor accuracies.

## C's files and event meanings

| File | Contents |
| --- | --- |
| truth_stack.csv | Attached stack, t = 0–34 s |
| truth_booster.csv | Separated first stage, t = 34–36 s |
| truth_upper_stack.csv | Separated upper stage plus fairing/payload, t = 34–36 s |
| truth_for_c.csv | Combined long table with vehicle_id, phase, event_side and imu_id |
| events.json | Shutdown times and full separation conservation audit |
| summary.json | Cutoff state, final state for each vehicle, assumptions and audit |
| manifest.json / run_config.json | Code/source hashes, configuration and sensor station metadata |

At t = 34 s the combined table deliberately has three rows: parent `stack/pre`,
`booster/post` and `upper_stack/post`. These are different objects at the same
event time. Filter by vehicle_id, and retain both event sides where appropriate.
The CG position/velocity columns track different physical points after splitting;
their change is not a teleportation of a fixed material point.

The original sensor's imu_id continues from stack to booster. The upper sensor
has a different imu_id and its own station datum. The upper base has original
station 42.6 m, recorded separately from its own zero station. Comparison-sensor
fields are blank for the upper stack because that virtual sensor remains on the
booster; blanks are not zero measurements. Continue the original sensor stream
using the stack followed by booster records, without feeding upper-stack rows
into the same estimator stream.

Gravity-only coast does not necessarily make every displaced rotating sensor's
output zero: rigid-body tangential/centripetal lever-arm terms still apply. Fluid
motion, drag changes and the separation actuator itself are outside this test.
An acceleration change at the event can occur when loads change; fixed-point
position and velocity remain continuous in the no-impulse map.

## Scope after this step

A can now schedule smooth engine shutdown, map one accepted stack state into two
conserving states and propagate/export both branches. The remaining physical
work is to select justified mission event conditions, supply B's separate-stage
loads and fluid behavior, model any actual separation impulse/discarded mass,
and implement restart/upper-stage ignition with suitable transient data.
This check does not validate real MECO, separation or recovery trajectories.
