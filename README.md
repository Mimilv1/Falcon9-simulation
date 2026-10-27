# Falcon9-simulation

Project 2 A-block with a first replaceable B model: powered dynamics, sourced
vehicle data, and an executable early-ascent segment exporting ideal truth for C.

## Current team convention

All shared position, velocity, force, moment, gravity and angular-rate vectors
use **current Body components**. **Body +x points toward the nose**. The frame is
right handed; engine-ring angle zero is +y and +90° is +z. Geographic heading,
ENU and ECEF transforms are deferred.

The quaternion is Hamilton [w,x,y,z], mapping current Body to the fixed **initial
Body basis**. Initial q=[1,0,0,0]. It records attitude relative to the beginning
of the run; no east/north/up direction is prescribed. This mathematical initial
basis is necessary for orientation bookkeeping, not a geographic interface.

`State.position_body_m` means the vector **from a fixed initial origin to CG,
expressed in current Body axes**. It is not CG position relative to itself.
Accordingly, both position and velocity derivatives include the rotating-basis
terms: ṙᴮ=vᴮ−ω×rᴮ and v̇ᴮ=aᴮ−ω×vᴮ. Dropping the first term's rotation correction
would make a stationary point move when only the rocket turns.

## Files

| File | Responsibility |
| --- | --- |
| `VehicleData.py` | Engine/Stage/Vehicle objects, staging, base/CG mapping, accepted-state synchronization |
| `vehicle_parameters.py` | Single evidence-led numeric input profile |
| `data/vehicle_parameters.json` | K/E/A metadata, source URLs, conditions and limitations |
| `vehicle_mass_model.py` | Fixed-O/F liquid inventory and mass-profile adapter |
| `mass_properties.py` | Components, CG/inertia and first/second geometric derivatives |
| `mass_flow.py` | Exhaust-port boundary data, exit gyration and balance validation |
| `load_interface.py` | Forces at structural stations -> net force/current-CG moment |
| `dynamics.py` | Open-system CG accelerations and RK4; Body-only shared states |
| `propagation.py` | Smooth-segment propagation, scheduled breakpoints and fuel-limit stopping |
| `frames.py` | Generic vector/quaternion operations relative to initial Body; no Earth conversion API |
| `Environement.py` | Standard troposphere, constant-magnitude gravity, prescribed wind |
| `propulsion.py` | Steady engine fuel flow and ambient-pressure-dependent thrust |
| `Aerodynamics.py` | Digitized Mach-dependent drag curve; drag through CG |
| `ascent.py` / `main.py` | Connect A+B, propagate one phase segment, export C truth |
| `prescribed_loads.py` | Explicit smooth force doublet for powered-attitude verification |
| `engine_schedule.py` | Analytic smooth shutdown amplitude and derivative |
| `staging.py` | Accepted-state separation map and conservation audit |
| `phase_sequence.py` | Shutdown, attached coast, split and two gravity-only coast branches |
| `config/powered_ascent.json` | Shared first-segment inputs and explicit assumptions |

Read [the shared interface](docs/core_api.md), [powered dynamics and assumptions](docs/powered_dynamics.md),
and [the English source table](docs/parameter_evidence.md).

## Use

Python 3.11+ and NumPy are required. Install with `python -m pip install -r requirements.txt`.
Run `python -m unittest discover -s tests -v` for analytical and integration checks.
The fixtures are diagnostic mathematical inputs, not measured flight data.

Run `python main.py` for **30 seconds of early vertical powered ascent**.
Results appear in `results/powered_ascent/`, including `truth_for_c.csv`, the
configuration, source metadata and summary. See [the first segment guide](docs/first_ascent_segment.md)
for assumptions, signal meanings and how C can use the output.

Run `python main.py --config config/ascent_attitude_check.json` for the connected
small-angle rotation check, writing `results/attitude_check/`. It adds an explicit
verification force and a second virtual accelerometer; see [the attitude check](docs/attitude_check.md)
for assumptions, moment signs and sensor lever-arm validation.

Run `python main.py --config config/shutdown_separation_check.json` to connect
smooth shutdown, zero-impulse separation and two tracked coast branches. The
[sequence guide](docs/shutdown_separation.md) distinguishes the explicit test
timing and gravity-only/frozen-liquid assumptions from real flight behavior.
Separate CSV streams and `events.json` appear in `results/shutdown_separation/`.

`build_falcon9()` creates independent team objects. `vehicle.mass_model()` evaluates
RK trial inventories without mutating Stage objects; use configuration=None.
`vehicle.to_state(consumption_by_stage=...)` includes CG migration when mapping
base position/velocity to CG. After an accepted state, `vehicle.apply_state(...)`
updates the vehicle and stage fuel inventories; do not synchronize RK trial states.
Pass consistent flow and derivative data to both mappings when flow is nonzero.

## Current scope

Powered integration now updates propellant mass, CG, inertia, translation and
attitude together with explicit exhaust momentum terms. Supply ExhaustPort data
for every positive consumption rate. Supported exhaust centroids and CG migration
are axial; off-axis unbalanced outflow is explicitly rejected.

Flow must be smooth with a consistent dq/dt inside each segment. Smooth ignition
and cutoff are supported. Abrupt mass-flow jumps require a separate impulse/state
map and are not silently approximated as smooth. `propagate` stops at a limiting
tank, within its stated tolerance, rather than clipping or burning negative fuel.

Settled flat cylindrical liquids, fixed density, fixed loaded ratio/consumption
O/F and omitted internal fluid swirl remain assumptions. Venting, variable O/F,
slosh, free-fall liquid behavior and full propellant feed-system dynamics are not
implemented. The first replaceable B provider supplies thrust, drag and gravity.
C must provide a fixed Body
sensor station to obtain physical IMU acceleration; CG coordinate acceleration
is a different quantity when the CG moves.

The shared branch includes the project phase configuration and runnable main;
generated results and separate diagnostic/plot scripts are not committed.
Actual-flight validation, complete GNC and geographic output conversion remain
future work. This is a mixed-source research model, not a reconstructed mission.
