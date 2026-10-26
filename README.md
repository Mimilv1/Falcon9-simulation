# Falcon9-simulation

Shared Project 2 core and evidence-led vehicle data, integrated with team commit
`5b61a45b07766fbd458d4410a1034977da824837`.

## Agreed coordinates

| Item | Convention |
| --- | --- |
| Fixed reference | ENU: +x east, +y north, +z up |
| Body longitudinal axis | +x toward the nose |
| Initial vertical Body axes | +x up, +y east, +z north |
| Quaternion | Hamilton, scalar-first [w,x,y,z], Body to ENU |
| Initial vertical quaternion | [0.5, −0.5, −0.5, −0.5] |
| Gravity input | ENU [0,0,−g], supplied by the environment model |
| Inertia order | Ixx axial; Iyy and Izz transverse |
| Structural stations | Body [axial station,0,0], measured from the fixed base datum |

The initial quaternion maps Body [1,0,0] to ENU [0,0,1]. Its four entries are
an orientation representation, not three angles or physical vehicle parameters.

## Files and ownership

| File | Purpose |
| --- | --- |
| `VehicleData.py` | Team EngineSpec/Engine/Stage/Vehicle objects, fresh factory and structural-base staging |
| `vehicle_parameters.py` | Single numeric input profile; declared K/E/A values and derived quantities |
| `data/vehicle_parameters.json` | Input attribute mappings, source URLs, locators, methods and limitations |
| `vehicle_mass_model.py` | Evidence-profile adapter and fixed-O/F first-stage inventories |
| `mass_properties.py` | Generic mass components, axial cylinders, inertia and CG combination |
| `load_interface.py` | Body application-point force to moment about current CG |
| `dynamics.py` | CG-based rigid-body derivatives, RK4, navigation truth |
| `frames.py` | Body/ENU quaternion operations and WGS-84 ENU/ECEF output transformations |

Read [the parameter evidence table](docs/parameter_evidence.md),
[the interface contract](docs/core_api.md), and
[the upstream review](docs/upstream_review.md).

## Use and checks

Python 3.11+ and NumPy are required. Install with `python -m pip install -r requirements.txt`.
Run `python -m unittest discover -s tests -v` for analytical/frame/integration checks.
All fixtures are synthetic diagnostic inputs, not measured flight data.

`from VehicleData import falcon9` preserves the team's entry point.
Use `build_falcon9()` to create independent mutable Stage/Vehicle objects.
`falcon9.get_mass_properties()` returns properties for the currently attached
stages and current `Stage.propellant_mass` values. Engine dry masses are already
included in stage dry masses; do not add them again.

`Stage.propellant_mass` is the authoritative stage-total inventory. First-stage
RP-1/LOX masses are derived with the fixed loaded inventory and consumption O/F.
TankSpec fields are reference geometry/property data, not a second mutable fuel state.
After changing construction inputs, rebuild the vehicle; dry CG/inertia snapshots
must agree with the component mass budget.

`Vehicle.position`/`velocity` refer to the fixed structural base in ENU.
The integrator's `State.position_enu_m` refers to instantaneous CG, and
`State.velocity_body_mps` is CG velocity expressed in Body. At fixed inventory,
`vehicle.to_state()` maps these definitions explicitly. `vehicle.mass_model()`
evaluates trial inventories without changing the shared vehicle. Pass
`configuration=None`; the attached Stage objects define the configuration.
The integrator returns a new State; it does not automatically overwrite Vehicle snapshots.

For B, application points are stage/vehicle structural stations in Body.
Subtract current CG once, then form the cross product with the force. Supply
gravity separately in ENU and do not duplicate it in nongravitational forces.
For C, navigation truth returns ENU CG position/velocity, Body-to-ENU attitude,
and Body rates under the zero-Earth-rotation approximation. A structurally fixed
sensor displaced from CG needs lever-arm treatment.

## Current scope

Demonstration scripts, synthetic flight scenarios and generated trajectory/plot
outputs are excluded from this shared change. The upstream `main.py` still only
advances its loop time; its Vehicle update call remains commented out. Running it
is not evidence of an implemented flight trajectory.

The dual-tank model is a settled-fluid geometric approximation using mixed-source
research inputs. It is not a reconstructed Block 5 mission. Dry-hardware placement,
payload geometry, actuator limits and some engine stations remain explicit assumptions.
The shared core rejects positive consumption with the dual-tank profile until
moving-CG kinematics and exhaust momentum are treated consistently. No powered
dual-tank flight, guidance/control loop, sensor estimator or real-flight validation
is claimed. Earth rotation and curved-Earth propagation remain outside this core.
