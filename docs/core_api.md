# Shared A/B/C contract — Body +x

All shared vectors use current Body axes; +x is axial toward the nose, Ixx is
axial inertia, and Iyy/Izz are transverse. No ENU/ECEF state fields or Earth
coordinate helpers are in the active API. Initial quaternion [1,0,0,0] maps current
Body to the fixed initial Body basis. Geographic heading remains unspecified.

## Origins and state meanings

| Quantity | Origin / meaning |
| --- | --- |
| State.position_body_m | Fixed initial origin -> instantaneous CG, expressed in CURRENT Body axes |
| State.velocity_body_mps | Inertial/geometric CG velocity, expressed in CURRENT Body axes |
| Vehicle.position / velocity | Same fixed origin and current Body axes, but tracking the structural base |
| CG/component stations | Distance from the vehicle structural base, along Body +x |
| AppliedForce.station_body_m | Fixed vehicle structural station, not a CG arm |
| ExhaustPort.station_body_m | Fixed nozzle-exit station from the same structural datum; not necessarily the engine pivot |
| Moment/inertia | About instantaneous CG, Body components |

Fixed-origin position does not become zero because Body axes rotate. Its stored
components change through ṙᴮ=vᴮ−ω×rᴮ. Velocity components similarly require the
transport term −ω×vᴮ. The quaternion records the changing basis.

## B inputs

The first implementation is `SteadyPoweredAscent` in `ascent.py`, connecting
`propulsion.py`, `Environement.py` and `Aerodynamics.py`. It uses the existing
vehicle profile; B can replace each model while keeping the Loads contract.
See [the executable first segment](first_ascent_segment.md).

`provider(time_s, trial_state) -> Loads` is evaluated at every RK4 stage.

| Loads field | Definition |
| --- | --- |
| force_body_n | Net thrust plus aerodynamic force, gravity excluded; thrust already includes relative exhaust momentum/pressure reaction |
| moment_body_nm | Net external moment about current CG |
| gravity_body_mps2 | Gravity in CURRENT Body, added once |
| consumption_kgps | Nonnegative [stage1,stage2] mass consumption |
| exhaust_ports | Each port's stage index, fixed structural station, q, dq/dt and per-unit-mass exit gyration |

Use `assemble_loads` to sum application-point forces and their r×F once.
Supply gravity explicitly. A constant vector in an initial basis generally
changes Body components as the rocket turns; do not hold gravity at [-g,0,0]
through arbitrary rotation. B may use the generic quaternion operations for this
physical direction calculation without introducing a geographic frame.

Port flows must sum to each stage's consumption. The weighted transverse exit
centroid and its derivative must vanish for this axial-flow model. Point ports
use zero local gyration; `uniform_disk_gyration` is an explicit uniform-flux disk
approximation. Swirl and moving/gimbaling exit-boundary kinematics are excluded.

`dq/dt` is the total derivative along the trajectory. Zero is an assumption of
locally constant flow, not an automatic derivative of a throttle function.
Use smooth flow at ignition/cutoff and split scheduled smooth changes at
`propagate(..., breakpoints_s=...)`. Abrupt jumps need an impulse/state map.

## Model and integration

`model.mass_properties(configuration, propellant_kg, consumption_kgps,
consumption_rate_kgps2=...)` returns MassProperties including CG migration rate,
CG migration acceleration and mass second derivative. With VehicleMassAdapter,
configuration=None and the attached Stage objects define the vehicle.
`model.propellant_bounds(configuration)` returns lower/upper fuel limits.

`propagate(start_s,state,end_s,step_s,model,configuration,provider,...)` returns
times, accepted states and a stop reason. At a limiting tank it stops within
inventory_tolerance_kg (default 1e-6 kg); it does not invent a negative inventory
or automatic engine command. Continue with an appropriate next provider/phase.

## C outputs

`navigation_truth` returns Body-component CG position/velocity, relative-initial
attitude, Body angular rate, nongravitational force and Body gravity.
`acceleration_truth(state,loads,properties,sensor_station_body_m=...)` returns
geometric CG acceleration, fixed-sensor acceleration and ideal specific force.
Pass the SAME structural sensor station every time for a physically fixed IMU.
The default evaluates a fixed-body point coincident with CG at that instant;
it does not describe a sensor that mechanically follows the moving CG.

Body-coordinate acceleration is not the time derivative of its velocity
components without the transport correction. Test fixtures exercise that
distinction and compare against fixed-initial-basis conservation balances.
