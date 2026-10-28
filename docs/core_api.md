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

`provider(time_s, trial_state) -> Loads` is evaluated at each RK4 stage.
The temporary `SteadyPoweredAscent` provider in `ascent.py` connects the retained
thrust, drag and environment models. It does not select a flight scenario.

| Loads field | Definition |
| --- | --- |
| force_body_n | Net nongravitational force, N; thrust already includes exhaust momentum/pressure reaction |
| moment_body_nm | External moment about current CG, N·m |
| gravity_body_mps2 | Gravity in current Body, added once |
| consumption_kgps | Nonnegative [stage1,stage2] consumption, kg/s |
| exhaust_ports | Stage index, fixed exit station, flow, analytic flow derivative and exit gyration |

`assemble_loads` sums station forces and r × F once. Do not add the same force
moment again. Gravity direction must follow attitude in Body components.
Port flows must match consumption. This core requires axial CG migration,
balanced transverse exhaust centroids and fixed exit stations; gimbaling exit
boundary motion and internal swirl are not implemented.

## A integration

`vehicle.mass_model()` evaluates trial fuel inventories without changing the
accepted Stage objects. Use configuration=None. `vehicle.to_state(...)` and
`vehicle.apply_state(...)` map structural-base and geometric-CG states, including
CG migration. Supply matching consumption and derivative values to both maps;
only synchronize accepted states.

`propagate(..., breakpoints_s=...)` integrates smooth segments and stops at a
limiting tank. It does not choose cutoff, separation or restart conditions.
Supply consistent flow and dq/dt. Abrupt flow jumps need a separate impulse map.

`SmoothShutdown(start_s,duration_s)` is an optional analytic envelope, not measured
engine behavior; the caller must justify its parameters. No event times are set.
`separate_coasting_vehicle(...)` requires zero flow/derivative and off engine
commands. It returns two child states and a conservation audit without modifying
the caller. The no-impulse/no-discard assumption must match the selected model.

## C truth

`navigation_truth(state,loads)` provides CG position/velocity, relative-initial
attitude, angular rate, nongravitational force and gravity in Body components.
`acceleration_truth(...,sensor_station_body_m=...)` provides ideal specific force
at a fixed structural sensor station. Use the same fixed station every time;
CG coordinate acceleration is not a physical IMU measurement.

Sensor errors, sampling and navigation estimation belong to C. Mission-specific
initial states, event conditions and result export remain team integration work.
