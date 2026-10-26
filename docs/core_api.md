# Reusable A-block mathematics

This commit supplies mass-property calculations, coordinate mathematics,
load assembly and an RK4 rigid-body core independently of `VehicleData.py`.
Importing the four core modules does not instantiate a Falcon 9 or load vehicle
parameter files. NumPy is required; Matplotlib is only for separate examples.

| Module | Responsibility |
| --- | --- |
| `mass_properties.py` | Axisymmetric component CG/inertia and geometric derivatives; parallel-axis combination |
| `frames.py` | Hamilton scalar-first quaternion operations, Body/ENU rotations, WGS-84 ENU/ECEF geometry |
| `load_interface.py` | Sum application-point forces and moments about the supplied current CG |
| `dynamics.py` | CG state derivatives, RK4 integration and navigation truth outputs |

## Coordinate and interface contract

- Body +z points toward the nose; initially +x east, +y north, +z up.
- Reference frame is fixed launch-site ENU: east, north, up.
- Quaternions use [w, x, y, z], Hamilton products and Body-to-ENU mapping.
- Initial vertical quaternion is [0.5, −0.5, −0.5, −0.5].
- Structural stations use a fixed base datum in Body axes. Inertias/moments use current CG.
- `State.position_enu_m` tracks CG; `State.velocity_body_mps` is CG ground-relative
  velocity expressed in Body. Public navigation velocity is rotated to ENU.
- Earth rotation is omitted by the propagator. WGS-84 constants belong only to
  coordinate output and do not introduce Earth-rotation dynamics.

`rk4_step(t, state, dt, model, configuration, provider)` accepts a model object
with `mass_model_name` and a callable
`mass_properties(configuration, propellant_kg, consumption_kgps)` returning
`MassProperties`. No particular Vehicle class is imported by the core.

The load provider is called at every RK4 evaluation and returns `Loads`:
Body nongravitational force, Body moment about current CG, ENU gravity, and
nonnegative consumption rates for [stage 1, attached stage 2].
`assemble_loads` requires gravity explicitly; it has no assumed local gravity.
Do not include gravity in both the force list and the separate gravity vector.
Do not add an application-point moment twice.

## Current physical scope

Mass components and inertia calculations use an axisymmetric geometric model.
The liquid-column model requires settled, frozen-density cylinders; it does not
resolve slosh or fluid behavior during free fall. Dataclass freezing prevents
field reassignment but does not make contaienu NumPy arrays read-only.

The preliminary variable-inertia balance omits exhaust angular-momentum flux.
The core explicitly rejects positive consumption with `mass_model_name="dual_tank"`
until moving-CG kinematics and open-system momentum are treated consistently.
Fixed-mass diagnostic runs are supported. A successful numerical check is not
validation against an actual launch or landing.

## Checks

Run `python -m unittest discover -s tests -p test_core_math.py -v`.
Axial inertia is Ixx; transverse inertias are Iyy and Izz. CG station vectors are [axial station, 0, 0].

The fixtures are synthetic numerical inputs, not Falcon 9 parameters. Checks
cover analytical mass combination, force-arm sign, ENU orientation, free fall,
and constant-moment rotation. Vehicle-specific parameters and scenarios belong
in a separate commit.
