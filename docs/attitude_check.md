# Connected powered-attitude verification

This extends the same early-ascent integration with a prescribed, small lateral
force. Fuel flow, mass properties, steady propulsion and atmosphere remain the
same as the first vertical segment. The purpose is to exercise rotation and C's
sensor interface while mass and CG change.

```text
python main.py --config config/ascent_attitude_check.json
```

Default output: `results/attitude_check/`. The vertical baseline still runs with
`python main.py` and writes `results/powered_ascent/`; the runs do not overwrite
one another by default. Each exports the configuration, sources and truth CSV.

## Explicit verification inputs (class A)

| Input | Adopted value | Reason / limitation |
| --- | --- | --- |
| Force window | 10–14 s in the 30 s ascent segment | Rotation begins after liftoff; this is not a mission pitch program |
| Peak Body force | [0, 10,000, 0] N | About 0.13% of the initial 7.606 MN axial thrust, chosen to produce small but nonzero rotation; no hardware or wind measurement is claimed |
| Fixed application station | [1, 0, 0] m from structural base | Below current CG, to make the moment sign easy to check; not a measured aerodynamic center or TVC pivot |
| Force shape | Fᵧ(t) = 10,000 sin³(2πu) N, u = (t − 10)/4, for 0 < u < 1; zero otherwise | Smooth positive then negative force; force and its first two derivatives join continuously to zero |
| Primary virtual IMU | [0, 0, 0] m | Fixed base station retained from the first segment |
| Comparison virtual IMU | [40, 0, 0] m | Fixed second station, chosen to expose a 40 m rotational lever arm; not actual surveyed hardware |
| Angle-of-attack check limit | 3° | Selected verification guard; exceeding it raises an error. This is not evidence that the source's zero-AoA drag model is validated up to 3° |

The sinusoidal doublet is an analytic verification load: it is explicitly
separate from aerodynamic forces and engine gimbal commands. Its Body-component
time integral is zero. Rotation, changing mass/inertia, exhaust angular momentum
and small engine ambient-pressure differences mean that zero final angular rate
or zero net inertial impulse is **not** imposed or guaranteed.

No new vehicle physical parameter is inferred from this exercise. The literature
inputs and sources remain those in [the early-ascent guide](first_ascent_segment.md).
Normal-force/moment aerodynamics, TVC and feedback control remain future models.
In particular, drag still acts through CG using the source zero-AoA Cᴅ curve;
the small tilted trajectory is a verification surrogate, not an aerodynamic
prediction validated against a real flight.

## Expected signs and sensor behavior

The arm from CG to [1,0,0] has negative x. During the positive-y lobe, the
cross-product M = r × F therefore has **negative Body z**. Starting from rest,
angular acceleration and rotation initially point along negative z, moving the
nose toward negative initial y. During the second lobe the applied moment
reverses. The solver computes the resulting attitude; it does not set a desired
angle directly. Body +x remains noseward throughout.

Gravity retains its physical downward direction while its Body components change
with attitude. The fixed initial Body basis is only internal orientation
bookkeeping, not ENU/ECEF. Shared vector outputs remain current Body components.

Both accelerometers have the same gravity subtraction. If Δr is the vector from
the primary fixed sensor to the comparison sensor, their ideal specific-force
difference is:

Δf = α × Δr + ω × (ω × Δr).

Here Δr = [40,0,0] m, α is Body angular acceleration in rad/s², and ω is Body
angular rate in rad/s. The first term is tangential acceleration and the second
is centripetal acceleration. Common CG migration and gravity terms cancel in
this sensor difference. This is why two sensors on the same rotating vehicle can
report different accelerations.

## Additional truth columns

| Column family | Meaning |
| --- | --- |
| prescribed_force_body_n_* | Injected verification force in current Body, N |
| force_nongrav_body_n_* | Total thrust, aerodynamic and prescribed force, excluding gravity, N |
| moment_cg_body_nm_* | Total external moment about the instantaneous CG, N·m |
| angular_acceleration_body_radps2_* | Angular acceleration truth in current Body, rad/s² |
| tilt_from_initial_axis_rad | Unsigned angle between current nose axis and the initial nose axis, rad; not an Euler angle |
| comparison_accelerometer_ideal_body_mps2_* | Ideal specific force at the second fixed station, m/s² |

The existing gyro, quaternion, primary accelerometer and gravity columns now
contain nontrivial attitude information. C can use these for sensor simulation
and debugging; an estimator still consumes measurements with its selected error
models, not the full truth columns as measurements.

## Verification and remaining work

Tests independently check r × F signs for positive y and positive z forces,
zero-rate angular acceleration, reversed-force symmetry, gravity sign, the two
sensor lever-arm equation, pulse endpoints and angle-of-attack guard. Time-step
refinement checks the expected improvement of the integrated attitude. Scheduled
pulse start, midpoint and end are integration breakpoints.

This finishes a connected A rotation/sensor verification case. It does not finish
all A work. Next: agree sensor mounting stations with C; replace verification
forcing with justified B/control loads as available; then connect phase timing,
smooth engine shutdown/restart and staging state maps with conservation checks.
