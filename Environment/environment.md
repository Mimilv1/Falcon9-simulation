# The Environment model

`Environement.py` is the world the Falcon 9 flies through: everything the
vehicle interacts with but does not carry. It is a pure service module —
dynamics calls it, it never calls back — and it covers exactly six services:

| Piece | Call | Output |
| --- | --- | --- |
| Atmosphere | `Environment.atmosphere(altitude_m)` | `Atmosphere`: T, p, ρ, speed of sound |
| Gravity | `env.gravity_body(state)` | acceleration vector, Body components, m/s² |
| Gravity gradient | `env.gravity_gradient_ecef(position_ecef_m)` | 3×3 Jacobian for the navigation filter |
| Frame acceleration | `env.frame_acceleration_body(state)` (also `_ecef`, gradients) | Coriolis + centrifugal acceleration, Body components, m/s² |
| Altitude over ground | `env.altitude(state, properties, station=None)` | scalar metres above the spherical ground |
| Wind | `env.wind_body(state)` | velocity vector, Body components, m/s |

## Thoughts baked into the model

**One frame, project-wide.** The ECEF origin is the Earth centre (0,0,0),
the launch site sits at `(0,0,earth_radius)` on +z. Earth rotation enters
translation: `frame_acceleration_ecef` / `frame_acceleration_body` evaluate
the Coriolis and centrifugal accelerations of the rotating frame from a
configurable `earth_rotation_rate_radps` and `launch_latitude_deg`
(defaults: the sidereal rate 7.2921150e-5 rad/s at 28.6°N; rate 0 disables
them). Dynamics adds this term once next to gravity, and the navigation
filter uses `frame_acceleration_gradients_ecef` in its prediction Jacobian.
The angular equations and gyro Earth-rate treatment retain the team
approximation; no complete inertial-rate coupling is claimed. All shared
vectors the module returns are *current Body* components; the quaternion
(`[w,x,y,z]`, Hamilton, Body → ECEF) in the state is the only attitude
input. At launch, Body +x points along ECEF +z (nose up).

**Gravity is a reduced Newtonian point mass.** F = G·M·m/r² is divided by
the vehicle mass once and for all: `a = μ/d²` toward the planet centre. μ
is built from the scenario configuration as `g_surface·R²`, so the gravity
at the launch site is exactly the configured value by construction — and
swapping in `g=3.7, R=10⁶` gives a consistent Mars-like planet for tests.
Gravity never becomes a force in newtons anywhere; `dynamics.py` adds it
as acceleration next to thrust/mass. The 3×3 gravity-gradient map
(`gravity_gradient_ecef`) exists for the navigation filter's prediction
Jacobian; it is not used for torques.

**Atmosphere: evidence up high, honest placeholder above.** Layers
0..84.852 km reproduce the USSA 1976 tables (validated against published
checkpoint values). Above the mesopause an isothermal tail decays
exponentially toward numerical vacuum (ρ < 1e-150 at 2000 km) — usable
drag numbers up there, deliberately no thermospheric physics. The 2000 km
envelope is the mission scope; anything outside raises a `ValueError`.
Never extrapolated. Missing regimes are rejected per project convention.

**Wind is a rotated surface reference vector with a stylized altitude
profile.** The scenario supplies a surface reference wind `wind_ecef_mps`,
optionally rotated about the launch-site vertical by `wind_rotation_rad`
so it spreads over the horizontal plane (crossflow, not one lone
coordinate, angles the vehicle). `wind_body` scales this vector by the
declared `WIND_PROFILE` fraction of the CG altitude: 1.0 at the surface,
a 2.0 jet-stream core at 11 km, decaying through the stratosphere and
exactly zero at the mesopause — there is no wind above the physically
modeled atmosphere. No gusts, no time dependence, no direction shear. It
enters physics only as the air's velocity subtracted from the vehicle
velocity in the aerodynamics relative-speed term; above the atmosphere it
is zero by construction, so its presence is harmless.

**Ground is a bare sphere**, radius `earth_radius`, used only as the
altitude datum; `ground_altitude_m` shifts its level (e.g. launch site
100 m above the ground sphere). No terrain, no contact models.

## Code basics

- One frozen dataclass (`Atmosphere`), the module helpers (`_layer_pressure`,
  `_wind_scale`, `_skew`), one class (`Environment`): scenario configuration
  plus pure evaluators. ~260 lines total.
- **Validation at construction**: an invalid scenario fails before
  flight, not mid-integration. Non-finite gravity, radius, ground offset,
  wind or wind rotation, non-positive radius, an Earth rotation rate below
  zero (zero disables it) and a launch latitude beyond ±90° all raise
  `ValueError`; parameters are keyword-only.
- **Evaluators are pure and stateless**: same state in, same values out,
  no memory between calls, no mutation of the caller's arrays.
- `altitude` optionally takes a station offset from the CG (used e.g.
  for engine-exit height); with `station=None` the properties argument is
  unused.

## Tests

`Environment/test_environment.py` — 41 dedicated tests, run with:

```bash
python3 -m unittest discover -s Environment
```

(in a shell with the project venv activated; without activation use
`.venv/bin/python -m unittest discover -s Environment`)

Covered: published USSA checkpoint values (0..80 km), layer-boundary
continuity, monotone pressure/density decay to the ceiling, hydrostatic
balance (dp/dh = −g·ρ) inside layers, exact tail scale height, envelope
rejection both sides of the boundary; gravity magnitude on pad, inverse-
square falloff, attitude round-trip via ECEF, downrange tilt, planet-centre
rejection; wind mapping at the pad, altitude decay and profile knots with
the jet-stream peak, calm above the modeled atmosphere, wind-rotation
spreading and its crossflow effect on the vehicle, attitude rotation with
preserved magnitude, constructor copy semantics; altitude zero point,
curved surface, attitude invariance, ground-level shift, station offsets;
constructor validation for every parameter; and simulator-plausibility
checks (max-Q dynamic pressure, boostback arc regime). The frame-
acceleration services are exercised through the ascent integration and
navigation-filter tests in `tests/`, not duplicated here. One real bug was
found by this suite and fixed in the module:
`surface_gravity_mps2=nan` used to be silently accepted.

The shared test suite in `tests/` owns the ascent-level integration and
cross-module tests (115 tests, all green at the time of writing).

## Limits, by design

Coriolis and centrifugal translation effects are included, but full
Earth-rate coupling in the angular equations and the gyro model is not;
J2 oblateness, wind gusts and time dependence, terrain and contact forces
are all out of scope. The launch-aligned (non-Greenwich) ECEF frame and
the approximate angular Earth-rate treatment are the binding fidelity
limits; the others refine accuracy inside them. Upgrade path: the frame
services are pure evaluators, so a fuller rotation model or a J2 gravity
gradient slots into `gravity_ecef` / `frame_acceleration_ecef` without
touching callers.
