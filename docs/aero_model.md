# Aerodynamics model: configurations, control and evidence

This note documents what `Aerodynamics.py` computes, which coefficients are
sourced and which are declared estimates, how moments and grid-fin control
enter the shared Body frame, and where the model's honest limits sit. It is
the reasoning behind the module docstring and the records in
`data/aero_configurations.json`.

Scope: one `Aerodynamics` model instance, four declared vehicle
configurations, optional grid-fin control. Not in scope: propulsion gimbal
(TVC), actuator rate/travel dynamics, a control law, and any mission-timing
decision about which configuration is active when — those belong to the
mission layer and are listed under "Open items".

## What the model computes

All forces and moments are current Body components, +x toward the nose,
stations in the same structural datum as the phase's mass properties. The
input is the relative velocity `v_rel = v_body − CG migration − wind_body`
(the convention shared with `DragOnlyAerodynamics`, which remains available
and byte-identical on zero-angle flight).

The total flow angle is `a = atan2(|v_rel_yz|, v_rel_x)` in **0..pi**: 0 is
nose-first, pi is engine-first (retro), and every direction in the y-z plane
is handled, not just a pitch plane.

| Term | Formula | Acts at |
| --- | --- | --- |
| Axial force | `-qbar*S_ref*CA0(M)*cos(a)*|cos(a)|` along +x | on the axis: no moment about an on-axis CG |
| Potential normal force | `qbar*S_ref*CNa*sin(a)*|cos(a)|` along `-(v_rel_yz direction)` | windward-end station (see moments) |
| Crossflow normal force | `qbar*CDc*S_side*sin(a)^2` per side segment, same direction | segment side-area centroid |
| Grid fin control (optional) | `qbar*S_fin*Cl_delta*delta_i*cos(a)*|cos(a)|` along the fin tangential direction | the fin's 3-D mount station |
| Rotary damping | `c * qbar*S_ref*L_ref^2/(2V) * omega`, pitch=yaw by symmetry, roll separate | moment about the CG |

Physics notes worth knowing:

- The `cos*|cos|` factor appears twice on purpose. Axial force scales with the
  axial dynamic-pressure share and flips direction for retro flight; fin
  control effectiveness scales with the through-flow share the same way, so
  the same deflection command reverses its sign in engine-first flight and
  vanishes at 90 deg total angle. The alternative — direction-independent
  effectiveness — would silently steer the wrong way during descent, the one
  phase grid fins exist for.
- The normal-force direction is always along `-(lateral relative-velocity
  direction)`; with the station moment this pair is self-consistent and
  reproduces the textbook stability result: a windward-end pressure center
  lying toward the nose is **destabilizing** in nose-first flight, deployed
  hardware near the trailing end is stabilizing in tail-first descent. The
  model does not stabilize the vehicle for you; Falcon-class vehicles are
  aerodynamically unstable and need active control.

## Configurations

| Name | Datum (station 0) | S_ref / L_ref | Axial curve | Grid fins | Phase |
| --- | --- | --- | --- | --- | --- |
| `full_stack` | stack base (S1 engine bay) | 21.08 / 62.4 m | Pezent 2024 readout, extended to M 10 | no | launch to separation |
| `booster_alone` | booster base (unchanged after `separate(1)`) | 10.52 / 42.6 m | declared estimate | no | post-separation coast, boostback |
| `booster_landing` | same as booster_alone | 10.52 / 42.6 m | declared estimate (+fin/leg increment) | yes | entry, landing burn (fins/legs present) |
| `second_stack` | stage-2 aft skirt base (renumbered by `separate(1)`) | 10.52 / 19.8 m | declared estimate | no | S2+Dragon ascent |

The provider that selects the active configuration owns mission timing;
selection is explicit (`evaluate(..., configuration_name)`) and an unknown
name raises. Default in `ascent.py` is `full_stack`. Which booster
configuration applies during a real descent (fins deploy between phases) is a
mission-layer switch, not an aero-model decision.

## Coefficients and evidence classes

Following the project convention, every entry in
`data/aero_configurations.json` carries `evidence`/`uncertainty`/`notes`
records:

| Coefficient | Class | Where the number comes from |
| --- | --- | --- |
| `full_stack` CA0, M 0..7 | sourced (`data/drag_curve_pezent.json`) | Pezent 2024 Fig. 2.13(b) raster readout, RasAeroII, 0 deg AoA, 21.08 m2 reference. Pixel readout, transonic rise worst. |
| `full_stack` CA0, M 7..10 | declared estimate (`drag_curve_pezent_extended.json`) | slow hypersonic decay; physical impact <0.5% of thrust above 55 km |
| booster/landing/second-stack CA0 | declared estimate | blunt-cylinder and sleeker-body engineering estimates, +/-50% |
| CNa = 2.0/rad | declared estimate | slender-body theory, referenced to the body cross-section |
| CDc = 1.1 (1.0 fins, 1.2 legs) | declared estimate | subcritical cylinder crossflow values |
| pitch/roll damping, cp stations | declared estimate | order-of-magnitude, static-margin geometric choices |

Only the full-stack axial curve is sourced; every other number is a
documented estimate with no flight evidence behind it. They are good enough
to close the simulation loop and produce plausible dynamics, **not** good
enough to predict quantitative trajectories. Replacing any of them is a data
change in the JSON, no code.

## Moments and the wildward-end convention

- Every station moment is computed once as `(station − CG) × F`, the same
  convention as `load_interface.arm_from_station`; the model returns the
  resultant as a free moment about the CG and the force through the CG, so
  `assemble_loads` cannot double-count.
- The potential normal force acts at the configuration's windward-end
  station: the declared station nose-first, **mirrored about the CG x** for
  retro flow, clamped to the body. The mirror is a declared approximation
  (equivalent windward-end geometry), validated by tests that pin both
  directions.
- Side segments are incremental: grid fins and legs add projected width on
  top of the body segment, and each carries its own centroid and crossflow
  drag. In tail-first descent the fin segments sit near the trailing end and
  therefore stabilize — this falls out of the geometry, it is not coded.

## Grid fin control

`evaluate(..., grid_fin_deflections_rad)` accepts one deflection per declared
fin (radians, instant commands):

- Fin i sits at `(station_m, R*cos(phi_i), R*sin(phi_i))`, azimuths about
  Body +x from +y toward +z; `booster_landing` declares the four intercardinal
  mounts.
- A positive deflection produces a force along `(0, -sin(phi), cos(phi))`,
  scaled by the through-flow share (see above), applied at the 3-D mount
  station. The axial lever carries pitch/yaw control; the radial offset
  carries roll. Because the station is three-dimensional, the moment
  combinatoric is automatic:

| Command pattern | Net force | Net moment |
| --- | --- | --- |
| equal deflection on all fins | zero | pure roll |
| opposite pair, differential deflection | side force | pure pitch, zero roll |
| single fin | side force | roll + pitch, exactly once each |

- Command contract: one value per declared fin, finite, within the declared
  travel (20 deg); violations raise. Configurations without declared fins
  reject nonzero commands. No deflection commands means byte-identical to
  the passive model.
- Authority check (declared values, Mach 3 descent): about 1.4 MNm of pitch
  and 385 kNm of roll at 15 km altitude at full deflection, an order of
  magnitude less above 30 km — consistent with grid fins mattering mainly
  below entry-burn altitudes.

Not modeled, documented in the module and JSON: Mach/alpha dependence of fin
effectiveness, interference with the passive crossflow segments, deflection
axial-drag increment, lattice state, and all rate/travel dynamics (they
belong in an actuator layer upstream, which currently does not exist in the
project).

## Domain and rejection behavior

The model follows the project rule — reject, never extrapolate — in four
places:

1. Mach outside the adopted axial curve (per configuration, all end at M 10).
2. Unknown configuration name.
3. Grid fin command shape or travel violations.
4. Construction-time station/range violations (configuration rejects, not a
   mid-flight break).

The caution runs the other way too: inside the domain, declared-estimate
coefficients do not crash — they produce plausible numbers with real
uncertainty (the booster axial curve alone is good to roughly +/-50%). The
tests pin conventions and invariants, not coefficient truth.

## Verification map

`tests/test_aerodynamics.py` pins, class by class:

| Test class | What is pinned |
| --- | --- |
| AxialEquivalence | zero-angle flight equals `DragOnlyAerodynamics` exactly |
| LateralForceAndMoment | normal-force direction, 90-deg crossflow CN value, segment fin increments |
| RetroFlow | engine-first axial opposition, windward-end mirror convention (both directions) |
| RotaryDamping | damping formula term by term, opposing signs |
| DomainAndRest | domain rejections, rest case zeros |
| ConfigurationData | JSON loads all four phases, sourced knot intact, fin/leg segments |
| HypersonicExtension | all configurations evaluate across M 7..10, reject at M 10.5, ascent propagates past Mach 7 |
| PhysicalInvariants | force extracts energy at every angle, damping dissipates, wind enters only through relative velocity, 90-deg continuity, qbar/Mach scaling, term-by-term formulas |
| GridFinControl | single-fin formula at all azimuths, roll/pitch combinatorics, retro reversal, axial-share scaling, additivity, boundary deflections, rest case, roll moment reaching `dynamics.accelerations` with the correct angular acceleration |
| ProviderIntegration | provider selection, unknown-name rejection, two-second vertical propagation through the A core |

## Open items

- **Control law** — nothing commands the fins or a future TVC model; the aero
  model accepts and validates commands only.
- **Provider wiring** — `SteadyPoweredAscent` does not forward fin commands;
  a descent-phase provider must pass them per step.
- **Actuator layer** — rate/travel/deadband dynamics are a declared upstream
  responsibility with no current implementation in the project.
- **Evidence upgrade** — every declared estimate in the JSON wants replacing
  with sourced data when available.
- **Configuration switch modeling** — fin deployment and leg deployment are
  binary today (`booster_alone` vs `booster_landing`); intermediate states
  (fins deployed, legs stowed) would be additional declared configurations.
