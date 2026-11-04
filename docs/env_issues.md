# Environment model: envelope, limits and the 175 km requirement

This note summarizes what `Environement.py` models today, why its ceiling is
120 km, what breaks if trajectories reach 175 km, and what the upgrade path
is. It is the reasoning behind the ceiling comment in `Environement.py`.

## What the model covers today

| Piece | Model | Envelope |
| --- | --- | --- |
| Atmosphere | USSA 1976 layers 0..84.852 km plus an explicit isothermal tail | 0..120 km, rejected outside |
| Gravity | Fixed direction along the initial up axis; constant magnitude, or inverse-square when an Earth radius is configured | flat-Earth convention |
| Wind | One fixed vector in the initial Body basis | scenario assumption |
| Ground | A single plane of constant altitude along the up axis | no terrain, no contact |

This covers launch, ascent, staging, coast, reentry and landing arcs that
stay below 120 km. Queries outside the envelope raise; the model never
extrapolates.

## Why the ceiling is 120 km and not higher

Three independent reasons:

1. **Atmosphere physics stops being valid above ~100 km.** The isothermal tail
   assumes constant temperature (187 K) and uniform air composition. Above
   ~100 km the thermosphere heats to 600..1000+ K and the constituents
   separate by mass, so a constant gas constant and heat-capacity ratio are
   no longer valid. A stretched exponential tail would produce numbers with
   no claim to correctness.
2. **No consumers.** Drag and ambient-pressure thrust are already negligible
   at 120 km. A taller table would refine quantities nobody uses.
3. **The frame model, not the atmosphere, is the real limit.** See below.

## The core issue: the simulation assumes a flat Earth

`gravity_body()` always returns gravity along one fixed direction (the
launch-pad up axis), and `altitude()` measures height against one infinite
flat plane. The whole simulation shares this flat-Earth convention: gravity
stays parallel everywhere, Earth rotation is absent, and there is a single
"up" direction for the entire run. This is an accepted approximation for the
current low-altitude, short-duration segments.

### Worked example: a 175 km booster arc

A booster coasting up to 175 km and back spends roughly
`2*sqrt(2h/g)` ≈ 6 minutes above the atmosphere, moving downrange at
1..2 km/s: several hundred kilometres.

| Quantity at 400 km downrange | Flat simulation | Round Earth (R = 6371 km) |
| --- | --- | --- |
| Gravity direction | parallel to pad up | tilted by 400/6371 ≈ 3.6° |
| Altitude reporting | height above a plane | surface curves away by `d²/2R` ≈ 12.5 km |

Integrated over six minutes, these errors move the predicted landing point by
tens of kilometres. Earth rotation (up to 465 m/s at the equator) is a third,
stacked error.

The important property: the simulation would **not crash** at 175 km. It would
run and produce a plausible-looking but wrong trajectory. Wrong-but-plausible
output is worse than an error, which is exactly why the module rejects
unmodeled altitudes instead of answering.

Analogy: the current frame is a flat map, adequate for the neighbourhood. The
175 km requirement is a transcontinental flight, where the flatness itself is
the dominant error. Better data painted onto the map does not fix it; the map
must become a globe.

## Consequences for the environment model

- The atmosphere is **not** the blocker. At 175 km the air is essentially
  vacuum; the column model only needs sourced tabulated values there.
- The blockers are **gravity direction** and the **altitude definition**.
  These live in the shared frame conventions (State bookkeeping,
  `frames.py`), not inside `Environement.py`. The environment can only
  describe the world in the coordinates the rest of the simulation speaks.

## Upgrade path

1. **Now:** keep the current model unchanged. Do not raise the ceiling before
   the frame upgrade; the 120 km rejection is the loud failure that prevents
   silent garbage. Preserve the existing split between frame-independent
   physics (`atmosphere(altitude)`, layer table) and frame-dependent adapters
   (`_height_along_up`, gravity direction, wind rotation).
2. **Near term (environment work, additive):** replace the analytic tail
   above the mesopause with sourced USSA 1976 tabulated values (density,
   temperature) as evidence JSON in `data/`, following the drag-curve
   convention. Document what temperature and sound speed mean above ~100 km,
   add checkpoint tests against the published table, and add a regression
   assertion that values below the mesopause are unchanged.
3. **Team decision (frame upgrade):** a geocentric or ECEF frame, at minimum a
   spherical Earth: gravity toward the planet centre, altitude above a curved
   surface, eventually Earth rotation. This means a new environment class
   alongside the current one with the same four method names
   (`atmosphere`, `altitude`, `gravity_body`, `wind_body`), plus an interface
   change the whole team must agree on (`gravity_body` and `altitude` need the
   state expressed in, or convertible to, an Earth-fixed frame). The current
   flat class stays as a verification baseline: the same low-altitude
   trajectory in both frames must agree.

## Open question for the team

Which 175 km is required?

- **Booster apogee on a suborbital arc:** needs spherical gravity and a
  curved-surface altitude, but no orbital mechanics. Moderate change.
- **Second stage to orbit at 175 km:** additionally needs ECEF, Earth
  rotation and orbital-capable dynamics. Substantially larger change.

The environment work (step 2, and the atmosphere half of step 3) is the same
in both cases, but the answer determines how much runway the rest of the
simulation needs.
