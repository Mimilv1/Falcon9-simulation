# Documentation flags

Status board for places where the prose documentation lags the code or
`data/`, and for pending changes that will make it lag. Anyone (human or
agent) who cites a documented statement that is flagged here must verify it
against the code and data first; the code and the JSON records are the
truth, the prose describes it.

When you find a new doc/data discrepancy: verify the code or data side
yourself, do not propagate the stale statement, and append a flag here.

## Open flags

docs/image.png include?

- **FL-1 — Configuration count in `docs/aero_model.md`** (opened
  2026-11-09). The doc states "four declared vehicle configurations" and
  lists fin/leg intermediate states under "Open items".
  `data/aero_configurations.json` now contains **five** configurations:
  `booster_return` was added (fins deployed, legs stowed; active during
  entry burn and aerodynamic descent). Its axial curve equals
  `booster_landing`'s; the difference is the deployed-leg crossflow
  segment. Action: update `aero_model.md` (config table, "Open items"
  section, configuration paragraphs).

- **FL-2 — Actuator layer documentation** (opened 2026-11-09).
  `Actuators.py` is an empty placeholder today; the actuator layer is
  planned and will be implemented before the report is submitted. When it
  lands: `docs/aero_model.md` ("Open items: actuator layer") and
  `docs/core_api.md` ("not implemented" statements) become stale and must
  be reconciled with the implementation; the mission-loop description and
  all test counts must be re-verified; the nominal full simulation and the
  report figures behave on it (`results/`-derived figures may need
  regeneration). The report (architecture, aerodynamics and
  return/actuation chapters) is written to cover it and must be updated to
  the final implementation.

- **FL-3 — `second_stack` aerodynamic configuration is declared but not
  wired into the active stage-2 continuation** (opened 2026-11-09,
  verified by code search: only `tests/test_aerodynamics.py` references
  `second_stack`). The upper-stage run (`upper_stage.py`) uses its own
  simple axial-drag estimate (declared $C_D = 0.3$ on the
  $10.52\,\mathrm{m^2}$ stage cross-section), not the
  `Aerodynamics`/`second_stack` model. `docs/aero_model.md` does not
  mention this. Action: state the selection explicitly where aero
  configurations and upper-stage propagation are documented; optionally
  switch the continuation onto the declared configuration.

- **FL-4 — `data/crs10_mission.json` `integration.attitude_source` is
  stale** (opened 2026-11-10, verified against `mission.py`/`kalman.py`).
  The record states "truth; current C filter does not estimate attitude",
  but `mission.py` now defaults `use_estimated_attitude` to true and the
  filter runs the 15/24-state MEKF with `GyroAttitude` and star-tracker
  aiding; `summary.json` overwrites the key with the runtime value
  (`mekf_gyro_gnss_lidar` in the current nominal run). The json string must
  not be quoted in the report; cite `mission.py` and the run record instead
  (the report Chapter 7 already does). Action: update the json record.

- **FL-6 — Navigation covariance and gate performance** (open,
  2026-11-10): The seed-1 full run rejects all 381 valid star fixes.
  Position-error time fractions within reported three-sigma envelopes are
  about 0.500/0.344/0.752. Chapters 7 and 9 now distinguish accuracy from
  consistency; runtime covariance/gating still requires investigation.
  The GNSS measurement covariance includes white-noise plus stationary
  Gauss-Markov marginal variances, but does not represent temporal
  correlation; the final Chapter 7 now states that distinction explicitly.

- **FL-7 — Jacobian certification** (open, 2026-11-10):
  `tests/test_gyrostartracker.py` prints finite-difference diagnostics without
  tolerance assertions. Chapter 9 now states that limitation explicitly.
  Quantified tolerance/step-size acceptance is still needed for a stronger
  verification claim.

- **FL-8 — Ascent heading field** (open, 2026-11-10):
  `ControlledAscent` stores `heading_ecef` from the configured angle but
  constructs the target using projected `EAST_ECEF`. The JSON's 45-degree
  heading explanation does not describe the active pointing direction.
  Chapter 8 documents the implemented eastward command without changing
  the runtime or reinterpreting the trajectory as an identified ISS plane.

- **FL-9 — Vehicle metadata source and scope strings** (open,
  2026-11-10): `data/vehicle_parameters.json` retains obsolete body-only,
  non-rotating and booster-only scope descriptions. Its cargo row points
  `source_ids` to S10 (the historical engine guide) instead of the named
  NASA CRS-10 overview. These strings do not override runtime behavior.
  The final report appendix resolves the cargo citation to the overview
  and uses the implemented rotating-frame, two-child and MEKF paths.
  Metadata cleanup remains separate from the unchanged R1 runtime inputs.

## Resolved flags

- **FL-0 — `docs/env_issues.md` references and Earth-rotation wording**
  (resolved 2026-11-09). `Environment/environment.md` previously claimed
  Earth rotation was not modeled and cited a nonexistent
  `docs/env_issues.md`; `docs/parameter_evidence.md` also cited it. Fixed
  in the same session: rotation is documented as supplied Coriolis and
  centrifugal accelerations, the wind profile is documented, and dead
  references were removed. Kept for history; no action needed.

- **FL-5 — Report frame and notation review** (resolved in the local
  Chapters 1-9 report, 2026-11-10): Corrected rotating-frame wording,
  lateral-flow norm, per-axis damping, active versus declared upper-stage
  drag, sensor period relation, A classification, gyro random-walk terms,
  nominal/error-state dimensions and strict monotonic-drift wording.
  Navigation performance moved to Chapter 9; statistical consistency is
  no longer asserted from mean NIS alone. Runtime limitations remain open
  under FL-6/7/8.
