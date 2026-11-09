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

## Resolved flags

- **FL-0 — `docs/env_issues.md` references and Earth-rotation wording**
  (resolved 2026-11-09). `Environment/environment.md` previously claimed
  Earth rotation was not modeled and cited a nonexistent
  `docs/env_issues.md`; `docs/parameter_evidence.md` also cited it. Fixed
  in the same session: rotation is documented as supplied Coriolis and
  centrifugal accelerations, the wind profile is documented, and dead
  references were removed. Kept for history; no action needed.
