# Report Guide (for AI agents)

Operating manual for writing and editing the technical report in this directory.
Read this file fully before touching any `.tex` file. The workflow and conventions
here are binding; skip them and the result will be rejected at the final check.

## Deliverable and hard requirements (read first)

1. **Scope**: technical report documenting the Falcon 9 / CRS-10 simulator in the
   repository root. It documents the system — models, implementation, results —
   not the process of writing the report.
2. **Length: at least 15 pages of body in the compiled `main.pdf`, figures
   included. There is no upper limit.** The gate covers the title page and
   table of contents through the end of the Conclusions chapter; the References
   list and the appendix are outside the gate and do not count toward it.
   Re-check the page gate (below) after every structural change; never
   discover the page count at the end.
3. **Language**: precise, technical English. Present tense, the system is the
   subject (not "we").
4. **Grounding**: every number, equation, parameter, and result must come from this
   repository (see source map below). Never invent values, citations, or results.
   Where evidence does not exist, state the assumption explicitly as an assumption.
5. **Evaluation focus**: the finished report is graded for technical completeness
   and precision. A top grade requires: governing equations with consistent
   symbols, quantified results (values, units, error percentages), explicit
   assumptions and limitations per model, validation against reference data, and
   every figure/table/equation referenced from the text. Placeholders, unexplained
   numbers, vague claims, and orphaned floats guarantee a bad grade.

## Workflow (mandatory, in this order)

1. **Research before writing.** For each chapter, read its sources (table below)
   and extract: the model equations, parameter names/values/units, where each
   parameter comes from, known simplifications, and test or validation results.
   Never start a chapter from memory or plausibility alone. Consult
   `../docs/flags.md` first: it lists known places where the prose docs lag
   the code or `data/*.json`. Flagged statements must be verified against
   the code or data record before the report cites them; when you find a new
   discrepancy, do not propagate it, and append it to `../docs/flags.md`.
2. **Outline first.** Fix the chapter list (binding outline below) and confirm
   the planned total clears the minimum of 15 body pages (there is no upper
   limit). Register all `\input` lines before writing content.
3. **Write one chapter at a time.** Finish a chapter including its equations,
   floats, and citations, then run the build and fix warnings before moving on.
4. **Generate figures from real data**, not sketches (figure rules below).
5. **Iterate with the verification loop** (below) after every chapter: build,
   page count, warning check, placeholder check.
6. **Final pass**: notation review across all chapters (same symbol = same
   meaning everywhere), then run the full definition-of-done checklist.

## Layout

- `main.tex` — master document: metadata, package loading, ordered `\input` list.
- `content/` — one file per chapter of the binding outline, starting with
  `\section{...}`; the appendix file starts with `\appendix`.
- `images/` — all figures (PDF preferred).
- `content/acronyms.tex` — all acronyms via `\DeclareAcronym{KEY}{...}`.
- `content/sample_content.tex` — syntax reference for tables, figures, headings.
- `PROPOSED_OUTLINE.md` — the adopted structure proposal this guide's binding
  outline implements; keep the two consistent when either changes.
- `tum/`, `tumuser.sty`, `latexmkrc` — class and build support. Never edit these.

## How to add a chapter

1. Create `content/<name>.tex`, starting with `\section{Chapter Title}`.
2. Register it in `main.tex` below the `% +++++ Main part` marker, in PDF order.
   Chapters are registered with `\include`, not `\input`: `\include` starts
   every chapter on its own page (automatic `\clearpage`) and writes a
   per-chapter `.aux` file:

   ```latex
   \include{content/introduction}
   %\include{content/next_section}
   ```

3. Build with `latexmk -pdf main.tex` from this directory.

## Source map: where the truth lives

Read these before writing the matching chapter. Numbers in the report must be
traceable to one of these (or to the code they describe).

| Report chapter | Sources |
| --- | --- |
| 1–3 Mission, vehicle, dynamics (incl. separation) | `docs/crs10_mission.md`, `data/crs10_mission.json`, `mission.py`, `Events.py`; `README.md`, `docs/parameter_evidence.md`, `data/vehicle_parameters.json`, `VehicleData.py`, `mass_properties.py`, `mass_flow.py`; `docs/powered_dynamics.md`, `docs/core_api.md`, `dynamics.py`, `frames.py`, `propagation.py`, `load_interface.py`, `staging.py` |
| 4 Environment (atmosphere, gravity, wind) | `Environment/environment.md`, `Environement.py`, `Environment/test_environment.py` |
| 5 Aerodynamics | `docs/aero_model.md`, `data/aero_configurations.json`, `data/drag_curve_pezent*.json`, `Aerodynamics.py`, `tests/test_aerodynamics.py` |
| 6 Sensor models and measurement generation | `sensors.py`, `mission.py` (GyroAttitude), `README.md` ("Coordinates" section) |
| 7 Navigation and state estimation | `kalman.py`, `README.md` ("Coordinates" section) |
| 8 Mission execution, control, actuation | `propulsion.py`, `engine_schedule.py`, `mass_flow.py`; `staging.py`, `ascent.py`, `recovery.py`, `mission.py`; `Actuators.py` once implemented; `tests/test_staging.py`, `tests/test_recovery.py` |
| 9–10 Validation, results, assessment | `data/crs10_ascent_reference.*`, `data/ascent_model_audit.json`, `docs/crs10_mission.md`, `tests/`, run outputs under `results/` |

The coordinate conventions from the root `README.md` (Body +x to the nose;
ECEF is launch-aligned, not Greenwich) are binding: do not silently rename them.

## Binding outline and page guidance

`PROPOSED_OUTLINE.md` is the adopted report structure; the table below is the
binding version (files, page guidance, responsibilities). The report follows the
simulator's data flow: define the vehicle and shared conventions, explain the
physical and measurement models, describe mission execution, then evaluate
results and their engineering implications.

| # | Chapter | File | Pages | Lead |
| --- | --- | --- | --- | --- |
| 1 | Introduction and Project Scope | `content/introduction.tex` | 0.75–1 | Kwanwoo |
| 2 | Simulation Architecture and Shared Conventions | `content/architecture.tex` | 1.5–2 | Kwanwoo |
| 3 | Vehicle Configuration, Dynamics, and Mission Definition | `content/vehicle_dynamics.tex` | 2.5–3 | Kwanwoo |
| 4 | Environment Models | `content/environment_models.tex` | 1.25 | Philipp |
| 5 | Aerodynamics | `content/aerodynamics.tex` | 1.5–1.75 | Philipp |
| 6 | Sensor Models and Measurement Generation | `content/sensor_models.tex` | 1.25 | Emilien |
| 7 | Navigation Architecture and State Estimation | `content/navigation.tex` | 1.25–1.5 | Emilien |
| 8 | Mission Execution and Numerical Implementation | `content/mission_execution.tex` | 2 | Kwanwoo |
| 9 | Verification, Validation, and Simulation Results | `content/results.tex` | 2–2.5 | Joint (assembled by Kwanwoo) |
| 10 | Engineering Assessment and Implications for Future GNC | `content/assessment.tex` | 1–1.25 | Joint |
| 11 | Conclusions | `content/conclusions.tex` | 0.5–0.75 | Joint |
| — | References (unnumbered) | `content/references.tex` | outside gate | All |
| A | Appendix: Reproducibility and Supporting Data | `content/appendix.tex` | outside gate | Shared |

Sum of midpoints ≈ 16 body pages, which clears the minimum-15 gate with
margin. The Per-chapter pages are rough planning guidance for distribution
only, not budgets: there is no upper limit, and no chapter is trimmed to fit
a per-chapter number. Roughly 450 words of body text per page in this
layout; floats and equations take the rest. Treat `texcount` (below) as the
word-total gauge and the compiled body page count as the criterion. The appendix carries the parameter/assumption/source register,
variable definitions and units, run configurations, long code excerpts and
additional results — anything pulled out of the body goes there, with a cross-reference.

### Chapter contents (per `PROPOSED_OUTLINE.md`)

Concrete content of each chapter; titles, files, page budgets and leads are
in the binding outline table above.

1.  **Introduction and Project Scope** — Project 2 objectives, why CRS-10 was
    selected, modeled flight scope, and validation objectives.
2.  **Simulation Architecture and Shared Conventions** — overall data flow
    (truth versus measurements versus estimates), coordinate frames, the
    quaternion convention, reference points, units, module interfaces, the
    mission loop and reproducibility.
3.  **Vehicle Configuration, Dynamics, and Mission Definition** — the vehicle
    and separation material stated once, in one place:
    1.  **Vehicle Configuration and Evidence Classification** — CRS-10 stage/Dragon
        configuration, geometry and mass budget, K/E/A definitions, provenance.
    2.  **Propulsion and Propellant Depletion** — thrust, ambient-pressure
        correction, specific impulse, mass flow, loaded vs consumed O/F, inventory.
    3.  **Mass, Center of Gravity, and Inertia** — dry-component and liquid models,
        component summation, property derivatives during fuel consumption.
    4.  **Equations of Motion and Load Interfaces** — translational and rotational
        equations, moving-CG/exhaust terms, force/moment conventions.
    5.  **Initial Conditions and Flight-Phase Definition** — initial state, mission
        sequence, transition conditions.
    6.  **Stage Separation and Child-State Mapping** — configuration change,
        inherited states, conservation-check definitions.
4.  **Environment Models** — atmosphere, gravity, Earth rotation, wind,
    spherical ground, model domains, sources, and assumptions.
5.  **Aerodynamics** — configuration-dependent aerodynamic coefficients, axial
    and crossflow forces, moments, damping, grid-fin aerodynamic response, and
    evidence limits.
6.  **Sensor Models and Measurement Generation** — IMU, GNSS and LiDAR mounting
    locations; ideal measurement equations; noise, bias, scale factors,
    latency, sampling rates, measurement validity. Sensor generation belongs
    here; the algorithm that processes the observations belongs in Chapter 7.
7.  **Navigation Architecture and State Estimation** — estimated states and
    initialization, attitude propagation/interface, filter prediction and
    measurement updates, process and measurement noise, covariance, estimation
    assumptions.
8.  **Mission Execution and Numerical Implementation** — the executed commands
    and numerics:
    1.  **Ascent Commands** — pitch program, attitude commands, engine commands, shutdown envelope, upper-stage continuation.
    2.  **Return Commands** — boostback prediction and steering, entry-burn logic, landing-ignition prediction, fuel budgets, restart conditions.
    3.  **Force and Moment Allocation and Actuator Representation** — engine selection, throttle, requested force direction, virtual/physical control moments, grid-fin commands, limits of the implemented actuator representation. Write this subsection last.
    4.  **Integration, Events, and Sensor Scheduling** — RK4, accepted-state synchronization, scheduled/state-triggered events, contact refinement, guidance clocks, sensor/output timing.
    5.  **Two-Vehicle Execution and Simulation Operation** — independent post-separation propagation, execution options, configuration/seed records, output files, failure handling.
9.  **Verification, Validation, and Simulation Results** — model tests,
    numerical checks, reference-data processing and comparison, flight and
    landing results, navigation errors and consistency, sensitivity analysis.
10. **Engineering Assessment and Implications for Future GNC** — dominant
    uncertainties, the most dangerous model/navigation errors, their
    consequences, demonstrated applicability, and improvements needed for
    future guidance and control.
11. **Conclusions** — achievement of the project objectives, key quantitative
    findings, verified scope, and remaining limitations.

References (unnumbered): combined, deduplicated bibliography of all authors,
verifiable sources only. Appendix
(`content/appendix.tex`, starts with `\appendix`): parameter/assumption/source
register, variable definitions and units, run configurations, long code
excerpts, additional results.

### Placement rules (definitions versus results)

- Chapter 3 defines phase transitions, separation and child-state mapping
  mathematics. Measured event times (including nominal-time phase rows),
  conservation residuals and landing metrics are results and belong in
  Chapter 9; return-planning propellant allocations belong to Chapter 8.
- The aerodynamic response to a given grid-fin deflection stays in Chapter 5;
  the choice of deflection belongs in Chapter 8. Chapter 3 defines thrust and
  flow for a supplied throttle; Chapter 8 explains how throttle is selected.
- Resulting flight profiles (flight-path plot, thrust profile) and the landing
  acceptance table belong in Chapter 9, referenced forward from the defining
  chapters.
- Chapter 10 (Engineering Assessment) is grounded in Chapter 9 numbers:
  dominant uncertainties, the most dangerous model/navigation errors, their
  consequences. Forward-looking improvements are labelled as recommendations
  (class A reasoning), never presented as measured findings; no new ungrounded
  numbers.
- Reference-attitude diagnostic results and gyro-driven navigation results
  must remain distinguishable in Chapter 9 and wherever navigation errors are
  quoted.

### Timing rules

- **Actuator layer** (`docs/flags.md` FL-2): once implemented, re-verify the
  loop chain (Chapter 2), the grid-fin aerodynamic response (Chapter 5) and
  actuation (Chapter 8), re-run the test suites, regenerate the nominal run if
  behavior changed, and write subsection 8.3 last, grounded in the final code.
- **Gyro-based attitude propagation** is implemented: `kalman.py` has an
  optional MEKF attitude/gyro-bias block (active with `q0`), `mission.py` runs
  `GyroAttitude` with the `use_estimated_attitude` setting default true.
  Prose still claiming the filter uses truth attitude
  only is stale: verify against `kalman.py`/`mission.py` and update Chapter 7
  before citing any navigation error figure.

## Conventions

- Headings: `\section` / `\subsection` / `\subsubsection` — no manual numbering.
- Acronyms: `\ac{KEY}` in text; define keys only in `content/acronyms.tex`. The
  file now holds only mission acronyms (the template leftovers AES, CPA, DPA,
  ... were removed); keep it that way and declare a key only when a chapter
  uses it.
- Figures:

  ```latex
  \begin{figure}[htbp]
    \centering
    \includegraphics[width=\linewidth]{images/<name>.pdf}
    \caption{Caption below the figure.}
    \label{fig:<name>}
  \end{figure}
  ```

  Always both `caption` and `label`; reference as `Figure~\ref{fig:<name>}`.
  Every float must be referenced from the text before it appears.
- Tables: `booktabs` (`\toprule`, `\midrule`, `\bottomrule`), see
  `content/sample_content.tex`. Label `tab:<name>`, reference
  `Table~\ref{tab:<name>}`.
- Math: inline `$...$` for short expressions; numbered equations only when
  referenced:

  ```latex
  \begin{equation}
    \label{eq:<name>}
    ...
  \end{equation}
  ```

  Reference with `Equation~\eqref{eq:<name>}`. Define every symbol in the text
  immediately after its first use, or in a notation/summary block.
- Units: SI, in math mode with a thin space, e.g. `$8\,\mathrm{MHz}$`,
  `$9\,\mathrm{km}$`, thousands as `$5\,461$`. Keep units identical across
  chapters ($\mathrm{km}$ vs $\mathrm{m}$: pick once).
- `amsmath`, `graphicx`, `booktabs`, `acro`, and `hyperref` are already
  available. Do not add packages without checking they exist in this TeX
  installation, and only for a documented need.

## Figures

- 8–14 floats across the report; every model chapter needs at least one figure
  or table. Plots of trajectory, loads, navigation errors belong in the results
  chapter.
- Generate plots from real runs, not hand sketches: run the simulator
  (`python main.py --case full --no-show` from the repo root), then plot from
  `results/*/truth.csv`, `measurements.csv`, `navigation.csv`, `events.json`.
  Reopen saved runs with `python main.py --plot-only results/<run>`.
  `docs/image.png` exists as a schematic candidate.
- **Runtime and run outputs.** The full case takes roughly 25–30 min wall
  time — start it in the background and poll its log before relying on it.
  Run it into a fresh output directory: an existing one raises
  `FileExistsError`. A completed run produces `truth.csv` (46 columns incl.
  `time_s`, `phase`, `cg_altitude_m`, `speed_mps`, `mach`,
  `dynamic_pressure_pa`, `angle_of_attack_rad`, `thrust_n`, `mass_kg`,
  `cg_station_m`, `inertia_kgm2_*`, propellant per stage),
  `stage2_truth.csv` + `stage2_events.json` + `event_timeline.csv`
  (separated upper stage), `navigation.csv` (state estimates, errors,
  sigma_0..8), `measurements.csv`, `events.json` (with separation
  conservation audit and landing metrics), `validation.csv` (ascent vs the
  video reference), `summary.json` (status, config, seed, input hashes).
  Sub-sampled CSVs overwrite nothing: truth rows appear every output
  period (0.1 s), not every integration step.
- Export vector PDF where possible (matplotlib: `savefig(..., format="pdf")`).
  PNG only when unavoidable. Store exported figures under `images/` (committed);
  the `results/` folders they were generated from are locally created and
  Git-ignored.
- Captions must be self-contained: a reader who sees only the caption and the
  figure must understand what is shown.

## References

- No bibliography engine (bibtex/biblatex) is configured. Use
  `content/references.tex` with a `thebibliography` environment, registered
  last, and cite with `\cite{key}`.
- References remain an unnumbered section. All authors maintain their own
  sources and merge them into the single deduplicated `thebibliography` for
  submission; each author retains the page, table or figure numbers for the
  material cited and the raw published-input/derived-estimate/assumption
  distinction.
- Cite only verifiable sources: the evidence documents above name the origin of
  each parameter (e.g. Falcon 9 user guide, public telemetry); reuse exactly
  those. Never invent authors, titles, DOIs, or page numbers. Repo-internal
  facts are attributed in text by module path (e.g. `Aerodynamics.py`), not in
  the bibliography.

## Writing style

- Professional, precise, direct. Present tense for how the simulator works;
  the system is the subject, not the authors.
- Be concrete: name the actual component, model, or parameter, give values with
  units, and show the governing equations.
- State assumptions, model limitations, and simplifications explicitly in the
  chapter that describes the affected model — as an explicit assumptions
  paragraph, not a footnote.
- Organize each chapter around one part of the simulator (environment,
  aerodynamics, propulsion, ...); motivational text stays under 3 sentences.
- Support every design decision with an equation, figure, table, or reference.
- Report validation numerically (error percentages, comparisons against
  `data/crs10_ascent_reference.csv` and published Falcon 9 performance), never
  qualitatively ("matches well" → "peak altitude within $1.2\%$ of reference").
- No filler or emphasis words ("very", "robust", "obviously", "seamlessly").
  Short paragraphs, one topic each, topic sentence first.
- After all chapters exist, do one notation review pass: every symbol used in
  more than one chapter must mean the same thing everywhere.

## Build and verification loop

Run from this directory after every chapter and at the end:

```bash
latexmk -pdf main.tex          # build (reruns pdflatex/bibtex as needed)
grep -i "references" main.toc  # -> page listed for References; that page number
                               #    minus 1 is the body page count: must be at
                               #    least 15 (no upper limit)
grep "Output written" main.log # -> total pages incl. references and appendix (info only)
grep -i "undefined" main.log   # -> must be empty (labels, acronyms)
grep -rn -i -e lipsum -e lorem content/ | grep -v sample_content
                              # -> must be empty
texcount -sum -inc content/*.tex
                              # rough word-total gauge
```

If the report quotes test counts, verify them first — from the repo root:

```bash
.venv/bin/python -m unittest discover -s tests -q
.venv/bin/python -m unittest discover -s Environment -p "test*.py" -q
```

`latexmk -c` removes auxiliary files (keeps `main.pdf`). A failed build prints
the error context in the `latexmk` output; fix the `.tex` source, never the
auxiliary files. If `pdfinfo` is used for the page count, run it as
`env -u LD_LIBRARY_PATH pdfinfo main.pdf`; in some launch environments an
injected `LD_LIBRARY_PATH` breaks poppler otherwise.

## Template traps

- The prose docs (`docs/*.md`) lag the code and `data/*.json` in places.
  `../docs/flags.md` is the tracked list of them; read it before citing any
  doc statement, and append flags for discrepancies you discover.
- `main.tex` no longer loads `lipsum` or `svg`; keep it that way. Do not use
  `\lipsum` anywhere in shipped content.
- `content/sample_content.tex` is a syntax reference only. Never register it in
  `main.tex`: its lorem text and crypto figures would end up in the PDF.
- `Actuators.py` is an empty placeholder **at the time of this writing, but
  the actuator layer is planned and will be implemented before the report is
  submitted**. Treat it as a first-class part of the system:
  - The report structure must cover it (chapters 2, 5 and 8 in the binding
    outline; force and moment allocation is subsection 8.3) —
    do not leave it out or mention it only in "future work".
  - Verify the actual implementation state of `Actuators.py` (and its tests)
    before writing the actuator subsections, and write those subsections
    last, when the final code exists.
  - Ground every actuator claim in the implemented code and tests, as for
    every other module. Never describe behavior that does not exist yet;
    scope to what the implementation actually covers (e.g. TVC gimbal
    dynamics, grid-fin rate/travel limits).
  - `docs/aero_model.md` and `docs/core_api.md` list the actuator layer
    under "Open items" / "not implemented" — reconcile those statements
    with the final code before citing them.
- Do not manually number sections, figures, tables, or equations.
- Do not edit `tum/`, `tumuser.sty`, or `latexmkrc`.

## Definition of done

All of the following must hold before the report is considered finished:

- [ ] `main.pdf` compiles with `latexmk -pdf main.tex`, zero errors.
- [ ] Body page count (title/ToC through Conclusions) is at least 15; there
      is no upper limit (gate above). References and appendix sit outside
      the gate.
- [ ] Chapter structure matches the binding outline: one file per chapter,
      `main.tex` input order matches the outline, no interim files are
      registered, every chapter is registered with `\include` (which starts
      it on a new page), and
      result material (measured times, residuals, flight profiles, landing
      metrics) lives in Chapter 9.
- [ ] Results chapter distinguishes reference-attitude diagnostics from
      gyro-driven navigation results, and all Chapter 10 claims trace to
      Chapter 9 numbers or are labelled recommendations.
- [ ] No undefined references, no lipsum/lorem remnants, no orphaned floats,
      no manual numbering.
- [ ] Every numeric claim is traceable to a source in the source map; every
      citation points to a real, verifiable source.
- [ ] Each model chapter states its assumptions and limitations.
- [ ] Results chapter contains quantitative validation against the CRS-10
      reference data.
- [ ] Notation is consistent across all chapters (final review pass done).
- [ ] Placeholder packages removed, acronym file holds only used acronyms,
      leftover build artifacts cleaned (`latexmk -c`).
