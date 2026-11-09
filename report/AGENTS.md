# Report Guide (for AI agents)

Operating manual for writing and editing the technical report in this directory.
Read this file fully before touching any `.tex` file. The workflow and conventions
here are binding; skip them and the result will be rejected at the final check.

## Deliverable and hard requirements (read first)

1. **Scope**: technical report documenting the Falcon 9 / CRS-10 simulator in the
   repository root. It documents the system — models, implementation, results —
   not the process of writing the report.
2. **Length: 14–20 pages in the compiled `main.pdf`, figures included.** The title
   page counts; so does a table of contents if one is added. Aim for 16–18 to leave
   margin on both sides. Re-check the page gate (below) after every structural
   change; never discover the page count at the end.
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
   Never start a chapter from memory or plausibility alone.
2. **Outline first.** Fix the chapter list and a per-chapter page budget
   (suggested outline below), sum the budget, and confirm the total lands in
   14–20 pages. Register all `\input` lines before writing content.
3. **Write one chapter at a time.** Finish a chapter including its equations,
   floats, and citations, then run the build and fix warnings before moving on.
4. **Generate figures from real data**, not sketches (figure rules below).
5. **Iterate with the verification loop** (below) after every chapter: build,
   page count, warning check, placeholder check.
6. **Final pass**: notation review across all chapters (same symbol = same
   meaning everywhere), then run the full definition-of-done checklist.

## Layout

- `main.tex` — master document: metadata, package loading, ordered `\input` list.
- `content/` — one file per chapter, starting with `\section{...}`.
- `images/` — all figures (PDF preferred). Create it; it does not exist yet.
- `content/acronyms.tex` — all acronyms via `\DeclareAcronym{KEY}{...}`.
- `content/sample_content.tex` — syntax reference for tables, figures, headings.
- `tum/`, `tumuser.sty`, `latexmkrc` — class and build support. Never edit these.

## How to add a chapter

1. Create `content/<name>.tex`, starting with `\section{Chapter Title}`.
2. Register it in `main.tex` below the `% +++++ Main part` marker, in PDF order:

   ```latex
   \input{content/introduction}
   %\input{content/next_section}
   ```

3. Build with `latexmk -pdf main.tex` from this directory.

## Source map: where the truth lives

Read these before writing the matching chapter. Numbers in the report must be
traceable to one of these (or to the code they describe).

| Report chapter | Sources |
| --- | --- |
| Mission, phases, settings, result | `docs/crs10_mission.md`, `data/crs10_mission.json`, `mission.py`, `Events.py` |
| Vehicle, mass, parameters | `README.md`, `docs/parameter_evidence.md`, `data/vehicle_parameters.json`, `VehicleData.py`, `mass_properties.py`, `mass_flow.py` |
| Environment (atmosphere, gravity, wind) | `Environment/environment.md`, `Environement.py`, `Environment/test_environment.py` |
| Aerodynamics | `docs/aero_model.md`, `data/aero_configurations.json`, `data/drag_curve_pezent*.json`, `Aerodynamics.py`, `tests/test_aerodynamics.py` |
| Propulsion, engine schedule | `docs/parameter_evidence.md`, `propulsion.py`, `engine_schedule.py`, `mass_flow.py` |
| Dynamics, frames, propagation, loads | `docs/powered_dynamics.md`, `docs/core_api.md`, `dynamics.py`, `frames.py`, `propagation.py`, `load_interface.py` |
| Staging, ascent control, recovery | `staging.py`, `ascent.py`, `recovery.py`, `tests/test_staging.py`, `tests/test_recovery.py` |
| Sensors and navigation filter | `sensors.py`, `kalman.py`, `README.md` ("Coordinates" section) |
| Validation and results | `data/crs10_ascent_reference.*`, `data/ascent_model_audit.json`, `docs/crs10_mission.md`, `tests/` |

The coordinate conventions from the root `README.md` (Body +x to the nose;
ECEF is launch-aligned, not Greenwich) are binding: do not silently rename them.

## Suggested outline and page budget

The report tells three stories in order: how the simulator is built (loop,
frames, models), what it simulates (vehicle, stages, flight paths), and how
well it does it (results and validation). Adjust depth to actual content, but
keep the total in 14–20 pages:

| # | Chapter | Pages | Core content |
| --- | --- | --- | --- |
| 1 | Introduction & scope | 0.75–1 | CRS-10 context, simulator purpose, the K/E/A evidence-class convention used throughout |
| 2 | Simulator architecture & loop | 1.5 | A/B/C contract table; the mission loop: event machine → step-size selection → RK4 → state sync; the control chain: guidance → actuator layer → applied forces/moments; step-cap priority list (dt 0.01 s, event boundaries, controller knots, output period, sensor sample times, landing clearance); guards (propellant limit, domain rejection, touchdown bisection); reproducibility (seed, input SHA-256 hashes, run outputs). Loop flow diagram |
| 3 | Vehicle, stages & mission | 1.5–2 | Stage/geometry/dry-mass table with evidence classes; propulsion inventory; propellant budgets; mission phase timeline table; the three flight paths: ascent arc, separated upper stage, booster return arc (trajectory sketch) |
| 4 | Reference frames & dynamics | 1.5 | Body +x, launch-aligned ECEF, quaternion transport terms; moving-CG linear/angular balances; sensor-station kinematics; RK4 Earth-centred safeguard |
| 5 | Environment models | 1 | USSA layers + tail, point-mass gravity + gradient, Coriolis/centrifugal frame acceleration, wind profile |
| 6 | Aerodynamics | 1.5 | Four configurations, axial CA0(M) curve figure, normal/crossflow/damping equations, grid fins, evidence classes |
| 7 | Propulsion & mass model | 1.5 | Thrust/Isp relations, sea-level conversion, mass flow, O/F, tank depletion, CG/inertia migration figure |
| 8 | Ascent, staging & upper stage | 1.25 | Attitude/pitch program, event conditions, smooth shutdown, separation conservation audit, stage-2 continuation |
| 9 | Booster return, control & actuation | 1.5–2 | Boostback predictor, entry/landing logic, PD gains table, moment limits, grid-fin interplay, landing criteria; **actuator layer** (command path from the return controller through actuator dynamics to TVC/grid-fin forces and moments — scope to the final implementation). Write the actuator subsection last |
| 10 | Sensors & navigation | 1.25–1.5 | IMU/GNSS/LiDAR spec table, 9-state Kalman filter, truth-attitude limitation |
| 11 | Results & validation | 2–2.5 | Ascent vs video reference, flight profiles (q̄, Mach, AoA), separated trajectories figure, landing metrics table, navigation errors |
| 12 | Conclusion & limitations | 0.5–0.75 | Evidence-class summary, declared estimates and open items |

Chapters 2 and 3 deliver the "how the simulation works" and "what the
stages and flight paths are" basics the report must teach before any model
detail; chapters 4–9 then go deeper per subsystem, and 10–11 show the
result. Sum of midpoints ≈ 16 pages — over budget chapters get trimmed, not
padded. One timing rule: the actuator layer is being implemented while this
report is written. Once it lands, re-verify the affected facts (loop chain
in chapter 2, grid-fin command path in chapter 6, actuation in chapter 9),
re-run the test suites, regenerate a nominal full run if the simulation
behavior changed, and write those subsections last.

Rule of thumb: roughly 450 words of body text per page in this layout, i.e.
5,000–7,000 words overall; floats and equations take the rest. Treat `texcount`
(below) as the gauge and the compiled page count as the criterion.

## Conventions

- Headings: `\section` / `\subsection` / `\subsubsection` — no manual numbering.
- Acronyms: `\ac{KEY}` in text; define keys only in `content/acronyms.tex`. The
  current acronym file still holds template leftovers (AES, CPA, DPA, ...) —
  replace them with mission acronyms; do not carry unused crypto acronyms into
  the final report.
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
grep "Output written" main.log # -> "(N pages, ...)": N must be 14-20
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

- `main.tex` loads `lipsum` (placeholder filler) and `svg` (needs Inkscape).
  Do not use `\lipsum` anywhere in shipped content; remove the unused
  `\usepackage{lipsum}` and `\usepackage{svg}` from `main.tex` before the final
  build once figures are in place.
- `content/sample_content.tex` is a syntax reference only. Never register it in
  `main.tex`: its lorem text and crypto figures would end up in the PDF.
- `content/introduction.tex` currently contains lorem ipsum — replace it fully.
- `Actuators.py` is an empty placeholder **at the time of this writing, but
  the actuator layer is planned and will be implemented before the report is
  submitted**. Treat it as a first-class part of the system:
  - The report structure must cover it (chapters 2, 6, 9 in the outline) —
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
- [ ] Page count is 14–20 (grep gate above).
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
