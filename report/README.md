# CRS-10 report

`main.pdf` is the compiled submission. Chapters 1–11 are registered with `\include` in `main.tex`; references and Appendix A follow the conclusions. `content/references.tex` contains the bibliography directly; no BibTeX step is required.

## Reproduce the submitted figures and checks

Run from the repository root, with Python 3.12 and the pinned dependencies:

```sh
python -m pip install -r requirements.txt
python report/reproduce.py --unpack --validate
python -m unittest discover -s tests -q
python -m unittest discover -s Environment -p "test*.py" -q
```

`reproduction_data.zip` contains the final nominal run, the paired earlier navigation baseline, and the coarse step/IMU-rate variant. The script verifies archive/member SHA-256 hashes, extracts to `results/submission-reproduction`, and regenerates flight/navigation figures, the architecture, mass sweep, mass-code excerpt, and numerical macros. It refuses to overwrite differing extracted files. Published image crops and unchanged environment/aerodynamic illustrations remain supplied assets with attribution.

`figure_sources.json` records runtime input fingerprints, result metrics, archive members and generated figure fingerprints. `navigation_validation.json` records independent accelerometer calibration, nonzero gyro calibration, 20-trial position ANEES, and isolated RK4 convergence checks. The full independent checks can also be run alone:

Original runtime byte hashes remain in the saved summary. Source verification also accepts the separately recorded LF-normalized text hash, so Git's LF/CRLF checkout conversion is distinguished from an actual content change. Archive result files retain exact byte-hash checks.

```sh
python report/validate_navigation.py
```

The final reviewed full run uses 30 error-state coordinates, seed 1, a 0.01 s maximum physical step and a 100 Hz IMU. Generate a new integration with:

```sh
python main.py --case full --no-show --seed 1 --dt 0.01 --imu-rate 100
```

New integrations use fresh result directories. The ascent-only case ends at 160 s. To inspect the archived full case without integration:

```sh
python main.py --plot-only results/submission-reproduction/nominal --no-show
```

## Build the PDF

From `report/`, with a XeLaTeX-capable TeX installation:

```sh
latexmk -xelatex main.tex
```

The supplied `latexmkrc` locates the TUM support directory; `-xelatex` selects the engine used for this report. An alternative is Tectonic:

```sh
tectonic -Z search-path=tum main.tex
```

Keep `tum/`, `tumuser.sty` and `latexmkrc` unchanged. Each chapter starts on a new page; internal figures use the available page space. The monospaced font disables smart-quote mapping so copied executable excerpts retain ASCII quotes.

The report distinguishes ascent calibration from independent flight validation, nominal flight navigation from controlled fixtures, and ideal truth-feedback command results from navigation-fed or hardware-actuator performance.
