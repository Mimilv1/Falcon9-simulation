# Falcon9-simulation

CRS-10 research simulator. The default run connects the vehicle, B environment
and multi-configuration aerodynamics, A moving-mass six-DOF dynamics and events,
and the existing C IMU/GNSS/LiDAR and navigation filter.

## Run

Python 3.11+:

```text
python -m pip install -r requirements.txt
python main.py --no-show
python main.py --case terminal --no-show
python main.py --case full --no-show
```

The default ascent runs to 160 s: powered ascent, smooth shutdown (140–141 s),
coast, stage separation (144 s), then booster-only coast. It stops before the
161 s boostback reference. The `full` case connects boostback, entry, descent
and leg contact. The revised nominal run with the default 100 Hz IMU passes
all five landing criteria: vertical speed 0.497 m/s, horizontal 0.902 m/s,
target error 3.306 m, tilt 0.535° and angular rate 0.013°/s. A 20 Hz
engineering run and a 0.025 s dynamics check also pass.
These are simulated acceptance results under the documented assumptions,
not validation of the actual CRS-10 landing.

Fast full integration check using the existing C models:

```text
python main.py --case full --imu-rate 20 --dt 0.05 --no-show
```

The terminal case is a separate 10 s assumed powered-descent fixture, not a
continuation of the ascent trajectory and not a landing success demonstration.

Outputs go to a new timestamped directory under `results/` (ignored by Git).
`overview.png` shows the result. To choose a location or reopen a saved run:

```text
python main.py --output results/team-ascent --no-show
python main.py --plot-only results/team-ascent
```

Use `--until 20` for a short check. `--dt` must be at most 0.01 s with the
configured 100 Hz IMU; dynamics-only checks may use `--no-navigation`.
An incomplete run saves its accepted results and returns exit code 2.
Existing run files are not overwritten.

## Team convention

Body +x points toward the nose. Shared forces, moments, angular rates and
state position/velocity use current Body components. Structural stations are
measured from the vehicle base; CG moments and inertia are about current CG.

The external frame is Earth-centred and launch-aligned, called ECEF in the
code: initial base position [0,0,6371000] m; local north +x, east −y, up +z.
It is not conventional Greenwich-axis ECEF. Initial Body +x maps to external
+z with Hamilton Body-to-ECEF quaternion [√0.5,0,−√0.5,0]. State position is
Earth-centre → CG in current Body axes, not position relative to CG itself.
Velocity is relative to the rotating Earth, expressed in Body axes.

B supplies inverse-square gravity, Coriolis and centrifugal accelerations
separately; A adds each once. The team's present rotational dynamics/gyro
model does not include a complete Earth-rate correction; see the model limits.

## Files

| File | Purpose |
| --- | --- |
| `main.py` | Single command-line entry point |
| `mission.py` | Phase runner, existing C-model adapter, CSV/JSON export, validation comparison and plots |
| `Events.py`, `data/crs10_mission.json` | Ordered events, guards, source-timed commands and explicit assumptions |
| `VehicleData.py`, `vehicle_parameters.py`, `data/vehicle_parameters.json` | Vehicle objects and the single CRS-10 parameter set |
| `vehicle_mass_model.py`, `mass_properties.py`, `mass_flow.py` | Fuel, moving CG, inertia and exhaust contract |
| `dynamics.py`, `frames.py`, `propagation.py` | Equations, quaternion transforms and integration utilities |
| `recovery.py` | Project-designed return guidance, force/moment allocation, fuel budgets and geometric contact |
| `ascent.py`, `engine_schedule.py` | B load adapter, retained ideal pitch command and smooth shutdown |
| `propulsion.py`, `Aerodynamics.py`, `Environement.py` | B thrust, aerodynamics, gravity, atmosphere and rotation effects |
| `staging.py` | Zero-impulse separation and conservation audit |
| `sensors.py`, `kalman.py` | Existing C sensor and navigation algorithms |

## Limits and next work

- The published approximate preflight timeline is adopted as simulation commands;
  these are not measured actual event times or flight-computer thresholds.
- Pitch tracking uses an unlimited pure couple, not a physical TVC actuator.
- Full-stack axial drag was obtained for a generic fairing geometry; Dragon
  applicability and booster/high-Mach coefficients remain uncertain.
- Navigation receives truth attitude; sensor settings remain provisional.
- Return guidance and contact detection are implemented as explicit A-class
  project designs. Landing performance, physical actuator realism and robustness
  still require validation; the controller uses truth-state feedback.
- Stage-1 nominal operating level is based on the 2015 guide (756 kN SL per
  engine). Applying it to CRS-10 is E/A; it does not identify actual rating/throttle.
- Pitch tau=100 s was calibrated against the same ascent data; improved agreement
  is calibration, not independent validation.

## Checks and evidence

```text
python -m unittest discover -s tests -v
python -m unittest discover -s Environment -p "test*.py" -v
```

- [Mission execution, assumptions and handoff](docs/crs10_mission.md)
- [Shared A/B/C interface](docs/core_api.md)
- [Aerodynamic model and evidence](docs/aero_model.md)
- [Vehicle parameter evidence](docs/parameter_evidence.md)
- [Moving-mass equations](docs/powered_dynamics.md)

Reference comparison uses `data/crs10_ascent_reference.csv` and its source
metadata JSON. It is published webcast digitization (E), not official onboard
telemetry. Simulation navigation error is an internal estimation metric,
separate from disagreement with the flight reference.
