# Falcon9-simulation

CRS-10 vehicle simulator with IMU, GNSS, LiDAR and a navigation filter.

## Run

Python 3.11+. In the repository folder:

```text
git switch main
git pull --ff-only origin main
python -m pip install -r requirements.txt
python main.py --case full --no-show
```

The full case runs from liftoff to first-stage landing. Shorter cases:

```text
python main.py --no-show
python main.py --case terminal --no-show
```

The default case ends at 160 s. The terminal case is an independent 10 s
powered-descent test, not a continuation of ascent.

## Results

The terminal prints the output directory under `results/`.

| File | Contents |
| --- | --- |
| `overview.png` | Flight and navigation-error plots |
| `truth.csv` | Simulated vehicle state and loads |
| `measurements.csv` | Sensor samples |
| `navigation.csv` | Estimated states and errors at the IMU station |
| `events.json` | Phase transitions, separation audit and landing metrics |
| `summary.json` | Outcome, configuration, sources and input hashes |
| `validation.csv` | Comparison with the ascent reference |

Successful full runs print `landing_success` and
`Full return flight completed: True`. Results are generated locally and ignored
by Git. Reopen a saved run with:

```text
python main.py --plot-only results/YOUR_RUN_FOLDER
```

## Coordinates

Body +x points toward the nose. Shared forces, moments, angular rates and
position/velocity use current Body components. Position is Earth centre to CG;
velocity is relative to the rotating Earth. Component stations are measured
from the structural base. Inertia and moments are about the current CG.

The external frame is Earth-centred and launch-aligned, called ECEF in the
code. Initial base: [0,0,6371000] m. Local north is +x, east −y, up +z.
It is not conventional Greenwich-axis ECEF. The initial Hamilton Body-to-ECEF
quaternion is [√0.5,0,−√0.5,0]. Gravity and rotating-frame accelerations are
supplied separately and added once each.

## Main files

| File | Role |
| --- | --- |
| `main.py` | Run a case or reopen results |
| `mission.py` | Connect models, events, sensors, navigation and output |
| `VehicleData.py` | Vehicle, stage and engine objects |
| `vehicle_parameters.py`, `data/vehicle_parameters.json` | CRS-10 parameter set |
| `vehicle_mass_model.py`, `mass_properties.py`, `mass_flow.py` | Mass, CG, inertia and propellant flow |
| `dynamics.py`, `frames.py`, `propagation.py` | Motion equations, transformations and integration |
| `load_interface.py` | Forces and moments about CG |
| `ascent.py`, `engine_schedule.py`, `recovery.py` | Flight commands and load allocation |
| `propulsion.py`, `Aerodynamics.py`, `Environement.py` | Thrust, aerodynamics, atmosphere, gravity and Earth rotation |
| `Events.py`, `staging.py`, `data/crs10_mission.json` | Phase conditions and separation |
| `sensors.py`, `kalman.py` | Sensor models and navigation filter |

## Model details and sources

- [Mission phases, settings and verified result](docs/crs10_mission.md)
- [A/B/C interface](docs/core_api.md)
- [Vehicle parameters and evidence](docs/parameter_evidence.md)
- [Aerodynamics and evidence](docs/aero_model.md)
- [Moving-mass dynamics](docs/powered_dynamics.md)

The navigation filter receives truth attitude. The return controller uses
truth-state feedback and ideal attitude moments. Sensor specifications,
landing geometry and several vehicle/aerodynamic values are research-model
assumptions, as identified in the linked documents.

## Tests

```text
python -m unittest discover -s tests -q
python -m unittest discover -s Environment -p "test*.py" -q
```
