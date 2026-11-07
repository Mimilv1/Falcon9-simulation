# Falcon9-simulation

Project 2 CRS-10 vehicle model, A dynamics core and temporary B load models.
Team mission conditions and the integrated flight runner are pending.

## Current team convention

All shared position, velocity, force, moment, gravity and angular-rate vectors
use **current Body components**. **Body +x points toward the nose**. The frame is
right handed; engine-ring angle zero is +y and +90° is +z. Geographic heading
remains unspecified.

The fixed reference frame is **ECEF**: the Earth centre is the origin (0,0,0) and
the vehicle launches from **(0,0,earth_radius)** on the +z axis. Earth rotation
is not modeled in v0, so the ECEF axes are treated as non-rotating (inertial).
The quaternion is Hamilton [w,x,y,z], mapping current Body to ECEF. At the pad
the vehicle is vertical, Body +x along ECEF +z: initial
q=[sqrt(0.5), 0, -sqrt(0.5), 0] (`frames.LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ`).

`State.position_body_m` means the vector **from the Earth centre to CG,
expressed in current Body axes**. It is not CG position relative to itself.
Position and velocity derivatives include rotating-basis transport terms.

## Files

| Files | Team use |
| --- | --- |
| `VehicleData.py` | Engine, Stage and Vehicle objects; base/CG mapping |
| `vehicle_parameters.py`, `data/vehicle_parameters.json` | Adopted values and K/E/A source records |
| `vehicle_mass_model.py`, `mass_properties.py` | Fuel inventory, mass, CG and inertia |
| `load_interface.py`, `mass_flow.py` | Force/moment assembly and exhaust data contract |
| `dynamics.py`, `propagation.py`, `frames.py` | Motion equations, integration and Body attitude |
| `propulsion.py`, `Aerodynamics.py`, `Environement.py` | Temporary thrust, multi-configuration aerodynamics (axial drag, normal force, station moments, rotary damping) and environment models for B to refine |
| `ascent.py` | Reusable temporary load provider connecting B models to A; no flight scenario |
| `engine_schedule.py`, `staging.py` | Parameterized shutdown-envelope utility and no-impulse separation map; no mission timing |
| `main.py` | Integration entry point; currently loads the vehicle without running a mission |
| `tests/` | Component and mathematical regression checks |

## Use

Python 3.11+ and NumPy are required:

```text
python -m pip install -r requirements.txt
python main.py
python -m unittest discover -s tests -v
```

`main.py` only loads the vehicle and reports integration status. There is no
default trajectory, injected side force or hardcoded flight-event timeline.
Private flight checks, configurations and generated results are not included.

## Integration still required

- A: connect the team's load provider and event conditions to `propagate`.
- B: refine thrust, aerodynamic and environment models; provide consistent
  Body forces, CG moments, gravity and exhaust/consumption data.
- C: use fixed-sensor truth to generate measurements and estimate states.
- Team: select representative phases and one reference flight/data set; define
  cutoff, separation and restart conditions from evidence or explicit assumptions.

CRS-10 is the default configuration; Dragon mass is estimated and CG/inertia are
assumed from an equivalent cylinder. Stage/engine values remain mixed-source. Remaining fuel at the numerical tank limit is not a known
landing reserve. Mission timing, separate-stage loads and free-flight liquid
behavior have not been agreed or validated.

## References

- [A/B/C interface](docs/core_api.md)
- [Mass and engine parameter evidence](docs/parameter_evidence.md)
- [Powered dynamics equations and limitations](docs/powered_dynamics.md)
- `data/drag_curve_pezent.json`: temporary drag-table evidence and readout method
- `data/drag_curve_pezent_extended.json`: the same curve with appended Mach 8/10
  declared-estimate knots so full-mission ascents pass the Mach 7 evidence ceiling
- `data/aero_configurations.json`: declared per-configuration aerodynamic records
  (full stack, booster alone, booster with grid fins and legs, second stack);
  only the full-stack axial curve is sourced, all other coefficients are documented estimates
