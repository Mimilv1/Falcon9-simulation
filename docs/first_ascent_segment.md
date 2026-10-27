# First connected A+B segment: early vertical powered ascent

The first project run covers **t = 0–30 s of vertical powered ascent**, starting
with nine first-stage engines already at steady full throttle. This is a short
part of the ascent phase, not the complete ascent or a reconstructed mission.
The full stack remains attached. The engine command is fixed throughout; t = 30 s
is the selected simulation endpoint, not engine cutoff.

The implementation is deliberately replaceable: B already supplied useful engine
and vehicle data. The new environment, propulsion and aerodynamic functions turn
those inputs into loads that A evaluates at every integration trial. A updates
fuel, mass, CG, inertia, position, velocity and attitude together. C receives ideal
sensor signals from the same accepted states.

## Run and inspect

From the repository folder, with Python 3.11+ and the requirements installed:

```text
python main.py
python main.py --step-s 0.01 --output-dir results/ascent_finer
python -m unittest discover -s tests -v
```

The default inputs are `config/powered_ascent.json`. Output is reproducible from
that file, the vehicle profile and drag table. `results/` is ignored by git.

| Output | Purpose |
| --- | --- |
| truth_for_c.csv | Accepted-state truth and ideal sensor signals, nominally every 0.02 s (50 Hz) |
| summary.json | Initial/final values, maxima, stop reason and assumptions |
| run_config.json | Exact effective phase inputs |
| vehicle_sources.json | Vehicle K/E/A evidence and source links |
| drag_curve_sources.json | Numerical Mach table, source PDF hash and digitization method |
| manifest.json | Git revision, dirty-worktree flag, model file hashes and frame/sensor convention |

The integration stops at the requested endpoint or the limiting propellant
inventory. An unsupported atmosphere height or Mach number raises an error;
the model does not silently extrapolate beyond its stated domain.

## Models and evidence

K = authoritative source value; E = literature estimate; A = explicit modeling
assumption. E/A indicates a literature-based model with stated approximations.

| Quantity / model | Adopted value or equation | Class and basis |
| --- | --- | --- |
| Segment and release state | 30 s; stationary structural base at ground; initial Body +x vertical; initial quaternion [1,0,0,0] | A: a selected integration exercise; nominal literature fuel loading, no ignition or hold-down transient |
| Engines | 9 first-stage engines, throttle = 1, zero gimbal; stage 2 off | A: steady command; engine parameters E from Pezent Table 2.16 |
| Per-engine vacuum data | Tᵥ = 914,000 N; Ispᵥ = 311 s; Aₑ = 0.68 m² | E: [Pezent, printed p.77 / PDF p.94, Table 2.16](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) |
| Propellant flow | q = τTᵥ / (Ispᵥg₀); constant q in this segment | E/A: fixed vacuum Isp approximation; τ is throttle and g₀ = 9.80665 m/s² |
| Thrust | T = τTᵥ − Aₑp, with local nozzle ambient pressure p | E/A: fixed exit area and vacuum performance; pressure effect changes as the rocket climbs |
| Atmosphere | T₀ = 288.15 K; p₀ = 101,325 Pa; lapse = 0.0065 K/m; R = 287.05287 J/(kg·K); γ = 1.4 | E/A: [US Standard Atmosphere 1976](https://ntrs.nasa.gov/citations/19770009539), dry troposphere over 0–11 km; geometric height approximates geopotential height |
| Atmosphere equations | T = T₀ − 0.0065h; p = p₀(T/T₀)^(g₀/(R·0.0065)); ρ = p/(RT); a = √(γRT) | Standard hydrostatic/ideal-gas model, not measured weather; checked against [published reference values](https://atmos.nmsu.edu/planetary_datasets/earth_temppres.html) |
| Drag coefficient | Piecewise-linear Cᴅ(M), Mach 0–7; approximately 0.1124 at M = 0.4 | E/A: raster readout of [Pezent Fig.2.13(b), printed p.76 / PDF p.93](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content); RasAeroII at zero AoA, 5.2 m fairing. Original numerical data unavailable. Pixel scales ≈ 0.009 Mach and 0.005 Cᴅ; transonic readout more uncertain. These are readout scales, not statistical error bounds |
| Drag area | S = 21.08 m², retained with the source Cᴅ convention | E: Pezent Table 2.18, printed p.77; do not change area without reconciling the coefficient |
| Aerodynamic load | D = ½ρV²SCᴅ, opposite air-relative velocity; applied through current CG | E/A: [drag equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/); applying this zero-AoA source to the mixed-source vehicle is A. No lift/normal force or aerodynamic moment |
| Gravity / wind | Constant magnitude 9.80665 m/s² / zero wind | A: gravity direction is rotated into current Body each trial; no Earth-rotation model |
| Exhaust boundary | Fixed coplanar exits at axial x = 0; uniform exit disks | A: engine-ring layout from the existing profile. Exit station is distinct from engine force pivot |
| Virtual IMU | Fixed structural station [0,0,0] m | A: base-mounted virtual sensor, not a surveyed Falcon 9 installation |
| Truth interval | 0.02 s, nominally 50 Hz | A: selected simulation interval; C chooses its actual sensor sampling rates |

The liquid, mass and inertia inputs retain their separate evidence in
[the parameter table](parameter_evidence.md). This combined model mixes public
research estimates and explicit assumptions; agreement with internal equations
does not establish agreement with Falcon 9 flight telemetry.

## C's signals

All vector columns use **current Body components, +x toward the nose**. Position
is the vector from the fixed initial origin to instantaneous CG. Velocity is CG
velocity expressed in the rotating Body basis. The quaternion maps current Body
to the fixed initial Body basis. No ENU/ECEF output is introduced. Scalar altitude
uses the scenario's initial vertical direction; attitude/gravity bookkeeping still
requires the internal initial basis.

| CSV column family | Meaning / units |
| --- | --- |
| position_body_m_* / velocity_body_mps_* | CG truth in m / m/s, not a direct GPS measurement |
| quaternion_body_to_initial_* | Hamilton [w,x,y,z] attitude truth relative to initial Body |
| gyro_ideal_body_radps_* | Ideal gyro angular rate in rad/s; Earth rotation omitted |
| accelerometer_ideal_body_mps2_* | Ideal specific force at the fixed IMU station, in m/s²; gravity excluded, CG migration and rigid-body lever-arm effects included |
| acceleration_body_mps2_* | Geometric CG inertial acceleration in Body components; different from accelerometer output |
| pressure_imu_pa | Ideal ambient pressure at the virtual sensor station, Pa; C may build a barometer model |
| base_altitude_m / cg_altitude_m | Scalar heights of two different points above the adopted zero datum |
| stage1_propellant_kg / rp1_kg / lox_kg | Remaining model propellant, kg |
| cg_station_x_m / inertia_body_kgm2_* | CG station from structural base, m / diagonal CG inertia, kg·m² |

C can add bias, noise, quantization and sampling to the ideal sensor columns now.
The navigation algorithm consumes the resulting measurements; continuous position,
velocity and attitude truth should be used only for evaluation or initialization
explicitly allowed by the project. This perfectly vertical, symmetric case has
almost zero angular rate, so it is useful for sensor plumbing and basic checks.
It is not sufficient to validate attitude estimation or full navigation.

## Next A/B milestones

1. Agree on the fixed sensor station and measurement-column contract with C.
2. Add a controlled attitude-changing case and B's normal-force/moment/gimbal
   models, then verify signs, moments and sensor lever arms.
3. Add a phase scheduler, justified ignition/cutoff flow histories and staging
   state maps; preserve the existing fuel-boundary checks.
4. Extend atmosphere/aerodynamic domains before extending the trajectory, then
   compare against suitable independent references and mission data.

The first connected segment removes the need to wait for every B model or all
flight phases before C can begin. Each refinement can replace a provider while
retaining A's integrator and C's interface.
