# Proposed CRS-10 Report Structure

This structure combines Philipp's existing draft with Kwanwoo's vehicle and dynamics material, while separating sensor measurement generation from navigation. The report follows the simulator's data flow: define the vehicle and shared conventions, explain the physical and measurement models, describe mission execution, then evaluate the results and their engineering implications.

## Chapter outline and proposed responsibilities

| Chapter | Title | Main content | Proposed lead |
| --- | --- | --- | --- |
| 1 | Introduction and Project Scope | Project 2 objectives, why CRS-10 was selected, modeled flight scope, and validation objectives. | Kwanwoo |
| 2 | Simulation Architecture and Shared Conventions | Overall data flow, truth versus measurements versus estimates, coordinate frames, quaternion convention, reference points, units, and module interfaces. | Kwanwoo |
| 3 | Vehicle Configuration, Dynamics, and Mission Definition | Vehicle configuration and K/E/A evidence; propulsion and propellant depletion; mass, CG and inertia; equations of motion; initial conditions and phase definitions; separation and child-state mapping. | Kwanwoo, combining both drafts |
| 4 | Environment Models | Atmosphere, gravity, Earth rotation, wind, spherical ground, model domains, sources, and assumptions. | Philipp |
| 5 | Aerodynamics | Configuration-dependent aerodynamic coefficients, axial and crossflow forces, moments, damping, grid-fin aerodynamic response, and evidence limits. | Philipp |
| 6 | Sensor Models and Measurement Generation | IMU, GNSS and LiDAR mounting locations; ideal measurement equations; noise, bias, scale factors, latency, sampling rates, and measurement validity. | Emilien |
| 7 | Navigation Architecture and State Estimation | Estimated states and initialization, attitude propagation/interface, filter prediction and measurement updates, process and measurement noise, covariance, and estimation assumptions. | Emilien |
| 8 | Mission Execution and Numerical Implementation | Ascent and return commands, force/moment allocation and actuator representation, integration, event handling, sensor scheduling, two-vehicle execution, and run/output procedures. | Kwanwoo, using Philipp's existing flight/control material |
| 9 | Verification, Validation, and Simulation Results | Model tests, numerical checks, reference-data processing and comparison, flight and landing results, navigation errors and consistency, and sensitivity analysis. | Joint, assembled by Kwanwoo |
| 10 | Engineering Assessment and Implications for Future GNC | Dominant uncertainties, the most dangerous model/navigation errors, their consequences, demonstrated applicability, and improvements needed for future guidance and control. | Joint |
| 11 | Conclusions | Achievement of the project objectives, key quantitative findings, verified scope, and remaining limitations. | Joint |
| — | References | Course material, official mission documents, papers, public flight data, and applicable sensor references. | Everyone maintains their own sources; combined bibliography |
| Appendix | Reproducibility and Supporting Data | Complete parameter/assumption/source register, variable definitions and units, run configurations, detailed code excerpts, and additional results. | Shared |

References remain unnumbered. Sensor generation belongs in Chapter 6; the algorithm that processes those observations belongs in Chapter 7.

## Chapter 3 integration

Combine the overlapping vehicle and dynamics material into one chapter with a single parameter table and one explanation of each model.

1. **Vehicle Configuration and Evidence Classification** — CRS-10 stage/Dragon configuration, geometry and mass budget, K/E/A definitions, and parameter provenance.
2. **Propulsion and Propellant Depletion** — thrust, ambient-pressure correction, specific impulse, mass flow, loading versus consumption O/F, and usable inventory.
3. **Mass, Center of Gravity, and Inertia** — dry-component and liquid models, component summation, and property derivatives during fuel consumption.
4. **Equations of Motion and Load Interfaces** — translational and rotational equations, moving-CG/exhaust terms, force/moment conventions, and structural-base/CG relationships.
5. **Initial Conditions and Flight-Phase Definition** — initial vehicle state, overall mission sequence, and transition conditions.
6. **Stage Separation and Child-State Mapping** — configuration change, inherited first-/second-stage states, and conservation-check definitions.

The phase definitions and separation mathematics are explained here. Their execution in the mission loop belongs in Chapter 8, and the measured simulation event times and conservation residuals belong in Chapter 9.

## Chapter 8 and the return-control material

Philipp's booster return and control material explains the commands that produce the simulated return. It fits within Chapter 8 rather than requiring a separate top-level chapter.

1. **Ascent Commands** — pitch program, attitude commands, engine commands, and shutdown envelope.
2. **Return Commands** — boostback prediction and steering, entry-burn logic, landing-ignition prediction and deceleration, fuel budgets, and restart conditions.
3. **Force and Moment Allocation and Actuator Representation** — engine selection, throttle, requested force direction, virtual/physical control moments, grid-fin commands, and the limits of the implemented actuator representation.
4. **Integration, Events, and Sensor Scheduling** — RK4, accepted-state synchronization, scheduled/state-triggered events, contact refinement, guidance clocks, and sensor/output timing.
5. **Two-Vehicle Execution and Simulation Operation** — independent post-separation propagation, execution options, configuration/seed records, output files, and failure handling.

The aerodynamic response to a given fin deflection belongs in Chapter 5. The choice of fin deflection belongs in Chapter 8. Similarly, Chapter 3 defines thrust and flow for a supplied throttle; Chapter 8 explains how throttle is selected. Actuator descriptions should follow the implemented command-to-force/moment path and its actual assumptions and limits.

## Reusing Philipp's existing files

| Existing file | Placement in the combined report |
| --- | --- |
| [vehicle_mission.tex](content/vehicle_mission.tex) | Vehicle and propulsion information in Chapter 3; return planning allocations in Chapter 8. Merge overlapping tables with Kwanwoo's material. |
| [environment_models.tex](content/environment_models.tex) | Chapter 4. |
| [aerodynamics.tex](content/aerodynamics.tex) | Chapter 5. Keep aerodynamic constitutive equations here. |
| [ascent_staging.tex](content/ascent_staging.tex) | Separation/state-mapping principles in Chapter 3; ascent commands, shutdown execution, and upper-stage execution in Chapter 8. |
| [return_control.tex](content/return_control.tex) | Return algorithms and allocation in Chapter 8; resulting flight profiles and landing metrics in Chapter 9. |
| [sensors_navigation.tex](content/sensors_navigation.tex) | Split into Chapter 6 for measurement generation and Chapter 7 for state estimation. |
| [references.tex](content/references.tex) | Combine and deduplicate references from all authors. |

## Final integration and evidence

- Combine the vehicle/dynamics drafts before finalizing Chapter 3, so parameters, equations, and separation assumptions are stated once.
- Integrate the gyro-based attitude propagation, agree on sensor settings, and check the attitude/frame convention across Chapters 2, 6 and 7.
- Run the final integrated simulation and identify the execution settings behind every quantitative result. Keep reference-attitude diagnostic results distinguishable from gyro-driven navigation results.
- Use that evidence to finalize Chapters 9–11 together. Detailed parameter tables and long code excerpts can go in the appendix.
- Each author should retain source links and relevant page, table or figure numbers for the material they use, and distinguish published inputs, derived estimates and modeling assumptions.
