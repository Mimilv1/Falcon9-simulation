# Documentation audit status

Reviewed 10 November 2026 against the final nominal run, runtime revision
`2d89a886a761407695eb26405aca78a9af6663a6`, and report reproduction ledger.
Code and active configuration define behavior; old prose is not an override.

## Implemented limits documented in the report

- **FL-2 — Ideal actuation:** `Actuators.py` is empty and nonzero engine
  gimbal commands are rejected. Ascent uses an unlimited ideal couple;
  return uses bounded ideal residual moments and static fin allocation.
  The report does not claim physical TVC/servo performance or navigation-fed
  landing control. Physical actuator dynamics remain a later GNC extension.
- **FL-3 — Upper-stage aerodynamics:** `second_stack` is a declared/tested
  configuration, while `upper_stage.py` actively uses its separate constant
  CD = 0.3 provider. Report Chapters 5/8 and `aero_model.md` state this choice.
- **FL-8 — Heading field:** the current ascent target uses projected launch
  east, rather than the stored 45-degree heading field. Chapter 8 states the
  implemented direction and does not claim reconstruction of the ISS plane.

## Resolved checks

| ID | Finding | Final disposition |
| --- | --- | --- |
| FL-0 | Environment prose and dead references | Central gravity, separate rotating-frame terms and altitude-limited wind are documented. |
| FL-1 | Aero configuration count and deployment states | `aero_model.md` now lists all five configurations, active fin commands and upper-stage provider selection. |
| FL-4 | Stale supplied-truth-attitude metadata | Mission configuration identifies the gyro-driven MEKF; active model includes 30 error coordinates. |
| FL-5 | Frame, notation and result placement | Body +x, launch-aligned Earth-centred frame and fixed-IMU comparison are shared; results remain in Chapter 9. |
| FL-6 | Rejected aiding and covariance mismatch | Accel-scale and GNSS colored-bias states, midpoint/SO(3) derivatives and epoch/history corrections are active. Final R1 accepts 381/381 valid star fixes at the unchanged 16.27 gate. Full covariance and innovations are saved. A separate 20-trial position ANEES check passes; full-flight ensemble consistency is not claimed. |
| FL-7 | Print-only navigation diagnostics | Fourteen assertion-based navigation tests and ten stored independent validation checks now cover selected derivatives, nonzero calibration, delay/startup, gating, covariance and RK4 order. The legacy print-only script is not counted as certification. |
| FL-9 | Metadata scope and cargo source | Cargo points to the official overview; active rotating-frame/two-child/estimated-attitude scope strings are corrected. |
| FL-10 | Smart quotes in copied PDF code | Monospaced quote mapping is disabled. Extracted command uses ASCII quotes; the extracted ValueError line compiles. |
| FL-11 | Missing raw-data/plot reproduction path | Tracked `report/reproduction_data.zip`, `reproduce.py`, `validate_navigation.py` and member hashes provide a checked fresh-extraction path. All ten regenerated figure PDFs match their saved SHA-256 hashes. |
| FL-12 | Return gain units | Return/landing direction-to-rate and rate-to-acceleration gains both use s^-1. Direct ascent/upper PD kp remains s^-2. |
| FL-13 | Landing guidance interval | Main text and appendix give the allowed 0.01–0.2 s interval. |
| FL-14 | Command clipping described as state bounds | Desired rates are clipped per axis; acceleration cap is on the command norm. Actual rate history and its 17.825 deg/s coast maximum are shown. |

## Evidence boundary

Final nominal navigation has position/velocity/attitude RMS of
0.670763 m / 0.0337753 m/s / 0.207582 deg. Calibrated instantaneous gyro
rate RMS is 0.0125122 deg/s versus raw 0.0115190 deg/s for zero nominal
scale/cross-axis truth. The separate deliberately excited nonzero-error
fixture improves rate RMS from 0.0375088 to 0.0100380 deg/s. These are
separate claims, not a guarantee that calibration improves every sample.
Ascent reference data are also used for calibration. Calm nominal landing,
ideal actuation and independent position ANEES do not establish physical
return validation, wind robustness or full-flight covariance certification.
