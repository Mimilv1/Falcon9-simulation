# Review of team update 5b61a45

## Findings and resulting changes

| Finding in upstream code | Resolution |
| --- | --- |
| EngineSpec requires isp_sl, but neither engine constructor supplies it; import fails | First-stage SL Isp is derived from the chosen literature engine set; MVac SL Isp is explicitly None. |
| ENU state fields coexist with Body-x engine/staging geometry and an identity quaternion | Retain ENU / Body+x; set initial axes x up, y east, z north with q=[0.5,−0.5,−0.5,−0.5]. |
| Prior local A code used NED / Body+z | Convert all four core files, state/output names, inertia ordering and analytical checks to ENU / Body+x. |
| Stage-1 dry CG is a midpoint assumption and inertia is a uniform solid-cylinder proxy | Use the sourced dry-component budget, with declared residual allocation and conditional CG bounds. |
| Engine thrust/exit areas/throttle limits combine unsupported or different-version numbers | Use Pezent Table 2.16 coherently with its Table 2.17 mass set; label engine geometry/actuator assumptions separately. |
| No documented tank-level evidence in Stage objects | Add immutable TankSpec references and a fixed-O/F liquid profile; derive inventories from the existing stage-total state. |
| Vehicle separation tracks a structural base, while the A state tracks CG | Preserve the team's base-origin separation, document it, add explicit fixed-inventory to_state mapping and trial-state mass adapter. |
| main.py advances time with its update call commented out | Preserve the team's entry file; document that it is not an implemented flight simulation. |

## Sharing scope

Shared: four reusable mathematical modules, team objects, evidence inputs/source
metadata, tests and English documentation. No A demonstration files, scenario
configuration, generated plots, trajectories, third-party PDF copies or runtime
folders are included. Synthetic test fixtures are explicitly diagnostic.

No powered dual-tank/open-system implementation is claimed. The active-depletion
guard is retained. B/C physics, navigation and flight-data validation remain
separate project work.
