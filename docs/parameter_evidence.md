# Adopted Vehicle Parameters and Supporting Evidence

K = authoritative known value; E = literature estimate; A = explicit modeling assumption.
E/A denotes a derived estimate that depends on a stated geometric or allocation assumption.
All axial distances follow **Body +x toward the nose**. All shared vectors are expressed in current **Body** axes; geographic conversion is deferred.

| Parameter | Adopted value | Class | Evidence and method |
| --- | --- | --- | --- |
| Body outer diameter | 3.66 m | K | [SpaceX guide, Table 2-1](https://www.spacex.com/assets/media/falcon-users-guide-2025-05-09.pdf) |
| Stage dry masses | Stage 1: 22,200 kg; Stage 2: 4,000 kg | E | [Pezent, Table 2.17](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) |
| Initial stage propellant | Stage 1: 410,900 kg; Stage 2: 107,500 kg | E | [Pezent, Table 2.17](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) |
| Fairing mass | 1,700 kg | E | [Pezent, Table 2.17](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content); separate from stage dry budgets |
| Initial RP-1 / LOX inventory | 123,500 kg / 287,400 kg | E | [Romero-Calvo, Table 1, p.215](https://hanspeterschaub.info/PapersPrivate/RomeroCalvo2023.pdf); sum = 410,900 kg |
| Loaded O/F mass ratio | 2.3271255 | E | Calculated from the initial inventory in [S5](https://hanspeterschaub.info/PapersPrivate/RomeroCalvo2023.pdf) |
| Consumption O/F mass ratio | 2.35 | E | [Bykerk, sec.2.2](https://www.eucass.eu/doi/EUCASS2023-018.pdf); constant research-engine operating model |
| Tank skin thickness | RP-1: 3.1 mm; LOX: 2.9 mm | E | [Ecker, Table 1](https://doi.org/10.2514/6.2017-4878); modeled structural thicknesses |
| Equivalent tank inner radius | RP-1: 1.8269 m; LOX: 1.8271 m | E/A | [S1](https://www.spacex.com/assets/media/falcon-users-guide-2025-05-09.pdf) outer radius minus [S3](https://doi.org/10.2514/6.2017-4878) modeled skin; internal hardware omitted |
| RP-1 density / temperature | 820 kg/m³ at 266 K | E | [Miquel, Table 5, PDF p.24](https://www.diva-portal.org/smash/get/diva2:1440090/FULLTEXT01.pdf); pressure unspecified, constant property |
| LOX density / temperature | 1,250 kg/m³ at 66 K | E | [Miquel, Table 5, PDF p.24](https://www.diva-portal.org/smash/get/diva2:1440090/FULLTEXT01.pdf); pressure unspecified, constant property |
| Equivalent tank section lengths | RP-1: 14.7 m; LOX: 22.1 m | E/A | [Ecker, Table 1](https://doi.org/10.2514/6.2017-4878) structural sections used as flat cylindrical envelopes |
| Tank bottom axial stations | RP-1: x=4.0 m; LOX: x=18.7 m | E/A | [S3](https://doi.org/10.2514/6.2017-4878) lower skirt + ordered section lengths; aligned source datums assumed |
| Stage 1 nominal dry CG | 17.52108 m above base | E/A | [Bykerk, Table 2](https://www.eucass.eu/doi/EUCASS2023-018.pdf) dry groups + [S2](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) total; 6,450 kg residual uniformly over x=4.0–40.8 m |
| Stage 1 dry inertia [Ixx,Iyy,Izz] | [61,368.59, 4,441,994.04, 4,441,994.04] kg·m² | E/A | Component shapes and parallel-axis calculation; Ixx is axial. Not measured inertia. |
| Vacuum engine thrust | Merlin 1D SL: 914 kN; MVac: 935 kN | E | [Pezent, Table 2.16](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content); per engine, selected literature set |
| Vacuum specific impulse | Merlin 1D SL: 311 s; MVac: 348 s | E | [Pezent, Table 2.16](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) |
| Nozzle exit area | Merlin 1D SL: 0.68 m²; MVac: 4.90 m² | E | [Pezent, Table 2.16](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) |
| Minimum throttle | Merlin 1D SL: 60%; MVac: 40% | E | [Pezent, Table 2.16](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) research model limits, not certified flight limits |
| Merlin 1D SL reference thrust / Isp | 845.099 kN / 287.556 s | E/A | [S2](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content) vacuum data with p=101,325 Pa: Tsl=Tvac−Ae·p, Isp_sl=Isp_vac·Tsl/Tvac |
| Engine dry-mass proxies | Merlin 1D: 470 kg; MVac: 700 kg | E | [Bykerk, Table 2](https://www.eucass.eu/doi/EUCASS2023-018.pdf); engine masses are already included in stage dry budgets |

## Retained assumptions and unresolved data

| Parameter | Adopted value | Class | Basis / limitation |
| --- | --- | --- | --- |
| Stage lengths / stage-2 base offset | 42.6 m / 12.6 m / 42.6 m | E/A | Legacy FT geometry plus the team's contiguous-stage model; actual interstage boundaries remain unresolved. |
| Stage-2 dry CG | 6.3 m from its own base; 48.9 m in the full vehicle | A | Uniform-cylinder midpoint and stage offset; not measured CG. |
| Fairing CG | 63.4 m from original base | A | Equivalent 13.2 m cylinder ending at 70 m. |
| Payload mass / CG / shape | 10,000 kg / 62.0 m / 4.0 m length, 1.2 m radius | A | Mission placeholder; replace after selecting a flight/payload. |
| Engine-ring radius | 1.3 m | A | Retained team layout; independent dimensional evidence unresolved. |
| Engine pivot axial stations | Stage 1: 1.0 m; stage 2: 0.5 m | A | Retained team force-application stations; these are not engine CGs. |
| Gimbal limit | 5° | A | Retained team assumption; verified actuator specification unavailable. |
| Residual dry-hardware placement | 6,450 kg uniformly over x=4.0–40.8 m | A | Remaining dry budget after the four RFZ dry groups; no surveyed location. |
| Liquid geometry | Settled flat cylinders, fixed density | A | Domes, stringers, transfer tube, slosh and reversed/free-fall liquid behavior omitted. |
| Earth model | Body-only states, Earth rotation omitted | A | Initial Body is the mathematical attitude reference; geographic outputs deferred. |

The dry-CG nominal value is not measured. If the remaining hardware's CG can lie
anywhere in x=4.0–40.8 m while all listed groups remain fixed, the conditional
dry-CG range is **12.17513–22.86703 m**. This is neither a confidence interval nor
a complete uncertainty bound on the literature estimates.

The selected consumption ratio depletes LOX before RP-1, leaving **1202.12766 kg**
RP-1 at the mathematical burn limit. This is a consequence of the adopted
inventory/ratio pair, not a known Falcon 9 landing reserve. Independent tank
states are required for variable O/F, venting or refilling.

The engine/mass set is from Pezent; geometry, dry-group distribution and fluid
conditions come from other research vehicles/versions. This is an explicitly
mixed-source Falcon-9-like model, not one identified Block 5 mission. RFZ full
tank rows include tank weight and were not used as pure propellant or empty-tank
masses. Source RFZ x-axial stations now match the agreed Body +x convention.

## Sources

The retained temporary drag model uses the approximate Mach table and source
locator in `data/drag_curve_pezent.json`, paired with reference area 21.08 m².
The temporary environment uses the dry 0–11 km troposphere from the
[US Standard Atmosphere 1976](https://ntrs.nasa.gov/citations/19770009539):
T₀ = 288.15 K, p₀ = 101,325 Pa, lapse = 0.0065 K/m,
R = 287.05287 J/(kg·K), γ = 1.4. Geometric height approximates geopotential
height (E/A); constant gravity magnitude and default zero wind are assumptions.
These component models do not define a team flight schedule.

- **S1** [SpaceX Falcon User Guide, May 2025](https://www.spacex.com/assets/media/falcon-users-guide-2025-05-09.pdf). Table 2-1, printed p.9. Official external geometry, not internal tank/CG data
- **S2** [Pezent dissertation, 2024](https://ir-api.ua.edu/api/core/bitstreams/cafcf02b-b6e7-47d8-a41b-24b1ba12a154/content). Tables 2.16-2.17, printed p.77 / PDF p.94. Consistent literature engine and mass set; not measured Block 5 flight data
- **S3** [Ecker et al., A Numerical Study on the Thermal Loads During a Supersonic Rocket Retro-Propulsion Maneuver, 2017](https://doi.org/10.2514/6.2017-4878). Table 1, paper p.3. Reverse-engineered Falcon 9 FT structural study, not manufacturer drawings
- **S4** [Tamas Bykerk, A standard model for the investigation of aerodynamic and aerothermal loads on a re-usable launch vehicle, 2023](https://www.eucass.eu/doi/EUCASS2023-018.pdf). sec.2.2 p.4; Table 2 p.5. Generic Falcon-9-like research vehicle; explicit estimates. Tank rows include tank weight; not pure liquids.
- **S5** [Romero-Calvo et al., Novel propellant settling strategies for liquid rocket engine restart in microgravity, 2023](https://hanspeterschaub.info/PapersPrivate/RomeroCalvo2023.pdf). Table 1, journal p.215 / PDF p.2. Public-literature capacities only. Its upper-stage mass and ambiguous CG/geometry are not adopted.
- **S6** [Valentin Miquel, Propellant Feeding System of a Liquid Rocket With Multiple Engines, KTH, 2020](https://www.diva-portal.org/smash/get/diva2:1440090/FULLTEXT01.pdf). Table 5, PDF p.24, Falcon 9 first-stage row. Comparative study: RP-1 266 K, 820 kg/m³; subcooled LOX 66 K, 1250 kg/m³. Pressure and original measurement provenance not provided; E reference conditions, not certified flight properties.
- **S7** [Huber, Lemmon and Bruno, Effect of RP-1 Compositional Variability on Thermophysical Properties, NIST, 2009](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=902380). Table 5, journal p.5553 / PDF p.4. Cross-check only: model density 808.1 kg/m³ at 288.7 K, 101.325 kPa. Not substituted for cold RP-1.
- **S8** Course Module 4 Simulation Examples. pp.10-13. NED/body state and rigid-body equations; numerical examples are not measured Falcon parameters
- **S9** [FAA Annual Compendium of Commercial Space Transportation, 2018](https://www.faa.gov/about/office_org/headquarters_offices/ast/media/2018_AST_Compendium.pdf). Falcon 9 FT fact sheet. Legacy FT geometry; stage/interstage boundary applicability remains unresolved. Not used for v3 internal tank section lengths.
