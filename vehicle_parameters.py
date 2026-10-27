"""Falcon-9-like research inputs; provenance is in data/vehicle_parameters.json.

K: official external geometry. E: literature estimates. A: declared model inputs.
This mixes source vehicles and is not a measured reconstruction of one mission.
"""
import math


class Falcon9Parameters:
    LBM_TO_KG = 0.45359237  # conversion
    G0 = 9.80665  # m/s², standard gravity

    def __init__(self):
        # base
        self.mass_ref_initial_nominal = 549000  # kg, reference only; not added to components
        self.mass = 0.0  # derived below from the selected components
        self.height =  70 # m
        self.body_diameter = 3.66  # m, official outer diameter
        self.radius = self.body_diameter / 2  # m, derived outer radius

        # Identity relative to the initial Body basis; no geographic orientation prescribed.
        self.quaternion_body_to_initial_body_wxyz = (1.0, 0.0, 0.0, 0.0)

        # geometry
        self.fairing_diameter = 5.2  # m, Df
        self.fairing_height = 13.2  # m, Lf
        self.stage1_length = 42.6  # m, L1, FT estimate
        self.stage2_length = 12.6  # m, L2, FT estimate
        self.body_area = math.pi * self.radius**2  # m², S
        self.fairing_area = math.pi * (self.fairing_diameter/2)**2  # m², Sf

        self.mass_ref_website = 549054  # kg, reference only
        self.stage1_dry_mass = 22200  # kg, md1
        self.stage1_propellant = 410900  # kg, mp1, initial
        self.stage2_dry_mass = 4000  # kg, md2
        self.stage2_propellant = 107500  # kg, mp2, initial
        self.fairing_mass = 1700  # kg, mf
        self.payload_mass = 10000  # kg, mPL, assumption
        self.total_propellant_faa = 1135925 * self.LBM_TO_KG  # kg, FAA Block 5
        self.stage1_propellant_sep = None  # kg, unresolved


        # A: geometry stations along +x from the vehicle base datum, not measured CGs.
        # RFZ Table 2 explicit dry groups; tank rows include liquid and are NOT dry masses.
        self.stage1_dry_groups = (
            ("shell", 8000.0, 21.0, 42.0),
            ("fins", 1520.0, 40.3, 0.0),
            ("engines", 4230.0, 0.452, 0.0),
            ("legs_folded", 2000.0, 6.66, 0.0),
        )  # name, mass kg, axial station m, equivalent axial length m
        self.stage1_dry_residual_mass = self.stage1_dry_mass - sum(g[1] for g in self.stage1_dry_groups)
        self.stage1_dry_residual_bottom = 4.0  # E: lower skirt, Ecker Table 1
        self.stage1_dry_residual_top = 4.0 + 14.7 + 22.1
        self.stage1_dry_residual_cg = (self.stage1_dry_residual_bottom + self.stage1_dry_residual_top)/2
        self.stage1_dry_cg = (
            sum(g[1]*g[2] for g in self.stage1_dry_groups)
            + self.stage1_dry_residual_mass*self.stage1_dry_residual_cg
        )/self.stage1_dry_mass  # E/A: sourced groups + uniform residual-hardware allocation
        self.stage2_cg_station = self.stage1_length+self.stage2_length/2  # m, A: consistent Stage offset + local midpoint
        self.fairing_cg_station = 63.4  # m, A
        self.payload_cg_station = 62.0  # m, A
        self.payload_length = 4.0  # m, A: equivalent cylinder
        self.payload_radius = 1.2  # m, A

        # Evidence-led research surrogate, NOT measured Falcon 9 internal geometry.
        self.stage1_oxidizer_fuel_ratio = 2.35  # E: RFZ 2023 sec.2.2 engine flow ratio
        self.stage1_rp1_initial = 123500.0  # E: Romero-Calvo 2023 Table 1
        self.stage1_lox_initial = 287400.0  # E: same inventory; sum = 410900 kg
        self.stage1_rp1_tank_skin = 0.0031  # m, E: Ecker 2017 Table 1
        self.stage1_lox_tank_skin = 0.0029  # m, E: same structural study
        self.stage1_rp1_tank_inner_radius = self.radius-self.stage1_rp1_tank_skin
        self.stage1_lox_tank_inner_radius = self.radius-self.stage1_lox_tank_skin
        self.stage1_rp1_temperature = 266.0  # K, E: Miquel 2020 Table 5
        self.stage1_lox_temperature = 66.0  # K, E: same comparative study
        self.stage1_rp1_density = 820.0  # kg/m³ at above T, E; pressure unspecified
        self.stage1_lox_density = 1250.0  # kg/m³ at above T, E; pressure unspecified
        self.stage1_rp1_tank_bottom = 4.0  # E/A: lower skirt length, equivalent flat bottom
        self.stage1_rp1_tank_length = 14.7  # E: Ecker 2017 Table 1 structural section
        self.stage1_lox_tank_length = 22.1  # E: same; domes not resolved
        self.stage1_lox_tank_bottom = self.stage1_rp1_tank_bottom+self.stage1_rp1_tank_length
        self.mass_model_name = "dual_tank"  # evidence profile; explicit balanced exhaust required

        self.mass = (self.stage1_dry_mass + self.stage1_propellant
                     + self.stage2_dry_mass + self.stage2_propellant
                     + self.fairing_mass + self.payload_mass)

        # S2 Table 2.16: use one literature engine set with the mass budget.
        self.stage1_engine_count = 9
        self.stage1_engine_thrust_vac = 914e3
        self.stage1_isp_vac = 311.0
        self.stage1_engine_exit_area = 0.68
        self.reference_sea_level_pressure_pa = 101325.0
        self.stage1_engine_thrust_sl = (self.stage1_engine_thrust_vac
                                      - self.stage1_engine_exit_area*self.reference_sea_level_pressure_pa)
        self.stage1_isp_sl = self.stage1_isp_vac*self.stage1_engine_thrust_sl/self.stage1_engine_thrust_vac
        self.stage1_engine_mass_flow_max = self.stage1_engine_thrust_vac/(self.stage1_isp_vac*self.G0)
        self.stage1_mass_flow_max = self.stage1_engine_count*self.stage1_engine_mass_flow_max
        self.stage1_thrust_max_sl = self.stage1_engine_count*self.stage1_engine_thrust_sl
        self.stage1_gimbal_max = math.radians(5.0)  # A, not validated actuator limit
        self.stage2_engine_count = 1
        self.stage2_engine_thrust_vac = 935e3
        self.stage2_isp_vac = 348.0
        self.stage2_engine_exit_area = 4.9
        self.stage2_mass_flow_max = self.stage2_engine_thrust_vac/(self.stage2_isp_vac*self.G0)
        self.stage2_gimbal_max = math.radians(5.0)  # A
        self.stage1_throttle_min = 0.60  # E, S2 Table 2.16
        self.stage2_throttle_min = 0.40  # E, S2 Table 2.16
        self.stage1_engine_dry_mass = 470.0  # E, S4 Table 2
        self.stage2_engine_dry_mass = 700.0  # E, S4 research-engine proxy
        self.stage1_engine_ring_radius = 1.3  # A, retained team layout
        self.stage1_engine_pivot_axial = 1.0  # A, pivot is not engine CG
        self.stage2_engine_pivot_axial = 0.5  # A
