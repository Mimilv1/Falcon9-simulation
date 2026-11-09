import math


class Falcon9Parameters:
    LBM_TO_KG = 0.45359237
    G0 = 9.80665

    def __init__(self):

        self.mass = 0.0
        self.body_diameter = 3.66
        self.radius = self.body_diameter / 2

        # Vertical on the pad: Body +x (nose) along ECEF +z at the launch site
        # (0,0,earth_radius); matches frames.LAUNCH_QUATERNION_BODY_TO_ECEF_WXYZ.
        _launch = math.sqrt(0.5)
        self.quaternion_body_to_ecef_wxyz = (_launch, 0.0, -_launch, 0.0)

        self.stage1_length = 42.6
        self.stage2_length = 12.6
        self.body_area = math.pi * self.radius**2

        self.stage1_dry_mass = 22200
        self.stage1_propellant = 410900
        self.stage2_dry_mass = 4000
        self.stage2_propellant = 107500
        self.payload_mass = 8430.0
        self.stage1_propellant_sep = None

        self.stage1_dry_groups = (
            ("shell", 8000.0, 21.0, 42.0),
            ("fins", 1520.0, 40.3, 0.0),
            ("engines", 4230.0, 0.452, 0.0),
            ("legs_folded", 2000.0, 6.66, 0.0),
        )
        self.stage1_dry_residual_mass = self.stage1_dry_mass - sum(
            g[1] for g in self.stage1_dry_groups
        )
        self.stage1_dry_residual_bottom = 4.0
        self.stage1_dry_residual_top = 4.0 + 14.7 + 22.1
        self.stage1_dry_residual_cg = (
            self.stage1_dry_residual_bottom + self.stage1_dry_residual_top
        ) / 2
        self.stage1_dry_cg = (
            sum(g[1] * g[2] for g in self.stage1_dry_groups)
            + self.stage1_dry_residual_mass * self.stage1_dry_residual_cg
        ) / self.stage1_dry_mass
        self.stage2_cg_station = self.stage1_length + self.stage2_length / 2
        self.payload_length = 7.2
        self.payload_radius = 1.83
        self.mission = "CRS-10"
        self.cargo_mass_kg = 2490.0
        self.payload_cg_station = (
            self.stage1_length + self.stage2_length + self.payload_length / 2
        )
        self.height = self.stage1_length + self.stage2_length + self.payload_length

        self.stage1_oxidizer_fuel_ratio = 2.35
        self.stage1_rp1_initial = 123500.0
        self.stage1_lox_initial = 287400.0
        self.stage1_rp1_tank_skin = 0.0031
        self.stage1_lox_tank_skin = 0.0029
        self.stage1_rp1_tank_inner_radius = self.radius - self.stage1_rp1_tank_skin
        self.stage1_lox_tank_inner_radius = self.radius - self.stage1_lox_tank_skin
        self.stage1_rp1_temperature = 266.0
        self.stage1_lox_temperature = 66.0
        self.stage1_rp1_density = 820.0
        self.stage1_lox_density = 1250.0
        self.stage1_rp1_tank_bottom = 4.0
        self.stage1_rp1_tank_length = 14.7
        self.stage1_lox_tank_length = 22.1
        self.stage1_lox_tank_bottom = (
            self.stage1_rp1_tank_bottom + self.stage1_rp1_tank_length
        )
        self.mass_model_name = "dual_tank"

        self.mass = (
            self.stage1_dry_mass
            + self.stage1_propellant
            + self.stage2_dry_mass
            + self.stage2_propellant
            + self.payload_mass
        )

        self.stage1_engine_count = 9
        self.stage1_engine_thrust_vac = 756000.0 + 0.68 * 101325.0
        self.stage1_isp_vac = 311.0
        self.stage1_engine_exit_area = 0.68
        self.reference_sea_level_pressure_pa = 101325.0
        self.stage1_engine_thrust_sl = (
            self.stage1_engine_thrust_vac
            - self.stage1_engine_exit_area * self.reference_sea_level_pressure_pa
        )
        self.stage1_isp_sl = (
            self.stage1_isp_vac
            * self.stage1_engine_thrust_sl
            / self.stage1_engine_thrust_vac
        )
        self.stage1_engine_mass_flow_max = self.stage1_engine_thrust_vac / (
            self.stage1_isp_vac * self.G0
        )
        self.stage1_mass_flow_max = (
            self.stage1_engine_count * self.stage1_engine_mass_flow_max
        )
        self.stage1_thrust_max_sl = (
            self.stage1_engine_count * self.stage1_engine_thrust_sl
        )
        self.stage1_gimbal_max = math.radians(5.0)
        self.stage2_engine_count = 1
        self.stage2_engine_thrust_vac = 935e3
        self.stage2_isp_vac = 348.0
        self.stage2_engine_exit_area = 4.9
        self.stage2_mass_flow_max = self.stage2_engine_thrust_vac / (
            self.stage2_isp_vac * self.G0
        )
        self.stage2_gimbal_max = math.radians(5.0)
        self.stage1_throttle_min = 0.60
        self.stage2_throttle_min = 0.40
        self.stage1_engine_dry_mass = 470.0
        self.stage2_engine_dry_mass = 700.0
        self.stage1_engine_ring_radius = 1.3
        self.stage1_engine_pivot_axial = 1.0
        self.stage2_engine_pivot_axial = 0.5
