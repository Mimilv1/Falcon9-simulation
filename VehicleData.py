import math


class Rocket():
    LBM_TO_KG = 0.45359237  # conversion
    G0 = 9.80665  # m/s², standard gravity

    def __init__(self):
        # base
        self.mass = 549000  # kg, current total
        self.height =  70 # m
        self.radius = 1.85  # m

        # position (ENU)
        self.x = 0.0  # 
        self.y = 0.0  # 
        self.z = 0.0  # up

        # attitude
        self.phi = 0.0  # roll, x axis
        self.theta = 0.0  # pitch, y axis
        self.psi = 0.0  # yaw, z axis

        # velocity
        self.vx = 0.0  #  m/s
        self.vy = 0.0  #  m/s
        self.vz = 0.0  # up, m/s

        # angular rate
        self.omega_x = 0.0  # rad/s
        self.omega_y = 0.0  # rad/s
        self.omega_z = 0.0  # rad/s

        # time
        self.t = 0.0  # s, start

        # geometry
        self.body_diameter = 3.66  # m, D
        self.fairing_diameter = 5.2  # m, Df
        self.fairing_height = 13.2  # m, Lf
        self.stage1_length = 42.6  # m, L1, FT estimate
        self.stage2_length = 12.6  # m, L2, FT estimate
        self.body_area = 10.520879637606859  # m², S
        self.fairing_area = 21.237166338267002  # m², Sf

        self.mass_ref_website = 549054  # kg, reference only
        self.stage1_dry_mass = 22200  # kg, md1
        self.stage1_propellant = 410900  # kg, mp1, initial
        self.stage2_dry_mass = 4000  # kg, md2
        self.stage2_propellant = 107500  # kg, mp2, initial
        self.fairing_mass = 1700  # kg, mf
        self.payload_mass = 10000  # kg, mPL, assumption
        self.total_propellant_faa = 1135925 * self.LBM_TO_KG  # kg, FAA Block 5
        self.stage1_propellant_sep = None  # kg, unresolved

        # mass, current
        self.stage1_propellant_current = self.stage1_propellant  # kg, remaining
        self.stage2_propellant_current = self.stage2_propellant  # kg, remaining
        self.stage1_mass_current = self.stage1_dry_mass + self.stage1_propellant  # kg
        self.stage2_mass_current = self.stage2_dry_mass + self.stage2_propellant  # kg
        self.propellant_fraction_stage1 = 1.0  # remaining, 0-1
        self.propellant_fraction_stage2 = 1.0  # remaining, 0-1
        self.mass_flow_stage1 = 0.0  # kg/s
        self.mass_flow_stage2 = 0.0  # kg/s
        self.stage1_attached = True  # separation flag
        self.fairing_attached = True  # separation flag

        # engines, stage 1 (9x Merlin 1D)
        self.stage1_engine_count = 9  # Merlin 1D
        self.stage1_engine_thrust_sl = 845e3  # N, per engine, sea level
        self.stage1_engine_thrust_vac = 914e3  # N, per engine, vacuum (8227 kN / 9)
        self.stage1_isp_sl = 282  # s, sea level
        self.stage1_isp_vac = 311  # s, vacuum
        self.stage1_engine_mass_flow_max = self.stage1_engine_thrust_sl / (self.stage1_isp_sl * self.G0)  # kg/s, per engine, full throttle
        self.stage1_mass_flow_max = self.stage1_engine_count * self.stage1_engine_mass_flow_max  # kg/s, total, full throttle
        self.stage1_thrust_max_sl = self.stage1_engine_count * self.stage1_engine_thrust_sl  # N, total, sea level
        self.stage1_gimbal_max = math.radians(5.0)  # rad, max deflection, assumption, no source

        # engines, stage 2 (1x Merlin 1D Vacuum)
        self.stage2_engine_count = 1  # MVac
        self.stage2_engine_thrust_vac = 981e3  # N, vacuum
        self.stage2_isp_vac = 348  # s, vacuum
        self.stage2_mass_flow_max = self.stage2_engine_thrust_vac / (self.stage2_isp_vac * self.G0)  # kg/s, full throttle
        self.stage2_gimbal_max = math.radians(5.0)  # rad, max deflection, assumption, no source

        # engines, current
        self.stage1_engines_active = self.stage1_engine_count  # count
        self.stage1_throttle = 0.0  # 0-1
        self.stage2_throttle = 0.0  # 0-1
        self.stage1_thrust_current = 0.0  # N
        self.stage2_thrust_current = 0.0  # N
        self.stage1_gimbal_pitch = 0.0  # rad
        self.stage1_gimbal_yaw = 0.0  # rad
        self.stage2_gimbal_pitch = 0.0  # rad
        self.stage2_gimbal_yaw = 0.0  # rad

        # center of gravity, initial
        self.cg_fraction = 0.5  # uniform cylinder
        self.stage1_cg = 21.3  # m, from base
        self.cg_full_stack = None  # m, unresolved

        # center of gravity, current
        self.stage1_cg_current = self.stage1_cg  # m, from base
        self.cg_full_stack_current = None  # m, from base

        # inertia, initial, stage 1, uniform cylinder
        self.Ixx = 725204.295  # kg·m², roll, full
        self.Iyy = 65860315.14750001  # kg·m², pitch, full
        self.Izz = 65860315.14750001  # kg·m², yaw, full
        self.Ixx_dry = 37172.79000000001  # kg·m², roll, dry
        self.Iyy_dry = 3375892.395000001  # kg·m², pitch, dry
        self.Izz_dry = 3375892.395000001  # kg·m², yaw, dry
        self.Ixy = 0.0  # product of inertia
        self.Ixz = 0.0  # product of inertia
        self.Iyz = 0.0  # product of inertia
        self.inertia_full_stack = None  # kg·m², unresolved

        # inertia, current
        self.Ixx_current = self.Ixx  # kg·m², roll
        self.Iyy_current = self.Iyy  # kg·m², pitch
        self.Izz_current = self.Izz  # kg·m², yaw
        self.Ixy_current = 0.0  # product of inertia
        self.Ixz_current = 0.0  # product of inertia
        self.Iyz_current = 0.0  # product of inertia
        self.inertia_full_stack_current = None  # kg·m², unresolved

        # options
        self.slosh_enabled = False  # excluded


Faclon9 = Rocket()