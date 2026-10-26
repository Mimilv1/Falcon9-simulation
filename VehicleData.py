from dataclasses import dataclass, field
import numpy as np
from copy import deepcopy


def quat_rotate(q, v):  # body to ENU
    w, qv = q[0], q[1:]
    t = 2.0 * np.cross(qv, v)
    return v + w * t + np.cross(qv, t)


@dataclass
class EngineSpec:  # shared, constant
    thrust_vac: float  # N
    isp_vac: float  # s
    isp_sl: float  # s
    exit_area: float = None  # m², unresolved
    throttle_min: float = None  # 0-1, unresolved
    gimbal_max: float = 0.0873  # rad, assumption 5°
    dry_mass: float = None  # kg, unresolved


@dataclass
class Engine:
    spec: EngineSpec
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))  # m, stage frame, unresolved
    running: bool = False  # on/off
    throttle: float = 0.0  # 0-1
    gimbal_pitch: float = 0.0  # rad
    gimbal_yaw: float = 0.0  # rad
    ignition_count: int = 0  # restarts


@dataclass
class Stage:
    dry_mass: float  # kg
    propellant_capacity: float  # kg, initial
    propellant_mass: float  # kg, remaining
    length: float  # m
    diameter: float  # m
    cg_dry: float = None  # m, from base, unresolved
    inertia_dry: np.ndarray = None  # kg·m², 3x3, unresolved
    engines: list = field(default_factory=list)  # Engine
    offset: float = 0.0  # m, base position in vehicle frame


@dataclass
class Vehicle:
    stages: list = field(default_factory=list)  # bottom to top
    fairing_mass: float = 0.0  # kg
    payload_mass: float = 0.0  # kg
    fairing_attached: bool = True  # flag

    # 6DOF state
    t: float = 0.0  # s
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))  # m, ENU
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))  # m/s, ENU
    quaternion: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.0, 0.0, 0.0]))  # body to ENU
    omega: np.ndarray = field(default_factory=lambda: np.zeros(3))  # rad/s, body frame

    def separate(self, index=1):  # lower = stages[:index], upper = stages[index:]
        if not 0 < index < len(self.stages):
            raise ValueError("index out of range")

        cut = self.stages[index].offset  # m, axial, vehicle frame
        r_body = np.array([cut, 0.0, 0.0])  # m, new origin, body frame

        lower_stages = deepcopy(self.stages[:index])
        upper_stages = deepcopy(self.stages[index:])
        for stage in upper_stages:
            stage.offset -= cut  # new origin

        lower = Vehicle(  # booster
            stages=lower_stages,
            fairing_mass=0.0,
            payload_mass=0.0,
            fairing_attached=False,
            t=self.t,
            position=self.position.copy(),
            velocity=self.velocity.copy(),
            quaternion=self.quaternion.copy(),
            omega=self.omega.copy(),
        )
        upper = Vehicle(  # upper stack
            stages=upper_stages,
            fairing_mass=self.fairing_mass,
            payload_mass=self.payload_mass,
            fairing_attached=self.fairing_attached,
            t=self.t,
            position=self.position + quat_rotate(self.quaternion, r_body),
            velocity=self.velocity + quat_rotate(self.quaternion, np.cross(self.omega, r_body)),  # v + w x r
            quaternion=self.quaternion.copy(),
            omega=self.omega.copy(),
        )
        return lower, upper

    def getmass(self):  # kg, current
        stages_mass = sum(s.dry_mass + s.propellant_mass for s in self.stages)  # dry + remaining propellant
        fairing = self.fairing_mass if self.fairing_attached else 0.0  # fairing dropped or not
        return stages_mass + fairing + self.payload_mass

# engine specs, shared
merlin_1d = EngineSpec(
    thrust_vac=914e3,  # N
    isp_vac=311,  # s
    exit_area=0.665,  # m², 0.92 m diameter
    throttle_min=0.57,  # 0-1
    gimbal_max=np.radians(5.0),  # rad, assumption
    dry_mass=470,  # kg
)
merlin_vac = EngineSpec(
    thrust_vac=981e3,  # N, alt 934e3
    isp_vac=348,  # s
    exit_area=8.55,  # m², 3.3 m diameter, alt 4.90
    throttle_min=0.64,  # 0-1
    gimbal_max=np.radians(5.0),  # rad, assumption
    dry_mass=500,  # kg, estimate
)

# stage 1 engines, octaweb
ring_radius = 1.3  # m, estimate
pivot_x = 1.0  # m, axial from base, estimate
stage1_engines = [Engine(spec=merlin_1d, position=np.array([pivot_x, 0.0, 0.0]))]  # center
stage1_engines += [
    Engine(spec=merlin_1d, position=np.array([pivot_x, ring_radius * np.cos(a), ring_radius * np.sin(a)]))
    for a in np.radians(np.arange(0, 360, 45))  # ring of 8
]

# stage 1
stage1 = Stage(
    dry_mass=22200,  # kg
    propellant_capacity=410900,  # kg
    propellant_mass=410900,  # kg, full
    length=42.6,  # m
    diameter=3.66,  # m
    cg_dry=21.3,  # m, uniform cylinder
    inertia_dry=np.diag([37172.79, 3375892.395, 3375892.395]),  # kg·m², roll, pitch, yaw
    engines=stage1_engines,
    offset=0.0,  # m, base
)

# stage 2
stage2 = Stage(
    dry_mass=4000,  # kg
    propellant_capacity=107500,  # kg
    propellant_mass=107500,  # kg, full
    length=12.6,  # m
    diameter=3.66,  # m
    cg_dry=6.3,  # m, uniform cylinder, estimate
    inertia_dry=np.diag([6697.8, 56268.9, 56268.9]),  # kg·m², roll, pitch, yaw, estimate
    engines=[Engine(spec=merlin_vac, position=np.array([0.5, 0.0, 0.0]))],  # estimate
    offset=42.6,  # m, on top of stage 1, no interstage
)

# vehicle
falcon9 = Vehicle(
    stages=[stage1, stage2],
    fairing_mass=1700,  # kg
    payload_mass=10000,  # kg, assumption
)

# For separation stage1, stage2 = falcon9.separate()