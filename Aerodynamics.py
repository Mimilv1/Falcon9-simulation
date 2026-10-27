"""First B aero model: drag only, sourced Cd(M), force through current CG."""
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np


class MachDragCurve:
    """Piecewise-linear figure readout; reject unmodeled Mach numbers."""
    def __init__(self,mach,coefficients):
        self.mach=np.asarray(mach,dtype=float)
        self.coefficients=np.asarray(coefficients,dtype=float)
        if (self.mach.ndim!=1 or len(self.mach)<2 or self.mach.shape!=self.coefficients.shape
            or not np.all(np.isfinite(self.mach)) or not np.all(np.isfinite(self.coefficients))
            or np.any(np.diff(self.mach)<=0) or self.mach[0]<0 or np.any(self.coefficients<0)):
            raise ValueError('Invalid Mach/drag-coefficient table.')

    @classmethod
    def from_json(cls,path,reference_area_m2):
        record=json.loads(Path(path).read_text(encoding='utf-8'))
        if not np.isclose(record['reference_area_m2'],reference_area_m2,rtol=0,atol=1e-9):
            raise ValueError('Cd curve and reference-area conventions must match.')
        return cls(record['mach'],record['drag_coefficient'])

    def __call__(self,mach):
        if not np.isfinite(mach) or mach<self.mach[0] or mach>self.mach[-1]:
            raise ValueError('Mach number is outside the adopted drag curve.')
        return float(np.interp(mach,self.mach,self.coefficients))


@dataclass(frozen=True)
class AeroOutput:
    force_body_n: np.ndarray
    relative_speed_mps: float
    mach: float
    dynamic_pressure_pa: float
    angle_of_attack_rad: float
    drag_coefficient: float


class DragOnlyAerodynamics:
    def __init__(self,reference_area_m2,drag_coefficient):
        if not np.isfinite(reference_area_m2) or reference_area_m2 <= 0:
            raise ValueError('Invalid reference area or drag coefficient.')
        if not callable(drag_coefficient) and (not np.isfinite(drag_coefficient) or drag_coefficient<0):
            raise ValueError('Invalid reference area or drag coefficient.')
        self.reference_area_m2=float(reference_area_m2)
        self.drag_coefficient=drag_coefficient

    def evaluate(self,state,properties,atmosphere,wind_body_mps):
        # Velocity of the fixed-body point coincident with CG at this instant,
        # not geometric-CG migration through the structure.
        migration=properties.cg_rate_body_mps
        relative=state.velocity_body_mps-migration-wind_body_mps
        speed=float(np.linalg.norm(relative))
        mach=speed/atmosphere.sound_speed_mps
        cd=float(self.drag_coefficient(mach) if callable(self.drag_coefficient) else self.drag_coefficient)
        if not np.isfinite(cd) or cd<0:
            raise ValueError('Drag coefficient must be finite and nonnegative.')
        dynamic_pressure=.5*atmosphere.density_kgm3*speed*speed
        drag=(-dynamic_pressure*self.reference_area_m2*cd*relative/speed
              if speed > 1e-10 else np.zeros(3))
        aoa=float(np.arctan2(np.linalg.norm(relative[1:]),relative[0])) if speed > 1e-10 else 0.0
        return AeroOutput(drag,speed,mach,dynamic_pressure,aoa,cd)
