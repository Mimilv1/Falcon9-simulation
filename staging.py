"""Zero-impulse separation map and conservation audit in common Body axes."""
from copy import deepcopy
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SeparationResult:
    lower: object
    upper: object
    lower_state: object
    upper_state: object
    audit: dict


def separate_coasting_vehicle(vehicle,state,time_s,loads, *, index=1,tolerance=1e-10):
    """Split an accepted zero-flow state; do not invent a separation kick.

    Engine commands must already be off. Positions and velocities of structural
    points are preserved; different child CG points naturally have different
    positions and rotation-dependent velocities. No discarded interstage mass.
    """
    if np.any(loads.consumption_kgps!=0) or np.any(loads.consumption_rate_kgps2!=0):
        raise ValueError('Separation requires zero mass flow and zero flow derivative.')
    if any(e.running or e.throttle!=0 for s in vehicle.stages for e in s.engines):
        raise ValueError('Turn engine commands off before zero-impulse separation.')
    if not np.isfinite(time_s) or not np.isfinite(tolerance) or tolerance<=0:
        raise ValueError('Require finite event time and positive audit tolerance.')
    parent=deepcopy(vehicle)
    parent.apply_state(state,time_s)
    p=parent.get_mass_properties()
    lower,upper=parent.separate(index)
    children=(lower,upper)
    states=tuple(child.to_state() for child in children)
    props=tuple(child.get_mass_properties() for child in children)
    attitude_error=max(min(np.linalg.norm(s.quaternion_body_to_initial_body_wxyz-state.quaternion_body_to_initial_body_wxyz),
                           np.linalg.norm(s.quaternion_body_to_initial_body_wxyz+state.quaternion_body_to_initial_body_wxyz)) for s in states)
    omega_error=max(np.linalg.norm(s.omega_body_radps-state.omega_body_radps) for s in states)
    if max(attitude_error,omega_error)>1e-12:
        raise ValueError('No-impulse separation must preserve attitude and angular rate.')
    mass=sum(c.mass_kg for c in props)
    cg=sum(c.mass_kg*s.position_body_m for c,s in zip(props,states))/mass
    momentum=sum(c.mass_kg*s.velocity_body_mps for c,s in zip(props,states))
    before_p=p.mass_kg*state.velocity_body_mps
    # Use the parent's instantaneous CG and translating velocity as the common
    # angular-momentum reference. Orbital child-CG terms must not be dropped.
    before_h=p.inertia_body_kgm2@state.omega_body_radps
    after_h=sum(c.inertia_body_kgm2@s.omega_body_radps+
                np.cross(s.position_body_m-state.position_body_m,
                         c.mass_kg*(s.velocity_body_mps-state.velocity_body_mps))
                for c,s in zip(props,states))
    kinetic_before=(.5*p.mass_kg*np.dot(state.velocity_body_mps,state.velocity_body_mps)
                    +.5*state.omega_body_radps@p.inertia_body_kgm2@state.omega_body_radps)
    kinetic_after=sum(.5*c.mass_kg*np.dot(s.velocity_body_mps,s.velocity_body_mps)
                       +.5*s.omega_body_radps@c.inertia_body_kgm2@s.omega_body_radps
                       for c,s in zip(props,states))
    cut=np.array([parent.stages[index].offset,0.,0.])
    position_error=max(np.linalg.norm(lower.position-parent.position),
                       np.linalg.norm(upper.position-parent.position-cut))
    velocity_error=max(np.linalg.norm(lower.velocity-parent.velocity),
                       np.linalg.norm(upper.velocity-parent.velocity-np.cross(parent.omega,cut)))
    mass_rel=abs(mass-p.mass_kg)/max(1.,p.mass_kg)
    linear_rel=float(np.linalg.norm(momentum-before_p)/max(1.,np.linalg.norm(before_p)))
    angular_rel=float(np.linalg.norm(after_h-before_h)/max(1.,np.linalg.norm(before_h)))
    energy_rel=abs(kinetic_after-kinetic_before)/max(1.,abs(kinetic_before))
    cg_error=float(np.linalg.norm(cg-state.position_body_m))
    if max(mass_rel,linear_rel,angular_rel,energy_rel)>tolerance or max(position_error,velocity_error,cg_error)>1e-8:
        raise ValueError('Separation state map failed its conservation audit.')
    audit={'time_s':float(time_s),'event':'stage_separation','model':'zero impulse; no discarded mass',
           'parent_mass_kg':p.mass_kg,'booster_mass_kg':props[0].mass_kg,'upper_stack_mass_kg':props[1].mass_kg,
           'mass_residual_kg':float(mass-p.mass_kg),
           'linear_momentum_residual_kgmps':float(np.linalg.norm(momentum-before_p)),
           'angular_momentum_residual_kgm2ps':float(np.linalg.norm(after_h-before_h)),
           'kinetic_energy_residual_j':float(kinetic_after-kinetic_before),
           'mass_relative_residual':float(mass_rel),'linear_relative_residual':linear_rel,
           'angular_relative_residual':angular_rel,'kinetic_energy_relative_residual':float(energy_rel),
           'recombined_cg_position_error_m':cg_error,
           'structural_position_continuity_error_m':float(position_error),
           'structural_velocity_continuity_error_mps':float(velocity_error),
           'quaternion_continuity_error':float(attitude_error),'angular_rate_continuity_error_radps':float(omega_error),
           'upper_base_original_station_m':upper.datum_shift_m,
           'note':'Child CG coordinates are different tracking points; shared fixed origin and initial attitude basis are retained.'}
    return SeparationResult(lower,upper,*states,audit)
