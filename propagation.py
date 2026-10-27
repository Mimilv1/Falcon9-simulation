"""Reusable event-bounded propagation, without a flight scenario or thrust law."""
from dataclasses import dataclass
import numpy as np
from dynamics import State, rk4_step
from mass_flow import PropellantLimitError


@dataclass(frozen=True)
class PropagationResult:
    times_s: np.ndarray
    states: tuple[State, ...]
    stop_reason: str


def propagate(start_s,state,end_s,step_s,model,configuration,provider, *,
              breakpoints_s=(),inventory_tolerance_kg=1e-6,max_steps=1000000):
    """Integrate smooth load segments; stop before a limiting tank is overdrawn.

    A provider must supply continuous flow and its derivative within a segment.
    Split scheduled transitions at breakpoints; an abrupt flow jump requires an
    explicit impulse/state-map model and is outside this smooth integrator.
    Fuel-limit trials are shortened, never clipped into unphysical inventories.
    The returned inventory is within inventory_tolerance_kg of the first limit.
    """
    if (not np.all(np.isfinite([start_s,end_s,step_s,inventory_tolerance_kg]))
            or end_s <= start_s or step_s <= 0 or inventory_tolerance_kg <= 0):
        raise ValueError('Require a finite forward interval, positive step and tolerance.')
    boundaries = np.asarray(breakpoints_s,dtype=float)
    if boundaries.ndim != 1 or not np.all(np.isfinite(boundaries)):
        raise ValueError('Breakpoints must be a finite one-dimensional sequence.')
    boundaries = sorted({float(t) for t in boundaries if start_s < t < end_s})+[float(end_s)]
    time = float(start_s)
    times,states = [time],[state]
    reason = 'end_time'
    lower = np.asarray(model.propellant_bounds(configuration)[0]) if hasattr(model,'propellant_bounds') else np.zeros(2)
    if np.any(state.propellant_kg < lower-1e-7):
        raise PropellantLimitError('Initial inventory is below its model limit.')
    for _ in range(max_steps):
        loads = provider(time,state)
        active = loads.consumption_kgps > 0
        available = state.propellant_kg-lower
        if np.any(active & (available <= inventory_tolerance_kg)):
            reason = 'propellant_limit'
            break
        boundary = next(t for t in boundaries if t > time)
        h = min(step_s,boundary-time)
        if np.any(active):
            h = min(h,float(np.min(available[active]/loads.consumption_kgps[active])))
        for _retry in range(80):
            if time+h == time:
                raise RuntimeError('Cannot resolve the fuel boundary at this tolerance.')
            try:
                candidate = rk4_step(time,state,h,model,configuration,provider)
                if np.any(candidate.propellant_kg < lower-1e-7):
                    raise PropellantLimitError('Final state crossed an inventory boundary.')
                break
            except PropellantLimitError:
                h *= .5
        else:
            raise RuntimeError('Fuel-boundary subdivision did not converge.')
        time += h
        if abs(time-boundary) < 1e-12*max(1,abs(boundary)):
            time = boundary
        state = candidate
        times.append(time); states.append(state)
        if time >= end_s:
            break
    else:
        raise RuntimeError('Maximum propagation steps exceeded.')
    return PropagationResult(np.asarray(times),tuple(states),reason)
