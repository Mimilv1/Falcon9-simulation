from dataclasses import dataclass
from pathlib import Path
import json
import math


class EventConditionError(ValueError):
    pass


@dataclass(frozen=True)
class Transition:
    name: str
    phase: str
    time_s: float | None
    guards: tuple[str, ...]


def load_mission(path=None):
    path = Path(path) if path else Path(__file__).parent / "data" / "crs10_mission.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    duration = data["shutdown"]["duration_s"]
    cutoff = next(e["time_s"] for e in data["transitions"] if e["name"] == "meco")
    if not math.isfinite(duration) or not 0 < duration < cutoff:
        raise ValueError("Shutdown duration must be positive and shorter than ascent.")
    next(e for e in data["transitions"] if e["name"] == "shutdown_start")["time_s"] = (
        cutoff - duration
    )
    times = [e["time_s"] for e in data["transitions"] if e["time_s"] is not None]
    if any(not math.isfinite(t) or t < 0 for t in times) or times != sorted(set(times)):
        raise ValueError("Mission times must be finite, nonnegative and ordered.")
    names = [e["name"] for e in data["transitions"]]
    if len(names) != len(set(names)):
        raise ValueError("Mission event names must be unique.")
    return data


class MissionSequence:
    def __init__(self, plan=None):
        self.plan = load_mission() if plan is None else plan
        self.transitions = tuple(
            Transition(e["name"], e["phase"], e["time_s"], tuple(e["guards"]))
            for e in self.plan["transitions"]
        )
        self.index = 0
        self.phase = self.plan.get("initial_phase", "powered_ascent")
        self.log = [
            {
                "event": (
                    "liftoff" if self.phase == "powered_ascent" else "scenario_start"
                ),
                "time_s": 0.0,
                "phase": self.phase,
            }
        ]

    @property
    def next_event(self):
        return (
            self.transitions[self.index] if self.index < len(self.transitions) else None
        )

    @property
    def next_clock_time(self):
        return next(
            (e.time_s for e in self.transitions[self.index :] if e.time_s is not None),
            None,
        )

    def due(self, time_s, signals=None):
        event = self.next_event
        if event is None:
            return None
        if event.time_s is None:
            key = "ground_contact" if event.name == "touchdown" else "burn_complete"
            return event if (signals or {}).get(key, False) else None
        return event if time_s >= event.time_s - 1e-9 else None

    def accept(self, name, time_s, conditions):
        event = self.next_event
        if event is None or event.name != name:
            raise EventConditionError("An event cannot skip the preceding phase.")
        if not math.isfinite(time_s) or time_s < self.log[-1]["time_s"]:
            raise EventConditionError("Event time must be finite and nondecreasing.")
        if event.time_s is not None and time_s < event.time_s - 1e-9:
            raise EventConditionError("The scheduled event is not due.")
        failed = [key for key in event.guards if not conditions.get(key, False)]
        if failed:
            raise EventConditionError("Unmet event conditions: " + ", ".join(failed))
        phase = event.phase
        if name == "touchdown":
            phase = (
                "landed"
                if conditions.get("landing_success", False)
                else "landing_failed"
            )
        record = {
            "event": name,
            "time_s": float(time_s),
            "reference_time_s": event.time_s,
            "phase": phase,
            "conditions": {key: bool(conditions[key]) for key in event.guards},
        }
        self.phase = phase
        self.index += 1
        self.log.append(record)
        return record


def restart_conditions(
    vehicle,
    state,
    loads,
    *,
    required_propellant_kg=None,
    attitude_ready=False,
    propulsion_ready=False,
    descending=False
):
    import numpy as np

    lower = vehicle.mass_model().propellant_bounds(None)[0][0]
    margin = state.propellant_kg[0] - lower
    budget = required_propellant_kg
    return {
        "detached": len(vehicle.stages) == 1 and vehicle.stages[0].name == "stage1",
        "zero_flow": bool(
            np.all(loads.consumption_kgps == 0)
            and np.all(loads.consumption_rate_kgps2 == 0)
        ),
        "engines_off": all(
            not e.running and e.throttle == 0 for s in vehicle.stages for e in s.engines
        ),
        "fuel_budget_available": bool(
            budget is not None
            and np.isfinite(budget)
            and budget > 0
            and margin >= budget
        ),
        "attitude_ready": bool(attitude_ready),
        "propulsion_ready": bool(propulsion_ready),
        "descending": bool(descending),
    }


def landing_verdict(metrics, criteria):
    checks = {"ground_contact": bool(metrics.get("ground_contact", False))}
    for key, limit in criteria.items():
        value = metrics.get(key)
        checks[key] = bool(
            value is not None and math.isfinite(value) and 0 <= value <= limit
        )
    checks["legs_deployed"] = bool(metrics.get("legs_deployed", False))
    checks["landing_success"] = all(checks.values())
    return checks
