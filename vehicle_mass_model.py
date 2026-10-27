"""Evidence profile adapter for the reusable core; Body +x is longitudinal."""
from enum import Enum
import json
from pathlib import Path
import numpy as np
from vehicle_parameters import Falcon9Parameters
from mass_properties import (
    MassProperties, cylinder_inertia_per_kg, fixed_cylinder, liquid_column,
    combine_components,
)


class Configuration(str, Enum):
    FULL_STACK = "full_stack"
    FIRST_STAGE = "first_stage"


class VehicleModel:
    """Literature-input adapter for mass properties, independent of Vehicle objects.

    Geometric derivatives are coupled to explicit exhaust data by dynamics.py.
    """

    def __init__(self, parameter_path=None, *, inputs=None, mass_model=None):
        self.inputs = inputs if inputs is not None else Falcon9Parameters()
        self.mass_model_name = mass_model or self.inputs.mass_model_name
        if self.mass_model_name not in ("uniform_cylinder", "dual_tank"):
            raise ValueError("Choose uniform_cylinder or dual_tank.")
        path = parameter_path or Path(__file__).parent / "data" / "vehicle_parameters.json"
        with open(path, encoding="utf-8") as handle:
            self.source = json.load(handle)
        self.parameters = self.source["parameters"]
        self.value = {}
        for key, record in self.parameters.items():
            value = (getattr(self.inputs, record["attribute"])
                     if "attribute" in record else record["value"])
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be a finite positive value.")
            self.value[key] = float(value)
        if not np.isclose(self.value["stage1_rp1_initial_kg"]+self.value["stage1_lox_initial_kg"],
                          self.value["stage1_propellant_initial_kg"], rtol=0, atol=1e-6):
            raise ValueError("Initial tank inventories must sum to the Stage 1 total.")
        self.dry_groups = tuple(tuple(g) for g in self.inputs.stage1_dry_groups)
        self.dry_residual_mass = self.value["stage1_dry_mass_kg"]-sum(g[1] for g in self.dry_groups)
        self.dry_residual_bottom = float(self.inputs.stage1_dry_residual_bottom)
        self.dry_residual_top = float(self.inputs.stage1_dry_residual_top)
        if self.dry_residual_mass < 0 or self.dry_residual_top <= self.dry_residual_bottom:
            raise ValueError("Invalid dry-mass budget or residual-hardware envelope.")
        derived_cg = (sum(g[1]*g[2] for g in self.dry_groups)
                      + self.dry_residual_mass*(self.dry_residual_bottom+self.dry_residual_top)/2
                      )/self.value["stage1_dry_mass_kg"]
        if not np.isclose(derived_cg, self.value["stage1_dry_cg_station_m"], rtol=0, atol=1e-8):
            raise ValueError("Dry CG must be recalculated from the selected component budget.")

    cylinder_inertia_per_kg = staticmethod(cylinder_inertia_per_kg)

    def stage1_tank_geometry(self):
        """Fixed literature section boundaries; liquid heights follow inventory/density."""
        v = self.value
        g = {}
        for fluid in ("rp1", "lox"):
            radius = v[f"stage1_{fluid}_tank_inner_radius_m"]
            initial = v[f"stage1_{fluid}_initial_kg"]
            length = v[f"stage1_{fluid}_tank_length_m"]
            height = initial/(v[f"stage1_{fluid}_density_kgm3"]*np.pi*radius**2)
            if radius > v["body_diameter_m"]/2 or height > length:
                raise ValueError(f"{fluid} inventory exceeds the equivalent tank envelope.")
            g.update({f"{fluid}_radius_m": radius, f"{fluid}_initial_kg": initial,
                      f"{fluid}_bottom_x_m": v[f"stage1_{fluid}_tank_bottom_m"],
                      f"{fluid}_tank_length_m": length, f"{fluid}_full_height_m": height,
                      f"{fluid}_unfilled_height_m": length-height})
        # Compatibility plotting radius only; calculations use the two separate radii.
        g["radius_m"] = g["rp1_radius_m"]
        return g

    def stage1_burn_limits(self):
        v = self.value
        ratio = v["stage1_oxidizer_fuel_ratio"]
        burnable = min(v["stage1_rp1_initial_kg"]*(1+ratio),
                       v["stage1_lox_initial_kg"]*(1+ratio)/ratio)
        return {"burnable_kg": burnable,
                "minimum_remaining_kg": v["stage1_propellant_initial_kg"]-burnable,
                "loaded_oxidizer_fuel_ratio": v["stage1_lox_initial_kg"]/v["stage1_rp1_initial_kg"],
                "flow_oxidizer_fuel_ratio": ratio}

    def stage1_inventory(self, propellant_kg):
        """Reconstruct inventory for one constant-O/F burn, without vents/refilling.

        State retains its two-stage total-mass interface. Variable O/F, venting,
        transfer, or a different initial loading require independent tank states.
        """
        v = self.value
        if not np.isfinite(propellant_kg):
            raise ValueError("Propellant must be finite.")
        burned = v["stage1_propellant_initial_kg"]-float(propellant_kg)
        limit = self.stage1_burn_limits()
        if burned < -1e-7 or propellant_kg < limit["minimum_remaining_kg"]-1e-7:
            raise ValueError("Stage 1 inventory outside constant-O/F burn limits; stop at the limiting tank.")
        ratio = v["stage1_oxidizer_fuel_ratio"]
        return (max(0.0, v["stage1_rp1_initial_kg"]-burned/(1+ratio)),
                max(0.0, v["stage1_lox_initial_kg"]-burned*ratio/(1+ratio)))

    def stage1_liquid_components(self, propellant_kg, consumption_kgps=0.0,consumption_rate_kgps2=0.0):
        g, v = self.stage1_tank_geometry(), self.value
        masses = self.stage1_inventory(propellant_kg)
        ratio = v["stage1_oxidizer_fuel_ratio"]
        rates = (-consumption_kgps/(1+ratio), -consumption_kgps*ratio/(1+ratio))
        rate_changes = (-consumption_rate_kgps2/(1+ratio),-consumption_rate_kgps2*ratio/(1+ratio))
        return tuple(liquid_column(f"stage1_{fluid}", mass, rate, g[f"{fluid}_bottom_x_m"],
                                   g[f"{fluid}_radius_m"], v[f"stage1_{fluid}_density_kgm3"],
                                   g[f"{fluid}_tank_length_m"],mass_acceleration_kgps2=change)
                     for fluid, mass, rate,change in zip(("rp1", "lox"), masses, rates,rate_changes))

    def stage1_dry_components(self):
        """RFZ dry-group proxy with explicit residual allocation and axisymmetry A.

        Shell is a thin cylindrical shell; localized groups are equivalent thin
        disks inside the body envelope (folded legs). These local inertias are
        idealizations; positions/masses alone do not determine actual inertia.
        """
        r = self.value["body_diameter_m"]/2
        groups = [fixed_cylinder("stage1_dry_"+name, mass, 0.0, z, length, r)
                  for name, mass, z, length in self.dry_groups]
        shell = groups[0]
        from dataclasses import replace
        shell_shape = np.diag([r*r]+[(6*r*r+shell.length_m**2)/12]*2)
        groups[0] = replace(shell, inertia_at_center_kgm2=shell.mass_kg*shell_shape)
        bottom = self.dry_residual_bottom
        top = self.dry_residual_top
        residual = fixed_cylinder("stage1_dry_unallocated", self.dry_residual_mass,
                                  0.0, (bottom+top)/2, top-bottom, r)
        shape = np.diag([r*r]+[(6*r*r+(top-bottom)**2)/12]*2)
        groups.append(replace(residual, inertia_at_center_kgm2=residual.mass_kg*shape))
        return tuple(groups)

    def stage1_dry_uncertainty(self):
        """Conditional CG bounds for fixed sourced groups and residual within tank region."""
        first_moment = sum(g[1]*g[2] for g in self.dry_groups)
        residual = self.dry_residual_mass
        mass = self.value["stage1_dry_mass_kg"]
        return {"unallocated_dry_mass_kg": residual,
                "nominal_cg_m": self.value["stage1_dry_cg_station_m"],
                "conditional_lower_cg_m": (first_moment+residual*self.dry_residual_bottom)/mass,
                "conditional_upper_cg_m": (first_moment+residual*self.dry_residual_top)/mass}

    def propellant_bounds(self,configuration):
        lower = self.stage1_burn_limits()['minimum_remaining_kg'] if self.mass_model_name == 'dual_tank' else 0.0
        upper2 = self.value['stage2_propellant_initial_kg'] if Configuration(configuration) is Configuration.FULL_STACK else 0.0
        return np.array([lower,0.0]),np.array([self.value['stage1_propellant_initial_kg'],upper2])

    def mass_properties(self, configuration, propellant_kg, consumption_kgps=(0.0, 0.0),
                        consumption_rate_kgps2=(0.0,0.0)):
        configuration = Configuration(configuration)
        prop = np.asarray(propellant_kg, dtype=float)
        rate = np.asarray(consumption_kgps, dtype=float)
        rate_change = np.asarray(consumption_rate_kgps2,dtype=float)
        if prop.shape != (2,) or rate.shape != (2,) or rate_change.shape != (2,) or not np.all(np.isfinite(np.r_[prop, rate,rate_change])):
            raise ValueError("Propellant and consumption arrays need two finite components.")
        if np.any(prop < -1e-7) or np.any(rate < 0):
            raise ValueError("Propellant and consumption must be nonnegative.")
        prop = np.maximum(prop, 0.0)
        if configuration is Configuration.FIRST_STAGE and (prop[1] > 1e-7 or rate[1] > 0 or rate_change[1] != 0):
            raise ValueError("A detached booster cannot retain or consume Stage 2 propellant.")
        v = self.value
        radius = v["body_diameter_m"]/2
        components = (list(self.stage1_dry_components()) if self.mass_model_name == "dual_tank"
                      else [fixed_cylinder("stage1_dry", v["stage1_dry_mass_kg"], 0.0,
                            v["stage1_length_m"]/2, v["stage1_length_m"], radius)])
        if self.mass_model_name == "dual_tank":
            components.extend(self.stage1_liquid_components(prop[0], rate[0],rate_change[0]))
        else:
            components.append(fixed_cylinder("stage1_propellant", prop[0], -rate[0],
                              v["stage1_length_m"]/2, v["stage1_length_m"], radius,
                              mass_acceleration_kgps2=-rate_change[0]))
        if configuration is Configuration.FULL_STACK:
            components += [
                fixed_cylinder("stage2_dry", v["stage2_dry_mass_kg"], 0.0,
                               v["stage2_cg_station_m"], v["stage2_length_m"], radius),
                fixed_cylinder("stage2_propellant", prop[1], -rate[1],
                               v["stage2_cg_station_m"], v["stage2_length_m"], radius,
                               mass_acceleration_kgps2=-rate_change[1]),
                fixed_cylinder("fairing", v["fairing_mass_kg"], 0.0,
                               v["fairing_cg_station_m"], v["fairing_length_m"], v["fairing_diameter_m"]/2),
                fixed_cylinder("payload", v["payload_mass_kg"], 0.0,
                               v["payload_cg_station_m"], v["payload_length_m"], v["payload_radius_m"]),
            ]
        return combine_components(components)
