#!usr/bin/env

#import isentropic functions and stagnation relations
from gas_dynamics.isentropic.isentropic import ( 
  sonic_velocity,
  stagnation_pressure,
  stagnation_temperature,
  stagnation_pressure_ratio,
  stagnation_temperature_ratio,
  stagnation_density_ratio,
  stagnation_ratio,
  stagnation_ratio_table,
  mach_from_pressure_ratio,
  mach_from_temperature_ratio,
  pressure_from_mach_ratio,
  temperature_from_mach_ratio,
  entropy_produced,
  mach_area_ratio_choked,
  mach_area_ratio,
  mach_from_area_ratio,
  mass_flux_max,
  mass_flux,
  plot_stagnation_ratios )

#import shock functions
from gas_dynamics.shocks.shocks import (
  shock_mach,
  shock_mach_before,
  shock_pressure_ratio,
  shock_mach_from_pressure_ratio,
  shock_temperature_ratio,
  shock_dv_a, 
  shock_stagnation_ratio, 
  shock_flow_deflection, 
  shock_angle, 
  shock_mach_given_angles, 
  prandtl_meyer_turn, 
  prandtl_meyer_mach, 
  shock_oblique_charts,
  shock_tables, 
  dirac_from_machs)

from gas_dynamics.fluids import fluid

#from .fanno import func1,
#from .rayleigh func 1