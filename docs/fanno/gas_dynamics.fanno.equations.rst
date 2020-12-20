####################
Equation Map - Fanno
####################

:py:func:`stagnation_enthalpy <gas_dynamics.fanno.fanno.stagnation_enthalpy>`

.. math::

   h_{t} = h + \frac{G^2}{\rho^2 2}


:py:func:`fanno_temperature_ratio <gas_dynamics.fanno.fanno.fanno_temperature_ratio>`

.. math::

   \frac{T_{2}}{T_{1}} = \frac{1 + \left[ (\gamma-1)/2 \right] M_{1}^2 }{ 1 + \left[ (\gamma-1)/2 \right] M_{2}^2 }


:py:func:`fanno_pressure_ratio <gas_dynamics.fanno.fanno.fanno_pressure_ratio>`

.. math::

   \frac{p_{2}}{p_{1}} = \frac{M_{1}}{M_{2}} \left( \frac{1 + \left[ (\gamma-1)/2 \right] M_{1}^2 }{ 1 + \left[ (\gamma-1)/2 \right] M_{2}^2 } \right) ^{1/2}


:py:func:`fanno_temperature_ratio <gas_dynamics.fanno.fanno.fanno_temperature_ratio>`

.. math::

   \frac{\rho_{2}}{\rho_{1}} = \frac{M_{1}}{M_{2}} \left( \frac{1 + \left[ (\gamma-1)/2 \right] M_{2}^2 }{ 1 + \left[ (\gamma-1)/2 \right] M_{1}^2 } \right) ^{1/2}


:py:func:`fanno_stagnation_pressure_ratio <gas_dynamics.fanno.fanno.fanno_stagnation_pressure_ratio>`

.. math::

   \frac{p_{t2}}{p_{t1}} = \frac{M_{1}}{M_{2}} \left( \frac{1 + \left[ (\gamma-1)/2 \right] M_{2}^2 }{ 1 + \left[ (\gamma-1)/2 \right] M_{1}^2 } \right) ^{\frac{\gamma+1}{2(\gamma-1)}}


:py:func:`fanno_temperature_star_ratio <gas_dynamics.fanno.fanno.fanno_temperature_star_ratio>`

.. math::

   \frac{T}{T^*} = \frac{(\gamma+1)/2 }{ 1 + \left[ (\gamma-1)/2 \right] M^2}


:py:func:`fanno_pressure_star_ratio <gas_dynamics.fanno.fanno.fanno_pressure_star_ratio>`

.. math::

   \frac{p}{p^*} = \frac{1}{M} \left( \frac{(\gamma+1)/2 }{ 1 + \left[ (\gamma-1)/2 \right] M^2} \right) ^{1/2}


:py:func:`fanno_density_star_ratio <gas_dynamics.fanno.fanno.fanno_density_star_ratio>`

.. math::

   \frac{\rho}{\rho^*} = \frac{1}{M} \left( \frac{ 1 + \left[ (\gamma-1)/2 \right] M^2} {(\gamma+1)/2 } \right) ^ {1/2}


:py:func:`fanno_velocity_star_ratio <gas_dynamics.fanno.fanno.fanno_velocity_star_ratio>`

.. math::

   \frac{V}{V^*} = \frac{M}{1} \left( \frac {(\gamma+1)/2 }{ 1 + \left[ (\gamma-1)/2 \right] M^2} \right) ^ {1/2}


:py:func:`fanno_parameter <gas_dynamics.fanno.fanno.fanno_parameter>`

.. math::

   \frac{ f(x_{2} - x_{1}) }{D_{e}} = \frac{\gamma + 1}{2\gamma} \ln \frac{ 1 + \left[ (\gamma-1)/2 \right] M_{2}^2 }{1 + \left[ (\gamma-1)/2 \right] M_{1}^2 } - \frac{1}{\gamma} \left( \frac{1}{M_{2}^2} - \frac{1}{M_{1}^2} \right) - \frac{\gamma + 1}{2\gamma} \ln \frac{M_{2}^2}{M_{1}^2}


:py:func:`fanno_parameter_max <gas_dynamics.fanno.fanno.fanno_parameter_max>`

.. math::

   \frac{f(x - x^*)} {D_{e}} = \frac{\gamma + 1}{2\gamma} \ln \frac{ 1 + \left[ (\gamma-1)/2 \right] M^2} {(\gamma+1)/2 } - \frac{1}{\gamma} \left( \frac{1}{M^2} - 1 \right) - \frac{\gamma + 1}{2\gamma} \ln M^2


:py:func:`mach_from_fanno <gas_dynamics.fanno.fanno.mach_from_fanno>` uses an equation solver to determine the mach number.