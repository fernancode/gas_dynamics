############
Equation Map
############

:py:func:`rayleigh_pressure_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_pressure_ratio>`

.. math::

   \frac{p_{2}}{p_{1}} = \frac{ 1 + \gamma M_{1}^2}{ 1 + \gamma M_{2}^2}


:py:func:`rayleigh_temperature_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_temperature_ratio>`

.. math::

   \frac{T_{2}}{T_{1}} = \left( \frac{ 1 + \gamma M_{1}^2}{ 1 + \gamma M_{2}^2} \right) ^2 \frac{M_{2}^2}{M_{1}^2}


:py:func:`rayleigh_density_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_density_ratio>`

.. math::

   \frac{\rho_{2}}{\rho_{1}} = \frac{M_{1}^2}{M_{2}^2} \left( \frac{ 1 + \gamma M_{2}^2}{ 1 + \gamma M_{1}^2} \right)


:py:func:`rayleigh_stagnation_temperature_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_stagnation_temperature_ratio>`

.. math::

   \frac{T_{t2}}{T_{t1}} = \left( \frac{ 1 + \gamma M_{1}^2}{ 1 + \gamma M_{2}^2} \right) ^2 \frac{M_{2}^2}{M_{1}^2} \left(\frac{1 + \left[ (\gamma-1)/2 \right] M_{2}^2 }{ 1 + \left[ (\gamma-1)/2 \right] M_{1}^2 } \right)


:py:func:`rayleigh_stagnation_pressure_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_stagnation_pressure_ratio>`

.. math::

   \frac{p_{t2}}{p_{t1}} = \frac{ 1 + \gamma M_{1}^2}{ 1 + \gamma M_{2}^2} \left(\frac{1 + \left[ (\gamma-1)/2 \right] M_{2}^2 }{ 1 + \left[ (\gamma-1)/2 \right] M_{1}^2 } \right) ^ {\gamma/(\gamma-1)}


:py:func:`rayleigh_mach_from_pressure_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_mach_from_pressure_ratio>`

.. math::

   M_{2} = \left[ \left( \frac{p_{2}}{p_{1}} \left( 1+\gamma M_{1}^2 \right) - 1 \right) \frac{1}{\gamma} \right]^{1/2}


:py:func:`rayleigh_mach_from_temperature_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_mach_from_temperature_ratio>` , :py:func:`rayleigh_mach_from_stagnation_temperature_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_mach_from_stagnation_temperature_ratio>` , :py:func:`rayleigh_mach_from_stagnation_pressure_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_mach_from_stagnation_pressure_ratio>` all use equation solvers to get the desired result.


:py:func:`rayleigh_pressure_star_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_pressure_star_ratio>`

.. math::

   \frac{p}{p^*} = \frac{1+\gamma}{1+\gamma M^2}


:py:func:`rayleigh_temperature_star_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_temperature_star_ratio>`

.. math::

   \frac{T}{T^*} = \frac{M^2 (1+\gamma)^2}{ (1+\gamma M^2)^2 }


:py:func:`rayleigh_density_star_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_density_star_ratio>`

.. math::

   \frac{\rho}{\rho^*} = \frac{1+\gamma M^2}{ (1+\gamma) M^2}


:py:func:`rayleigh_stagnation_temperature_star_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_stagnation_temperature_star_ratio>`

.. math::

   \frac{T_{t}}{T_{t}^*} = \frac{2 (1+\gamma)^2 M^2}{ (1+\gamma M^2)^2 } \left( 1+ \frac{\gamma-1}{2}M^2\right)


:py:func:`rayleigh_stagnation_pressure_star_ratio <gas_dynamics.rayleigh.rayleigh.rayleigh_stagnation_pressure_star_ratio>`

.. math::

   \frac{p_{t}}{p_{t}^*} = \frac{1+\gamma}{1+\gamma M^2} \left( \frac{1+ \left[ (\gamma -1)/2 \right] M^2 }{ (\gamma +1 )/ 2} \right) ^ {\frac{\gamma}{\gamma -1}}


:py:func:`rayleigh_heat_flux <gas_dynamics.rayleigh.rayleigh.rayleigh_heat_flux>`

.. math::

   q = c_{p} (T_{t2} - T_{t1})