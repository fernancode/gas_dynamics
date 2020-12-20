#####################
Equation Map - Shocks
#####################

:py:func:`shock_mach <gas_dynamics.shocks.shocks.shock_mach>`

.. math::

   M_{2} = \left[ \frac{ M_{1}^2 + 2 / (\gamma -1) }{ \left[ 2 \gamma /( \gamma -1 ) \right] M_{1}^2 -1 } \right]^ {\frac{1}{2}}


:py:func:`shock_mach_before <gas_dynamics.shocks.shocks.shock_mach_before>`

.. math::

   M_{2} = \left[ \frac{ M_{2}^2 + 2 / (\gamma -1) }{ \left[ 2 \gamma /( \gamma -1 ) \right] M_{2}^2 - 1 } \right]^ {\frac{1}{2}}

   
:py:func:`shock_pressure_ratio <gas_dynamics.shocks.shocks.shock_pressure_ratio>`

.. math::

   \frac{p_{2}}{p_{1}} = \frac{ 2 \gamma }{ \gamma + 1} M_{1}^2 - \frac{ \gamma - 1 }{ \gamma + 1}


:py:func:`shock_mach_from_pressure_ratio <gas_dynamics.shocks.shocks.shock_mach_from_pressure_ratio>`

.. math::

   M = \left[\frac{\gamma+1}{2\gamma} \left(\frac{p_{2}}{p_{1}}\right)+\frac{\gamma-1}{\gamma+1}\right]^{\frac{1}{2}}


:py:func:`shock_temperature_ratio <gas_dynamics.shocks.shocks.shock_temperature_ratio>`

.. math::

   \frac{T_{2}}{T_{1}} = \frac{\left( 1 + \left[ \left( \gamma - 1 \right) /2 \right] M_{1}^2 \right) \left( \left[ 2 \gamma / \left( \gamma -1 \right) \right] M_{1}^2 -1 \right)}{ \left[ \left( \gamma + 1 \right)^2 / 2 \left(\gamma - 1 \right) \right] M_{1}^2 }


:py:func:`shock_dv_a <gas_dynamics.shocks.shocks.shock_dv_a>`

.. math::

   \frac{dV}{a_{1}} = \left( \frac{2}{\gamma + 1} \right) \left( \frac{ M_{1}^2 -1} {M_{1}} \right)


:py:func:`shock_stagnation_pressure_ratio <gas_dynamics.shocks.shocks.shock_stagnation_pressure_ratio>`

.. math::

   \frac{p_{t2}}{p_{t1}} = \left( \frac{\left[ \left( \gamma + 1 \right) /2 \right] M_{1}^2} { 1 + \left[ \left( \gamma - 1 \right) /2 \right] M_{1}^2 } \right)^ { \frac{\gamma}{\gamma -1}} \left[ \frac{2\gamma}{\gamma+1} M_{1}^2 - \frac{\gamma -1}{ \gamma +1}\right] ^ {\frac{1}{1-\gamma}}


:py:func:`shock_flow_deflection <gas_dynamics.shocks.shocks.shock_flow_deflection>`

.. math::

   \delta = \arctan \left[ 2 \cot(\theta) \left( \frac{ M_{1}^2 \sin^2 (\theta) - 1}{ M_{1}^2 (\gamma + \cos 2\theta) + 2 } \right) \right]


:py:func:`shock_angle <gas_dynamics.shocks.shocks.shock_angle>` , :py:func:`shock_mach_given_angles <gas_dynamics.shocks.shocks.shock_mach_given_angles>` , and :py:func:`shock_flow_deflection_from_machs <gas_dynamics.shocks.shocks.shock_flow_deflection_from_machs>` all employ equation solvers with combinations of the above functions to return the desired values.