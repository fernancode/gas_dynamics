############################
Equation Map - Prandtl-Meyer
############################

:py:func:`prandtl_meyer_turn <gas_dynamics.prandtl_meyer.prandtl_meyer.prandtl_meyer_turn>`

.. math::

   \nu = \left( \frac{\gamma + 1}{\gamma -1} \right)^{\frac{1}{2}} \tan^{-1} \left[ \frac{\gamma-1}{\gamma+1} (M^2 -1) \right] ^{\frac{1}{2}} - \tan^{-1}(M^2 - 1)^{\frac{1}{2}}


:py:func:`prandtl_meyer_mach <gas_dynamics.prandtl_meyer.prandtl_meyer.prandtl_meyer_mach>` employs an equation solver to return the Mach number in the equation above.


:py:func:`mach_wave_angle <gas_dynamics.prandtl_meyer.prandtl_meyer.mach_wave_angle>`

.. math::

   \mu = \tan^{-1} \left( \frac{1} {(M^2 -1)^{1/2}} \right)
