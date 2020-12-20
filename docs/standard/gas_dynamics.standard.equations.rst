#######################
Equation Map - Standard
#######################

:py:func:`sonic_velocity <gas_dynamics.standard.standard.sonic_velocity>`

.. math::

   a = \sqrt{\gamma R T}


===============================
Stagnation Relations and Ratios
===============================

:py:func:`stagnation_pressure <gas_dynamics.standard.standard.stagnation_pressure>` will return :math:`P`, :math:`P_{t}`, or :math:`M` depending on parameters given .

.. math::

   p_{t} = p\left(1+\frac{\gamma-1}{2} M^{2}\right)^{\frac{\gamma}{\gamma-1}}


:py:func:`stagnation_temperature <gas_dynamics.standard.standard.stagnation_temperature>` will return :math:`T`, :math:`T_{t}`, or :math:`M` depending on parameters given.

.. math::

   T_{t} = T\left(1 + \frac{\gamma-1}{2} M^{2}\right)


:py:func:`stagnation_pressure_ratio <gas_dynamics.standard.standard.stagnation_pressure_ratio>` returns :math:`\frac{P}{P_{t}}`

.. math::
   \frac{p}{p_{t}} = \frac{1}{\left(1 + \frac{\gamma-1}{2}M^2 \right)^\frac{\gamma}{\gamma-1}}



:py:func:`stagnation_temperature_ratio <gas_dynamics.standard.standard.stagnation_temperature_ratio>` returns :math:`\frac{T}{T_{t}}`

.. math::

   \frac{T}{T_{t}} = \frac{1}{\left(1 + \frac{\gamma-1}{2} M^{2}\right)}


:py:func:`stagnation_density_ratio <gas_dynamics.standard.standard.stagnation_density_ratio>` returns :math:`\frac{\rho}{\rho_{t}}`.

.. math::

   \frac{\rho}{\rho_{t}} = \left( \frac{1}{1+\frac{\gamma-1}{2} M^{2}} \right)^{\frac{1}{\gamma-1}}


==================================
Mach Number and Property Relations
==================================

:py:func:`mach_from_pressure_ratio <gas_dynamics.standard.standard.mach_from_pressure_ratio>` solves for :math:`M_{2}` from the following equation, while :py:func:`pressure_from_mach_ratio <gas_dynamics.standard.standard.pressure_from_mach_ratio>` will solve for :math:`p_{2}`.

.. math::

   \frac{p_{2}}{p_{1}} = \left( \frac{ 1 + \frac{\gamma-1}{2}M_{1}^2}{1 + \frac{\gamma-1}{2}M_{2}^2} \right)^{\frac{\gamma}{\gamma-1}}e^{\frac{-\Delta s}{R}}


:py:func:`mach_from_temperature_ratio <gas_dynamics.standard.standard.mach_from_temperature_ratio>` solves for :math:`M_{2}` from the following equation, while :py:func:`temperature_from_mach_ratio <gas_dynamics.standard.standard.pressure_from_mach_ratio>` will solve for :math:`T_{2}`.

.. math::

   \frac{T_{2}}{T_{1}} = \frac{1 + \frac{\gamma-1}{2}M_{1}^2}{1 + \frac{\gamma-1}{2}M_{2}^2}


:py:func:`mach_area_star_ratio <gas_dynamics.standard.standard.mach_area_star_ratio>` returns the ratio of :math:`\frac{A}{A*}`.

.. math::

   \frac{A}{A*} = \frac{1}{M} \left( \frac{1 + \frac{\gamma-1}{2} M^2}{ \frac{\gamma+1}{2}} \right)^{\frac{\gamma+1}{2(\gamma-1)}}


:py:func:`mach_area_ratio <gas_dynamics.standard.standard.mach_area_ratio>` returns the ratio of :math:`\frac{A_{2}}{A_{1}}` given two Mach numbers, whi;e :py:func:`mach_from_area_ratio <gas_dynamics.standard.standard.mach_from_area_ratio>` will return the possible mach numbers that satisfy the area ratio.

.. math::

   \frac{A_{2}}{A_{1}} = \frac{M_{1}}{M_{2}} \left( \frac{1+\frac{\gamma-1}{2}M_{2}^2}{1+\frac{\gamma-1}{2}M_{1}^2}\right)^{\frac{\gamma+1}{2(\gamma-1)}}


=========
Mass Flux
=========

:py:func:`mass_flux <gas_dynamics.standard.standard.mass_flux>` returns the flow rate per unit area while :py:func:`mass_flux_max <gas_dynamics.standard.mass_flux_max>` will return the maximum flow rate per unit area, where :math:`M=1`.

.. math::

   \frac{\dot{m}}{A}=M\left(1+\frac{\gamma-1}{2}M^2\right)^{\frac{-(\gamma+1)}{2(\gamma-1)}}\sqrt{\left(\frac{\gamma}{R}\right)}\frac{p_{t}}{\sqrt{T_{t}}}


.. math::

   \frac{\dot{m}}{A^*} = \sqrt{\frac{\gamma}{R}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}}\frac{p_{t}}{\sqrt{T_{t}}} 
