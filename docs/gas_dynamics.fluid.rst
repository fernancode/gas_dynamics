###########################
Notes about the fluid class
###########################

The fluid class is used to hold all the information about a fluid. Importing air from the fluids module, we can see some of the built in properties.

.. code-block:: python

    >>> from gas_dynamics.fluids import air, air_us
    >>> for attr in air.__dict__:
    ...     print(attr, ':', air.__dict__[attr])
    ...
    name : Air
    gamma : 1.4
    cp : 1000
    cv : 716
    R : 286.9
    gc : 1
    units : J / kg-K
    rho : None
    temperature : None
    velocity : None
    a : None
    mach : None
    mass_veolcity : None
    >>>
    >>> for attr in air_us.__dict__:
    ...     print(attr, ':', air_us.__dict__[attr])
    ...
    name : Air
    gamma : 1.4
    cp : 0.24
    cv : 0.171
    R : 53.3
    gc : 32.174
    units : Btu / lbm-R
    rho : None
    temperature : None
    velocity : None
    a : None
    mach : None
    mass_veolcity : None
    >>>

When creating your own fluid, the properties that must be set when initiated are the fluid name, ratio of specific heats, gas constant, and a string of the units being used. Currently the string is meant to serve as a reminder to the user as to what units are being used, and no unit checking is done by the module in any way.

.. code-block:: python

    >>> foobar = fluid(name='foobar', gamma=1.4, R=53.3, units='Btu / lbm-R')

If we run out of the gates and immediately try and calculate the speed of sound of this fluid (which is remarkably similar to air) at standard temperature 491.67 Rankine, we notice the answer is considerably off from what we expect, ~~1100 ft/s.

.. code-block:: python

    >>> a = gd.sonic_velocity(gas=foobar, T=491.67)
    >>> a
    191.54220266040588
    >>>

The reason for this is we have yet to set the proportionality factor for our fluid which uses the US standard system. Currently it is set to 1 as a default, as for the metric system the conversion is unity.

.. code-block:: python

    >>> foobar.gc = 32.174
    >>> a = gd.sonic_velocity(gas=foobar, T=491.67)
    >>> a
    1086.4681666204492
    >>>

Currently supported fluids in metric and standard are

* Air
* Argon
* Carbon Dioxide
* Carbon Monoxide
* Helium
* Hydrogen
* Methane
* Nitrogen
* Oxygen
* Water Vapor

.. code-block:: python

    >>> from gas_dynamics import nitrogen, nitrogen_us


.. automodule:: gas_dynamics.fluid
   :members:
   :undoc-members:
   :show-inheritance: