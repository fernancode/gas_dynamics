################################
Thermo and Fluid Sciences Review
################################

A brief review of thermodynamics.

The Laws
########

Before discussing the 0th, 1st, and 2nd laws, it is important to acknowledge the existence of a relation betweeen the properties of state. This is sometimes referred to as the 00 law. Behold the equation of state

.. math::

   p = \rho R T

0
*
The zeroth law states that two systems in thermal equilibrium with a third are in thermal equilbrium with each other.

1
*
The first law of thermodynamics deals with the conservation of energy and can be stated in a variety of ways. 
The simplest statement for a closed system is by far

.. math::
   \Sigma Q = \Sigma W

where :math:`\Sigma Q` is the heat transferred into the system, and :math:`\Sigma W` is the work transferred from the system.

On a unit mass basis, one can also say

.. math::
   \delta q = \delta w + de

We take special care to note :math:`\delta q` and :math:`\delta w` as inexact differentials, or path dependent, while :math:`de` is an exact differential and represents the change only between an initial and final state.

We can continue to expand the expression for the first law knowing that :math:`e \equiv u + \frac{v^2}{2} + gz` and say

.. math::
   \delta q = \delta w + u + \frac{v^2}{2} + gz

Finally we arrive at

.. math::
   0 = \delta q + \delta w + \Delta u + \frac{{\Delta v}^2}{2} + g \Delta z

All the equations above are different ways to say the same fundamental thing - energy is conserved. Heat and work are two types of energy in transit, heat being the transfer of energy due to a temperature gradient (and always from the higher temperature system to the lower), and work being the effect of elevating a mass in a gravity field.

2
*
The second law of thermodynamics points to the very direction of all processes in our universe. Carnot, Clausius, and Lord Kelvin all came to this significant conclusion while working to understand things like steam engines and machines. Kelvin and Planck stated that it is impossible for an engine operating in a cycle to produce net work when exchanging heat with only one temperature source. In laymans terms, this can be thought as a machine using energy from a thermal reservoir cannot operate unless there is another, lower, thermal reservoir.

Clausius determined that heat can never pass from a cold body to a warmer body without the aid of an external process, giving recognition of the irreversibility of certain processes such as fluid friction and heat transfer. All real processes are irreversible. For a fictitious reversible process, we can quantify this degredation of energy quality by defining entropy.

.. math::
   \Delta S \equiv \int_{}^{}\frac{\delta Q_r}{T}

.. math::
   ds \equiv \frac{\delta q_r}{T}

This defines a change in entropy in terms of a reversible addition of heat :math:`\delta q` (which we know is an impossible process) for a system temperature. We can write the definition of entropy as such because entropy is a state, and is independent of the process. Regardless, it may be more prudent to write the definition of entropy as 

.. math::
   ds \geq \frac{\delta q_r}{T}

or

.. math::
   ds = \frac{\delta q}{T} + ds

Where the first expression shows we understand that we are always trending towards greater entropy, and the second stating that the change in entropy is the heat we add in addition to irreversible contributions from friction, viscosity, thermal conductivity, and so on.


Property Relations, Entropy, and Perfect Gases
##############################################



Example Derivations
###################

Pumps
*****