.. gas_dynamics documentation master file, created by
   sphinx-quickstart on Sat Nov  7 14:41:31 2020.

Home
========================================

.. image:: https://badge.fury.io/py/gas-dynamics.svg
    :target: https://badge.fury.io/py/gas-dynamics

Package containing functions for working with compressible flow.

.. code-block:: python

   pip install gas-dynamics
   

.. figure:: README_images/shockwave.png
   :width: 800
   :alt: Shock waves forming across a T-38 Talon

   Credits: NASA Images

Last update: |today|

.. toctree::
   :caption: Introduction
   :hidden:

   gettingstarted

.. toctree::
   :caption: Review
   :hidden:

   review/gas_dynamics.tfsreview
   review/gas_dynamics.compressible_flow

.. toctree::
   :caption: Gas Dynamics
   :hidden:

   standard/gas_dynamics.standard   
   shocks/gas_dynamics.shocks
   prandtl_meyer/gas_dynamics.prandtl_meyer
   fanno/gas_dynamics.fanno
   rayleigh/gas_dynamics.rayleigh
   gas_dynamics.fluid
   gas_dynamics.extra

.. toctree::
   :caption: Misc
   :hidden:
   :glob:

   Misc/*