#!usr/bin/env
#Equations, plots, and tables for working with constant area flow
#subject to losses from friction, otherwise known as Fanno flow
#
#
#Copyright 2020 by Fernando A de la Fuente
#All rights reserved

from gas_dynamics.fluids import fluid, air
from numpy import log
from scipy.optimize import fsolve



#==================================================
#stagnation enthalpy
#==================================================
def stagnation_enthalpy(enthalpy: float, gas=air) -> float:
    """Return the stagnation enthalpy

    Notes
    -----
    Given the fluid state and a given enthalpy, return its stagnation
    enthalpy

    Parameters
    ----------
    enthalpy : `float`
        The enthalpy of the fluid
    gas : `fluid`
        A user defined fluid object

    Examples
    --------

    """

    ht = enthalpy + gas.mass_velocity**2 / (gas.rho**2 * 2 * gas.gc)
    return ht



#==================================================
#fanno temperature
#==================================================


#==================================================
#fanno mach from temperature
#==================================================



#==================================================
#fanno pressure
#==================================================


#==================================================
#fanno mach from pressure
#==================================================



#==================================================
#fanno temperature ratio
#==================================================
def fanno_temperature_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the temperature ratio for a fanno flow given two Mach numbers

    Notes
    -----
    Given the two Mach numbers of a constant area adiabatic duct under the influence of
    friction alone, return the temperature ratio of region two over region one. 
    Default fluid is air.

    Parameters
    ----------
    M1 : `flaot`
        The mach number at region 1 \n
    M2 : `float`
        The mach number at region 2 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M1, M2 = 1.2, 1
    >>> T2_T1 = gd.fanno_temperature_ratio(M2,M1)
    >>> T2_T1
    0.9316770186335405
    >>>
    """

    gamma = gas.gamma
    T2_T1 = ( 1+ (gamma-1)/2 * M1**2)/( 1+ (gamma-1)/2 * M2**2)
    return T2_T1



#==================================================
#fanno pressure ratio
#==================================================
def fanno_pressure_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the pressure ratio for a fanno flow given two Mach numbers

    Notes
    -----
    Given the two Mach numbers of a constant area adiabatic duct under the influence of
    friction alone, return the pressure ratio of region two over region one. 
    Default fluid is air.

    Parameters
    ----------
    M1 : `flaot`
        The mach number at region 1 \n
    M2 : `float`
        The mach number at region 2 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M1, M2 = 1.2, 1
    >>> p2_p1 = gd.fanno_pressure_ratio(M1,M2)
    >>> p2_p1
    1.243221621433604
    >>>
    """

    gamma = gas.gamma
    p2_p1 = M1/M2 * (( 1+ (gamma-1)/2 * M1**2)/( 1+ (gamma-1)/2 * M2**2))**.5
    return p2_p1



#==================================================
#fanno density ratio
#==================================================
def fanno_density_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the density ratio for a fanno flow given two Mach numbers

    Notes
    -----
    Given the two Mach numbers of a constant area adiabatic duct under the influence of
    friction alone, return the density ratio of region two over region one. 
    Default fluid is air.

    Parameters
    ----------
    M1 : `flaot`
        The mach number at region 1 \n
    M2 : `float`
        The mach number at region 2 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M1, M2 = 1.2, 1
    >>> rho2_rho1 = gd.fanno_density_ratio(M2,M1)
    >>> rho2_rho1
    0.8633483482177806
    >>>
    """

    gamma = gas.gamma
    rho2_rho1 = M1/M2 * ((1+(gamma-1)/2*M2**2)/(1+(gamma-1)/2*M1**2))**.5
    return rho2_rho1



#==================================================
#fanno stagnation star ratio
#==================================================
def fanno_stagnation_pressure_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the stagnation pressure ratio pt2/pt1 for a fanno flow given two mach numbers

    Notes
    -----
    Given the two Mach numbers of a constant area adiabatic duct under the influence of
    friction alone, return the stagnation pressure ratio of region two over region one. 
    Default fluid is air.

    Parameters
    ----------
    M1 : `flaot`
        The mach number at region 1 \n
    M2 : `float`
        The mach number at region 2 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M1, M2 = 1.2, 1
    >>> pt2_pt1 = gd.fanno_stagnation_pressure_ratio(M2,M1)
    >>> pt2_pt1
    1.0304397530864196
    >>>
    """
    
    gamma = gas.gamma
    pt2_pt1 = M1/M2 * (( 1+ (gamma-1)/2 * M2**2)/( 1+ (gamma-1)/2 * M1**2))**((gamma+1)/(2*(gamma-1)))
    return pt2_pt1



#==================================================
#fanno temperature star ratio
#==================================================
def fanno_temperature_star_ratio(M: float, gas=air) -> float:
    """Return the ratio of temperature over temperature where Mach equals one

    Notes
    -----
    Given a Mach number of a constant area adiabatic duct under the influence of
    friction alone, return the temperature ratio of region two over region one where 
    Mach in region one equals one. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The mach number
    gas : `fluid`
        The user defined fluid object

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 1.2
    >>> gd.fanno_temperature_star_ratio(M)
    0.9316770186335405
    >>>
    """

    gamma = gas.gamma
    T_Tstar = ((gamma+1)/2)/(1 + (gamma-1)/2 * M**2)
    return T_Tstar



#==================================================
#fanno pressure star ratio
#==================================================
def fanno_pressure_star_ratio(M: float, gas=air) -> float:
    """Return the ratio of pressure over pressure where Mach equals one

    Notes
    -----
    Given a Mach number of a constant area adiabatic duct under the influence of
    friction alone, return the pressure ratio of region two over region one where 
    Mach in region one equals one. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The mach number
    gas : `fluid`
        The user defined fluid object

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 1.2
    >>> gd.fanno_pressure_star_ratio(M)
    0.8043618151097336
    >>>
    """

    gamma = gas.gamma
    p_pstar = 1/M * (((gamma+1)/2)/(1 + (gamma-1)/2 * M**2))**.5
    return p_pstar



#==================================================
#fanno density star ratio
#==================================================
def fanno_density_star_ratio(M: float, gas=air) -> float:
    """Return the ratio of density over density where Mach equals one

    Notes
    -----
    Given a Mach number of a constant area adiabatic duct under the influence of
    friction alone, return the density ratio of region two over region one where 
    Mach in region one equals one. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The mach number
    gas : `fluid`
        The user defined fluid object

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 1.2
    >>> gd.fanno_density_star_ratio(M)
    0.8633483482177806
    >>>
    """

    gamma = gas.gamma
    rho_rhostar = 1/M * ((1 + (gamma-1)/2 * M**2)/((gamma+1)/2))**.5
    return rho_rhostar



#==================================================
#fanno velocity choked ratio
#==================================================
def fanno_velocity_star_ratio(M: float, gas=air) -> float:
    """Return the ratio of velocity over velocity where Mach equals one

    Notes
    -----
    Given a Mach number of a constant area adiabatic duct under the influence of
    friction alone, return the velocity ratio of region two over region one where 
    Mach in region two equals one. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The mach number
    gas : `fluid`
        The user defined fluid object

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 1.2
    >>> v_vstar = gd.fanno_velocity_star_ratio(M)
    >>> v_vstar
    1.1582810137580164
    >>>
    """

    gamma = gas.gamma
    v_vstar = M/1 * (((gamma+1)/2)/(1 + (gamma-1)/2 * M**2))**.5
    return v_vstar



#==================================================
#fanno
#==================================================
def fanno_parameter(M1: float, M2: float, gas=air) -> float:
    """Return the product of friction factor and length divided by diameter for two Mach numbers

    Notes
    -----
    Given the two Mach numbers of a constant area adiabatic duct under the influence of
    friction alone, return the fanno parameter that describes that system where fanno 
    parameter is the product of friction factor and length over diameter.
    Default fluid is air. 

    Parameters
    ----------
    M1 : `flaot`
        The mach number at region 1
    M2 : `float`
        The mach number at region 2
    gas : `fluid`
        The user defined fluid object
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M1, M2 = 3, 2 
    >>> fanno = gd.fanno_parameter(M1, M2)
    >>> fanno
    0.21716290559704166
    >>>
    """

    gamma = gas.gamma
    fanno = (gamma+1)/(2*gamma) * log((( 1+ (gamma-1)/2 * M2**2)/( 1+ (gamma-1)/2 * M1**2))) - 1/gamma * (1/(M2**2) - 1/(M1**2)) - (gamma+1)/(2*gamma) * log((M2**2)/(M1**2))
    return fanno



#==================================================
#fanno parameter max
#==================================================
def fanno_parameter_max(M: float, gas=air) -> float:
    """Return the maximum product of friction factor and length divided by diameter

    Notes
    -----
    Given a Mach number of a constant area adiabatic duct under the influence of
    friction alone, determine the maximum length to diameter ratio for a fluid to 
    reach a Mach number of 1. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The starting Mach number
    gas : `fluid`
        The user defined fluid object

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 2
    >>> fanno_max = gd.fanno_parameter_max(M)
    >>> fanno_max
    0.3049965025814798
    >>>    
    """

    gamma = gas.gamma
    fanno_ratio_max = (gamma + 1)/(2*gamma) * log(((gamma+1)/2 * M**2) / (1 + (gamma-1)/2 * M**2)) + 1/gamma * (1/(M**2)-1)
    return fanno_ratio_max



#==================================================
#mach from fanno parameter
#==================================================
def mach_from_fanno(fanno: float, M1: float, gas=air) -> float:
    """Return the Mach number that would result from the fanno parameter and initial mach number

    Notes
    -----
    Given the Mach number and fanno parameter that describes that system, return the resulting
    mach number. Default fluid is air.


    Parameters
    ----------
    fanno : `float`
        The fanno parameter for the system
    M : `float`
        The starting Mach number
    gas : `fluid`
        The user defined fluid object

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> fanno, M1 = .3, 2.64
    >>> M2 = gd.mach_from_fanno(fanno=fanno, M1=M1)
    >>> M2
    1.567008305615555
    >>>
    """

    def mach_solve(M, M1=M1, fanno=fanno, gas=gas):
        zero = fanno_parameter(M1=M1, M2=M, gas=gas) - fanno
        return zero

    if M1 < 1:
        x0 = .5
    elif M1 > 1:
        x0 = 1.5
    else:
        x0=1

    sol = fsolve(mach_solve, args=(M1, fanno, gas), x0=x0)
    return sol[0]