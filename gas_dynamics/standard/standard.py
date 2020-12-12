#!usr/bin/env
#Equations, plots, and tables for solving compresible flow problems.
#Functions range from stagnation relations to determining the mach
#number from changes in local properties to determining the maximum mass
#flux. Tables can be made for any gas and its respective ratio of specific
#heats, as well as plots and charts of relationships.

#  Typical usage example:
#  Generate isentropic flow stagnation relations for methane
#  >>> import gas_dynamics as gd
#  >>> gd.stagnation_ratio_table(range=[0,2], step=.2, gas=methane)
#  M: 0.000   |   P/Pt: 1.000    |    T/Tt: 1.000    |    A/A*: inf    |   rho/rho_t: 1.000
#  M: 0.200   |   P/Pt: 0.974    |    T/Tt: 0.994    |    A/A*: 2.988    |   rho/rho_t: 0.980
#  M: 0.400   |   P/Pt: 0.901    |    T/Tt: 0.975    |    A/A*: 1.600    |   rho/rho_t: 0.924
#  M: 0.600   |   P/Pt: 0.794    |    T/Tt: 0.946    |    A/A*: 1.192    |   rho/rho_t: 0.839
#  M: 0.800   |   P/Pt: 0.669    |    T/Tt: 0.907    |    A/A*: 1.039    |   rho/rho_t: 0.737
#  M: 1.000   |   P/Pt: 0.542    |    T/Tt: 0.862    |    A/A*: 1.000    |   rho/rho_t: 0.629
#  M: 1.200   |   P/Pt: 0.425    |    T/Tt: 0.813    |    A/A*: 1.032    |   rho/rho_t: 0.523
#  M: 1.400   |   P/Pt: 0.325    |    T/Tt: 0.761    |    A/A*: 1.121    |   rho/rho_t: 0.426
#  M: 1.600   |   P/Pt: 0.243    |    T/Tt: 0.709    |    A/A*: 1.267    |   rho/rho_t: 0.342
#  M: 1.800   |   P/Pt: 0.179    |    T/Tt: 0.659    |    A/A*: 1.474    |   rho/rho_t: 0.271
#  M: 2.000   |   P/Pt: 0.130    |    T/Tt: 0.610    |    A/A*: 1.754    |   rho/rho_t: 0.213
#
#Copyright 2020 by Fernando A de la Fuente
#All rights reserved

#TODO: everything needs a 
#returns
#   `type`


import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from gas_dynamics.fluids import fluid, air, methane, argon



#==================================================
#sonic_velocity
#fluid class and us standard option 11/6/2020
#==================================================    
def sonic_velocity(temperature=273.15, gas=air) -> float:
    """Returns the local speed of sound.
    
    Notes
    -----
    Given a ratio of specific heats, gas constant, and temperature
    this function returns the locoal speed of sound. Default fluid is air.
    
    Parameters
    ----------
    gas : `fluid`
        A user defined fluid object. Default is air \n
    temperature : `float`
        The temperature \n

    Returns
    -------
    float
        The local speed of sound\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.sonic_velocity(air, temperature=500)
    448.2186966202994
    >>> gd.sonic_velocity(gas=Argon, temperature=300)

    """

    gamma, R, gc = gas.gamma, gas.R, gas.gc
    a = ( gc*gamma*R*temperature)**.5
    return a



#==================================================
#stagnation_pressure
#implemented fluid class and string for output returned
#==================================================    
def stagnation_pressure(stagnation_pressure=None, mach=None, pressure=None, gas=air, output=False) -> float:
    """Returns the stagnation pressure given pressure and Mach number.

    Notes
    -----
    Given a pressure, Mach number, and a ratio of specific heats return
    the stagnation pressure. Alternatively, provided two arguments
    the function will return the missing one. Default fluid is air.

    Parameters
    ----------
    stagnation_pressure : `float`
        The stagnation pressure.\n
    pressure : `float`
        The pressure.\n
    mach : `float`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    output : `bool`
        Print out a string to verify the output is the parameter desired \n
    
    Returns
    -------
    float
        The stagnation pressure, static pressure, or mach number\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> pt = gd.stagnation_pressure(pressure=10, mach=1)
    >>> pt
    18.92929158737854
    >>> M = gd.stagnation_pressure(pressure=10, stagnation_pressure=pt)
    >>> M
    1.0
    >>>
    """

    pt = stagnation_pressure
    p = pressure
    gamma = gas.gamma
    if pt == None:
        pt = p* ( 1 + (gamma-1)/2 * mach**2)** (gamma/(gamma-1))
        if output == True:
            print('Returned stagnation pressure')
        return pt
   
    if mach == None:
        mach = (((pt/p)**((gamma-1)/gamma) -1 ) * 2/(gamma-1) ) ** .5
        if output == True:
            print('Returned Mach')
        return mach

    if p == None:
        p = pt / ( 1 + (gamma-1)/2 * mach**2)** (gamma/(gamma-1))
        if output == True:
            print('Returned pressure')
        return p



#==================================================
#stagnation_temperature
#implemented output string, fluid class, 
#==================================================    
def stagnation_temperature(temperature=None, stagnation_temperature=None , mach=None, gas=air, output=False) -> float :
    """Returns the stagnation temperature given temperature and Mach number.
    
    Notes
    -----
    Given a temperature, Mach number, and a ratio of specific heats 
    this function returns the stagnation temperature. Alternatively,
    provided two arguments the function will return the missing one.
    Default fluid is air.

    Parameters
    ----------
    stagnation_temperature : `float`
        The stagnation temperature\n
    temperature : `float`
        The temperature\n
    mach : `float`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    output : `bool`
        Print out a string to verify the output is the parameter desired \n
    
    Returns
    -------
    float
        The stagnation temperature, temperature, or mach number\n

    Examples
    --------
    >>> Tt = gd.stagnation_temperature(temperature=300, mach=1)
    >>> Tt
    360.0
    >>> M = gd.stagnation_temperature(temperature=300, stagnation_temperature=Tt)
    >>> M 
    1.0
    >>>
    """

    Tt = stagnation_temperature
    T = temperature
    gamma = gas.gamma
    if Tt == None:
        Tt = T * ( 1 + (gamma-1)/2 * mach**2)
        if output == True:
            print('Returned stagnation temperature')
        return Tt

    if mach == None:
        mach = ((Tt /T - 1) * 2/(gamma-1))**.5
        if output == True:
            print('Returned Mach')
        return mach
        
    if T == None:
        T = Tt/( 1 + (gamma-1)/2 * mach**2)
        if output == True:
            print('Returned temperature')
        return T



#==================================================
#stagnation_density
#implemented output string, fluid class, 
#TODO: FIX ME
#==================================================    
def stagnation_density(density=None, stagnation_density=None , mach=None, gas=air, output=False) -> float :
    """Returns the stagnation density given density and Mach number.
    
    Notes
    -----
    Given a density, Mach number, and a ratio of specific heats 
    this function returns the stagnation density. Alternatively,
    provided two arguments the function will return the missing one.
    Default fluid is air.

    Parameters
    ----------
    stagnation_density : `float`
        The stagnation density\n
    density : `float`
        The density\n
    mach : `float`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    output : `bool`
        Print out a string to verify the output is the parameter desired \n
    
    Returns
    -------
    float
        The stagnation temperature, temperature, or mach number\n

    Examples
    --------

    """

    rho_t = stagnation_density
    rho = density
    gamma = gas.gamma
    if rho_t == None:
        rho_t = density *  ( 1 + (gamma-1)/2 * mach**2) ** (1/(gamma-1))
        if output == True:
            print('Returned stagnation density')
        return Tt

    if mach == None:
        
        if output == True:
            print('Returned Mach')
        return mach
        
    if rho == None:
        
        if output == True:
            print('Returned temperature')
        return T



#==================================================
#stagnation_pressure_ratio
#added fluid class, removed unnnecessary metric argument
#==================================================
def stagnation_pressure_ratio(mach: float, gas=air) -> float:
    """Returns the pressure ratio of p / p_t
    
    Notes
    -----
    Given a Mach number and ratio of specific heats return the relation of
    pressure over stagnation pressure. Default fluid is air.

    Parameters
    ----------
    mach : `float`
        The Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The stagnation pressure ratio\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> p_pt = gd.stagnation_pressure_ratio(mach=3)
    >>> p_pt
    0.027223683703862817
    >>>
    """

    gamma = gas.gamma
    denom = 1 + (gamma-1)/2 * mach**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio



#==================================================
#stagnation_temperature_ratio
#added fluid class
#==================================================
def stagnation_temperature_ratio(mach: float, gas=air) -> float:
    """Returns the temperature ratio of T / T_t
    
    Notes
    -----
    Given a Mach number and ratio of specific heats return the relation of
    temperature over stagnation temperature.  Default fluid is air.

    Parameters
    ----------
    mach : `float`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air\n
    
    Returns
    -------
    float
        The stagnation temperature ratio\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> T_Tt = gd.stagnation_temperature_ratio(mach=1.5)
    >>> T_Tt
    0.6896551724137931
    >>>
    """

    gamma = gas.gamma
    Tt_ratio = 1 / (1+(gamma-1)/2 *mach**2)
    return Tt_ratio



#==================================================
#stagnation_density_ratio
#added fluid class
#==================================================
def stagnation_density_ratio(mach: float, gas=air) -> float:
    """Returns the density ratio rho / rho_t
    
    Notes
    -----
    Given a Mach number and ratio of specific heats, return the relation
    of density over stagnation density. Default fluid is air.
    
    Parameters
    ----------
    mach : `float`
        The Mach # \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The stagnation density ratio\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> rho_rho_t = gd.stagnation_density_ratio(mach=1.5) 
    >>> rho_rho_t
    0.39498444639115327
    >>>
    """

    gamma = gas.gamma
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * mach**2 )) ** (1 / (gamma-1))
    return rho_t_ratio



#==================================================
#stagnation ratio 
#==================================================
def stagnation_ratio(mach: float, gas=air) -> list:
    """Return stagnation pressure, temperature, density, and choked area ratio for a mach number
        
    Notes
    -----
    Given a mach number and the fluid, return the three stagnation ratios
    and the ratio of the area to the choked area. Default fluid is air.

    Parameters
    ----------
    mach : `float`
        The mach number
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    list
        The stagnation pressure, stagnation temperature, stagnation
        density ratio, and choked area ratio.\n

    Examples
    --------
    >>> import gas_dynamics as gd

    """

    p_pt = stagnation_pressure_ratio(mach=mach, gas=gas)
    T_Tt = stagnation_temperature_ratio(mach=mach, gas=gas)
    rho_rhot = stagnation_density_ratio(mach=mach, gas=gas)
    A_Astar = mach_area_star_ratio(mach=mach, gas=gas)
    return p_pt, T_Tt, rho_rhot, A_Astar



#==================================================
#stagnation ratio tables
#added fluid object
#==================================================
def stagnation_ratio_table(range=[0,5], step=.1, gas=air) -> str:
    """Returns the isentropic flow tables in the given range.
    
    Notes
    -----
    Given a ratio of specific heats, print out the stagnation
    temperature ratio, stagnation pressure ratio, the area to 
    choked area ratio, and the stagnation density ratio for every
    incremental Mach number.
    
    Parameters
    ----------
    range : `list`
        The starting and ending Mach numbers in a list, ex: [0, 5] \n
    step : `float`
        The step size between min and max mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    str
        The isentropic flow table\n
        
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.stagnation_ratios(range=[0,2], step=.2, gas='nitrogen')  
    M: 0.000   |   P/Pt: 1.000    |    T/Tt: 1.000    |    A/A*: inf  
    M: 0.200   |   P/Pt: 0.972    |    T/Tt: 0.992    |    A/A*: 2.964
    M: 0.400   |   P/Pt: 0.896    |    T/Tt: 0.969    |    A/A*: 1.590
    M: 0.600   |   P/Pt: 0.784    |    T/Tt: 0.933    |    A/A*: 1.188
    M: 0.800   |   P/Pt: 0.656    |    T/Tt: 0.887    |    A/A*: 1.038
    M: 1.000   |   P/Pt: 0.528    |    T/Tt: 0.833    |    A/A*: 1.000
    M: 1.200   |   P/Pt: 0.412    |    T/Tt: 0.776    |    A/A*: 1.030
    M: 1.400   |   P/Pt: 0.314    |    T/Tt: 0.718    |    A/A*: 1.115
    M: 1.600   |   P/Pt: 0.235    |    T/Tt: 0.661    |    A/A*: 1.250
    M: 1.800   |   P/Pt: 0.174    |    T/Tt: 0.607    |    A/A*: 1.439
    M: 2.000   |   P/Pt: 0.128    |    T/Tt: 0.556    |    A/A*: 1.688
    >>>
    """

    gamma = gas.gamma

    Mach_min = range[0]
    Mach_max = range[1]

    mach_nums = [i for i in np.arange(Mach_min,Mach_max+step,step)]
    t_list = [stagnation_temperature_ratio(mach=i, gas=gas)for i in mach_nums]
    p_list = [stagnation_pressure_ratio(mach=i, gas=gas) for i in mach_nums]
    a_list = [mach_area_star_ratio(mach=i, gas=gas) for i in mach_nums]
    rho_list = [stagnation_density_ratio(mach=i, gas=gas) for i in mach_nums]

    labl = '\u03B3 = ' + str(gamma)
    print("Isentropic Flow Parameters for " + gas.name + ", "+ labl)
    for index, num in enumerate(mach_nums):
        print('M: %0.3f' % num, '  |   P/Pt: %0.5f' % p_list[index], '   |    T/Tt: %0.4f' % t_list[index],  '   |    A/A*: %0.3f' % a_list[index],  '   |   rho/rho_t: %0.3f ' % rho_list[index])
    print("\n \n \n")



#==================================================
#mach_from_presure_ratio
#added fluid object
#==================================================    
def mach_from_pressure_ratio(pressure_initial: float, pressure_final: float, mach_initial: float, entropy=0, gas=air) -> float:
    """Return the Mach number given a Mach number and the local pressures
    
    Notes
    -----
    Given the local pressure in two regions and the Mach number in one,
    return the Mach number in the second region. Default arguments
    are for air and isentropic flow.

    Parameters
    ----------
    pressure_initial : `float`
        Pressure in region 1 \n
    pressure_final : `float`
        Pressure in region 2 \n
    mach_initial : `float`
        Mach number in region 1 \n
    entropy : `float`
        Change in entropy, if any \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The Mach number\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M2 = gd.mach_from_pressure_ratio(pressure_initial=10, pressure_final=2, mach_initial=1) 
    >>> M2
    2.1220079294384067
    >>>
    """

    gamma, R = gas.gamma, gas.R
    mach_final = (((pressure_initial/pressure_final * np.exp(entropy/R))**((gamma-1)/gamma) * (1 + (gamma-1)/2 * mach_initial**2) - 1) * 2/(gamma-1))**0.5
    return mach_final



#==================================================
#mach_from_temperature_ratio
#added fluid class
#==================================================
def mach_from_temperature_ratio(temperature_initial: float, temperature_final: float, mach_initial: float, gas=air) -> float:
    """Return the Mach number given a Mach number and two local temperatures
    
    Notes
    -----
    Given the local temperatures in two regions and the mach number in one,
    return the Mach number in the second region. Default fluid is air.

    Parameters
    ----------
    temperature_initial : `float`
        Temperature in region 1 \n
    temperature_final : `float`
        Temperature in region 2 \n
    mach_final : `float`
        Mach number in region 1 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The mach number\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_final = gd.mach_from_temperature_ratio(temperature_initial=300, temperature_final=150, mach_final=1)
    >>> mach_final
    2.6457513110645907
    >>>
    """

    gamma = gas.gamma
    mach_final = (( temperature_initial/temperature_final * (1 + (gamma-1)/2 * mach_initial**2) - 1) * 2/(gamma-1))**0.5
    return mach_final



#==================================================
#pressure_from_mach_ratio
#added fluid class
#==================================================
def pressure_from_mach_ratio(mach_initial: float, mach_final: float, pressure_initial: float, entropy=0, gas=air) -> float:
    """Return the pressure given a pressure in one region and the two Mach numbers
    
    Notes
    -----
    Given the Mach numbers in two regions and the pressure in one,
    return the missing pressure from the second region. Default arguments
    are for air and isentropic flow.
    
    Parameters
    ----------
    mach_initial : `float`
        Mach number in region 1 \n
    mach_final : `float`
        Mach number in region 2 \n
    pressure_initial : `float`
        Pressure in region 1 \n
    entropy : `float`
        Change in entropy, if any \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The local pressure\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> p_final = gd.pressure_from_mach_ratio(mach_initial=1, mach_final=2, pressure_initial=10)  
    >>> p_final
    2.4192491286747444
    >>>
    """

    gamma, R = gas.gamma, gas.R
    p_final = pressure_initial*((1+((gamma-1)/2)*mach_initial**2)/(1+((gamma-1)/2)*mach_final**2))**(gamma/(gamma-1))*np.exp(-entropy/R)
    return p_final



#==================================================
#temperature_from_mach_ratio
#added fluid class
#==================================================
def temperature_from_mach_ratio(mach_initial: float, mach_final: float, temperature_initial: float, gas=air) -> float:
    """Return the temperature given a temperature in one region and the two Mach numbers
    
    Notes
    -----
    Given the local Mach number in two regions and the temperature in one,
    return the missing temperature from the second region. Default fluid
    is air.

    Parameters
    ----------
    mach_initial : `float`
        Mach number in region 1 \n
    mach_final : `float`
        Mach number in region 2 \n
    temperature_initial : `float`
        Temperature in region 1 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The local temperature\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> T_final = gd.temperature_from_mach_ratio(mach_initial=1, mach_final=2, temperature_final=297.15) 
    >>> T_final
    198.10000000000002
    >>>
    """

    gamma = gas.gamma
    temperature_final = temperature_initial*(1+((gamma-1)/2)*mach_initial**2)/(1+((gamma-1)/2)*mach_final**2)
    return temperature_final



#==================================================
#entropy_produced
#added fluid class
#add examples
#==================================================
def entropy_produced(stagnation_pressure_initial: float, stagnation_pressure_final: float, gas=air) -> float:
    """Return the change in specific entropy from the stagnation pressure ratio
    
    Notes
    -----
    Given two stagnation pressures and the fluid, determine the entropy
    produced per unit mass.

    Parameters
    ----------
    stagnation_pressure_initial : `float`
        Stagnation pressure in region 1 \n
    stagnation_pressure_final : `float`
        Stagnation pressure in region 2 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The specific entropy\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.entropy_produced(pt_initial=10, pt_final=9, gas='air')
    30.238467993796142 #J / kg K
    >>> 
    """
    
    R = gas.R   
    ds = -R*np.log(stagnation_pressure_final/stagnation_pressure_initial)
    return ds



#==================================================
#mach_area_ratio_choked
#added fluid class 
#==================================================
def mach_area_star_ratio(mach: float, gas=air) -> float:
    """Returns the ratio of A / A* given the Mach number.
    
    Notes
    -----
    Given the Mach number and the ratio of specific heats, return the area
    ratio of the Mach number given to the area where Mach number is equal 
    to 1. Default fluid is air.
    
    Parameters
    ----------
    mach : `float`
        Mach Number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The ratio of area over choked area\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> A_Astar =gd.mach_area_ratio_choked(mach=3)
    >>> A_Astar
    4.23456790123457
    >>>
    """
    if mach == 0:
        return float('inf')

    gamma = gas.gamma
    a_star_ratio = 1/mach*((1+(gamma-1)/2*mach**2)/((gamma+1)/2))**((gamma+1)/(2*(gamma-1)))
    return a_star_ratio



#==================================================
#mach_area_ratio
#added fluid class
#==================================================
def mach_area_ratio(mach_initial: float, mach_final: float, gas=air, entropy=0) -> float:
    """Return the area ratio given the two Mach numbers

    Notes
    -----
    Given two mach numbers, return the area ratio required to accelerate
    or deaccelerate the flow accordingly. Default fluid is air.

    Parameters
    ----------
    mach_initial : `float` 
        Mach number in region 1 \n
    mach_final : `float`
        Mach number in region 2 \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    ds : `float`
        Entropy produced, if any \n
    
    Returns
    -------
    float
        The area ratio\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> A2_A1 = gd.mach_area_ratio(mach_initial=1.5, mach_final=2.5)
    >>> A2_A1
    2.241789331255894   #area ratio
    >>>
    """

    gamma, R = gas.gamma, gas.R
    A2_A1 = mach_initial/mach_final * ((1 + (gamma-1)/2 * mach_final**2 )/(1 + (gamma-1)/2 * mach_initial**2 ))**((gamma+1)/(2*(gamma-1))) * np.exp(entropy/R)
    return A2_A1



#==================================================
#mach_from_area_ratio

#==================================================    
def mach_from_area_ratio(area_ratio: float, gas=air) ->list:
    """Return the possible mach numbers given a choked area ratio A / A*
    
    Notes
    -----
    Given a ratio of area over an area where Mach = 1, return the subsonic and supersonic
    Mach numbers for the change area.

    Parameters
    ----------
    area_ratio : `float`
        The ratio of area over choked area \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    list
        The subsonic and supersonic mach numbers for the area ratio\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.mach_from_area_ratio(2)
    [0.30590383418910816, 2.197198121652187]
    >>>
    """

    def zero(mach, gas):
        return mach_area_star_ratio(mach=mach, gas=gas) - area_ratio

    subsonic = fsolve(zero, args=(gas), x0=0.1)
    supersonic = fsolve(zero, args=(gas), x0=5)
    sols = [subsonic[0], supersonic[0]]
    return sols



#==================================================
# mass_flux_max
# added fluid class
#==================================================
def mass_flux_max(stagnation_pressure: float, stagnation_temperature: float, gas=air) -> float:
    """Returns the maximum flow rate per unit choked area
    
    Notes
    -----
    Given stagnation pressure, stagnation temperature, and the fluid, 
    return the flow rate for a Mach number equal to 1. Default fluid 
    is air.

    **Units**:\n
    J / kg-K and Pa return kg/m^2 \n    
    kJ / kg-K and kPa returns kg/m^2 \n  
    ft-lbf / lbm-R and psi returns lbm/in^2 \n  
        
    Parameters
    ----------
    stagnation_pressure : `float`
        The stagnation pressure.\n
    stagnation_temperature : `float`
        The stagnation temperature.\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    metric : `bool`
        Use metric or US standard.\n
    
    Returns
    -------
    float
        The maximum mass flux\n
    
    Examples
    --------
    >>> mdot = 5 #kg/s
    >>> mdot_per_area = gd.choked_mdot(1000000, 300) #units are in Pascals
    >>> mdot_per_area
    2333.558560606226
    >>> throat_area = mdot / mdot_per_area
    >>> throat_area             #units are in meters squared
    0.0021426503214477164
    >>>
    #alternatively, we can use the english system; psi, rankine and
    get lbm/s/in^2
    >>> from gas_dynamics.fluids import air_us  
    >>> air_us.units
    'Btu / lbm-R'
    >>> flux = gd.mass_flux_max( stagnation_pressure=500, stagnation_temperature=500, gas=air_us)  
    >>> flux
    2.097208828890205
    >>>
    """

    gamma, R, gc = gas.gamma, gas.R, gas.gc
    
    mdot_a_star = (((gc*gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * stagnation_pressure/(stagnation_temperature**.5))
    return mdot_a_star



#==================================================
#mass_flux
#added fluid class
#==================================================
def mass_flux(mach: float, stagnation_pressure: float, stagnation_temperature: float, gas = air) -> float:
    """Determine mass flow rate for a mach number up to 1

    Notes
    -----
    Given stagnation pressure, stagnation temperature, and the fluid, 
    return the flow rate per unit area for the given Mach number. Default
    fluid is air.
    
    **Units**:\n
    J / kg-K and Pa return kg/s/m^2 \n    
    kJ / kg-K and kPa returns kg/s/m^2 \n  
    ft-lbf / lbm-R and psi returns lbm/s/in^2 \n  
    Btu / lbm-R and psi returns lbm/s/in^2 \n  
    
    Parameters
    ----------
    mach : `float`
        The mach number. Should not exceed 1 \n
    stagnation_pressure : `float`
        The stagnation pressure \n
    stagnation_temperature : `float`
        The stagnation temperature \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The mass flux\n

    Examples
    --------
    >>> #metric, input units are Pa, K, output is kg/s/m^2
    >>> flux = gd.mass_flux(mach=.8, stagnation_pressure=1e6, stagnation_temperature=500) 
    >>> flux
    1741.3113452036841
    >>> #us standard, input units are psi and Rankine, output units are lbm/s/in^2
    >>> from gas_dynamics.fluids import air_us  
    >>> air_us.units
    'Btu / lbm-R'
    >>> flux = gd.mass_flux(mach=.8, stagnation_pressure=500, stagnation_temperature=500, gas=air_us) 
    >>> flux
    2.01998480961849
    >>>
    """

    gamma, R, gc = gas.gamma, gas.R, gas.gc
    
    term1 = mach * (1 + mach**2 * (gamma-1)/2 )**(-(gamma+1)/(2*(gamma-1)))
    mass_flux = term1 * (gamma*gc/R)**.5 * stagnation_pressure/(stagnation_temperature**.5)
    return mass_flux



#==================================================
#plot_stagnation_ratios
#added fluid object, plot dark and light mode added
#added nice color scheme for lines.
#==================================================
def plot_stagnation_ratios(range=[.1,5], step=.01, gasses=[air, methane, argon], dark=True):
    """Plot the isentropic stagnation relationships for different gasses
    
    Notes
    -----
    Plots Mach number vs T/T, P/Pt, A/A*, rho/rho_t, for a list of
    specific heat ratios.
    
    Parameters
    ----------
    range : `list`
        The starting and ending Mach # in a list, ex: [.01,5] \n
    step : `float`
        The increment between min and max \n
    gasses : `list`
        A list of the user defined gas objects to be plotted \n
        ex: gasses = [air, methane, argon] \n
    dark : `bool`
        Use a dark mode plot. Default true.\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.plot_stagnation_ratios()
    >>>
    """

    if dark == True:
        plt.style.use('dark_background')
        gridcolor = 'w'
    else:
        plt.style.use('default')
        gridcolor = 'k'


    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    fig, axs = plt.subplots(2,2)
    title = "Isentropic Stagnation Relations"
    fig.suptitle(title)

    Mach_min = range[0]
    Mach_max = range[1]

    for n, f in enumerate(gasses):
        mach_nums = [i for i in np.arange(Mach_min,Mach_max+step,step)]
        t_list = [stagnation_temperature_ratio(mach=i, gas=f) for i in mach_nums]
        p_list = [stagnation_pressure_ratio(mach=i, gas=f) for i in mach_nums]
        a_list = [mach_area_star_ratio(mach=i, gas=f) for i in mach_nums]
        rho_list = [stagnation_density_ratio(mach=i, gas=f) for i in mach_nums]
        labl = '\u03B3 ' + f.name 

    #T/Tt
        axs[0,0].plot(mach_nums,t_list,label=labl, color=colors[n])
        axs[0,0].set_xlabel('Mach Number')
        axs[0,0].set_ylabel('T / Tt')
        axs[0,0].grid(b=True, which='major', color = gridcolor, linestyle = '--', alpha = .1)
        axs[0,0].grid(b=True, which='minor', color = gridcolor, linestyle = '--', alpha = .1)
        axs[0,0].minorticks_on()
        axs[0,0].legend()

    #P/Pt
        axs[0,1].plot(mach_nums,p_list,label=labl, color=colors[n])
        axs[0,1].set_xlabel('Mach Number')
        axs[0,1].set_ylabel('P / Pt')
        axs[0,1].grid(b=True, which='major', color = gridcolor, linestyle = '--', alpha = .1)
        axs[0,1].grid(b=True, which='minor', color = gridcolor, linestyle = '--', alpha = .1)
        axs[0,1].minorticks_on()
        axs[0,1].legend()
    #A/A*
        axs[1,0].plot(mach_nums,a_list,label=labl, color=colors[n])
        axs[1,0].set_xlabel('Mach Number')
        axs[1,0].set_ylabel('A / A*')
        axs[1,0].grid(b=True, which='major', color = gridcolor, linestyle = '--', alpha = .1)
        axs[1,0].grid(b=True, which='minor', color = gridcolor, linestyle = '--', alpha = .1)
        axs[1,0].minorticks_on()
        axs[1,0].legend()

    #rho/rho_t
        axs[1,1].plot(mach_nums,rho_list,label=labl, color=colors[n])
        axs[1,1].set_xlabel('Mach Number')
        axs[1,1].set_ylabel('rho / rho_t')
        axs[1,1].grid(b=True, which='major', color = gridcolor, linestyle = '--', alpha = .1)
        axs[1,1].grid(b=True, which='minor', color = gridcolor, linestyle = '--', alpha = .1)
        axs[1,1].minorticks_on()
        axs[1,1].legend()
    plt.show()