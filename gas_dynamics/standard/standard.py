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
def sonic_velocity(gas=air ,metric=True, T=273.15) -> float:
    """Returns the local speed of sound.
    
    Notes
    -----
    Given a ratio of specific heats, gas constant, and temperature
    this function returns the locoal speed of sound. Default fluid is air.
    
    Parameters
    ----------
    gas : `fluid`
        A user defined fluid object. Default is air \n
    metric : `bool`
        Use metric or US standard \n
    T : `float`
        The temperature \n

    Returns
    -------
    float
        The local speed of sound\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.sonic_velocity(air,T=500)
    448.2186966202994
    >>> gd.sonic_velocity(gas=Argon, T=300)

    """

    if metric == False:
        gc = 32.174
    else:
        gc = 1
    gamma, R = gas.gamma, gas.R
    a = ( gc * gamma * R *T )**.5
    return a



#==================================================
#stagnation_pressure
#implemented fluid class and string for output returned
#==================================================    
def stagnation_pressure(pt=None, M=None, p=None, gas=air, output=False) -> float:
    """Returns the stagnation pressure given pressure and Mach number.

    Notes
    -----
    Given a pressure, Mach number, and a ratio of specific heats return
    the stagnation pressure. Alternatively, provided two arguments
    the function will return the missing one. Default fluid is air.

    Parameters
    ----------
    pt : `float`
        The stagnation pressure.\n
    p : `float`
        The pressure.\n
    M : `float`
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
    >>> pt = gd.stagnation_pressure(p=10, M=1)
    >>> pt
    18.92929158737854
    >>> M = gd.stagnation_pressure(p=10, pt=pt)
    >>> M
    1.0
    >>>
    """

    gamma = gas.gamma
    if pt == None:
        pt = p* ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        if output == True:
            print('Returned stagnation pressure')
        return pt
   
    if M == None:
        M = (((pt/p)**((gamma-1)/gamma) -1 ) * 2/(gamma-1) ) ** .5
        if output == True:
            print('Returned Mach')
        return M

    if p == None:
        p = pt / ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        if output == True:
            print('Returned pressure')
        return p



#==================================================
#stagnation_temperature
#implemented output string, fluid class, 
#==================================================    
def stagnation_temperature(T=None, Tt=None , M=None, gas=air, output=False) -> float :
    """Returns the stagnation temperature given temperature and Mach number.
    
    Notes
    -----
    Given a temperature, Mach number, and a ratio of specific heats 
    this function returns the stagnation temperature. Alternatively,
    provided two arguments the function will return the missing one.
    Default fluid is air.

    Parameters
    ----------
    Tt : `float`
        The stagnation temperature\n
    T : `float`
        The temperature\n
    M : `float`
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
    >>> Tt = gd.stagnation_temperature(T=300, M=1)
    >>> Tt
    360.0
    >>> M = gd.stagnation_temperature(T=300, Tt=Tt)
    >>> M 
    1.0
    >>>
    """

    gamma = gas.gamma
    if Tt == None:
        Tt = T * ( 1 + (gamma-1)/2 * M**2)
        if output == True:
            print('Returned stagnation temperature')
        return Tt

    if M == None:
        M = ((Tt /T - 1) * 2/(gamma-1))**.5
        if output == True:
            print('Returned Mach')
        return M
        
    if T == None:
        T = Tt/( 1 + (gamma-1)/2 * M**2)
        if output == True:
            print('Returned temperature')
        return T



#==================================================
#stagnation_pressure_ratio
#added fluid class, removed unnnecessary metric argument
#==================================================
def stagnation_pressure_ratio(M: float, gas=air) -> float:
    """Returns the pressure ratio of p / p_t
    
    Notes
    -----
    Given a Mach number and ratio of specific heats return the relation of
    pressure over stagnation pressure. Default fluid is air.

    Parameters
    ----------
    M : `float`
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
    >>> p_pt = gd.stagnation_pressure_ratio(M=3)
    >>> p_pt
    0.027223683703862817
    >>>
    """

    gamma = gas.gamma
    denom = 1 + (gamma-1)/2 * M**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio



#==================================================
#stagnation_temperature_ratio
#added fluid class
#==================================================
def stagnation_temperature_ratio(M: float, gas=air) -> float:
    """Returns the temperature ratio of T / T_t
    
    Notes
    -----
    Given a Mach number and ratio of specific heats return the relation of
    temperature over stagnation temperature.  Default fluid is air.

    Parameters
    ----------
    M : `float`
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
    >>> T_Tt = gd.stagnation_temperature_ratio(M=1.5)
    >>> T_Tt
    0.6896551724137931
    >>>
    """

    gamma = gas.gamma
    Tt_ratio = 1 / (1+(gamma-1)/2 *M**2)
    return Tt_ratio



#==================================================
#stagnation_density_ratio
#added fluid class
#==================================================
def stagnation_density_ratio(M: float, gas=air) -> float:
    """Returns the density ratio rho / rho_t
    
    Notes
    -----
    Given a Mach number and ratio of specific heats, return the relation
    of density over stagnation density. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
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
    >>> rho_rho_t = gd.stagnation_density_ratio(M=1.5) 
    >>> rho_rho_t
    0.39498444639115327
    >>>
    """

    gamma = gas.gamma
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M**2 )) ** (1 / (gamma-1))
    return rho_t_ratio



#==================================================
#stagnation ratio 
#TODO: make me
#added fluid object
#==================================================
def stagnation_ratio(M: float, gas=air) -> list:
    """Return stagnation pressure, temperature, density, and choked area ratio for a mach number
        
    Notes
    -----
    Given a mach number and the fluid, return the three stagnation ratios
    and the ratio of the area to the choked area. Default fluid is air.

    Parameters
    ----------
    M : `float`
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
    t_list = [stagnation_temperature_ratio(M=i, gas=gas)for i in mach_nums]
    p_list = [stagnation_pressure_ratio(M=i, gas=gas) for i in mach_nums]
    a_list = [mach_area_star_ratio(M=i, gas=gas) for i in mach_nums]
    rho_list = [stagnation_density_ratio(M=i, gas=gas) for i in mach_nums]

    labl = '\u03B3 = ' + str(gamma)
    print("Isentropic Flow Parameters for " + gas.name + ", "+ labl)
    for index, num in enumerate(mach_nums):
        print('M: %0.3f' % num, '  |   P/Pt: %0.3f' % p_list[index], '   |    T/Tt: %0.3f' % t_list[index],  '   |    A/A*: %0.3f' % a_list[index],  '   |   rho/rho_t: %0.3f ' % rho_list[index])
    print("\n \n \n")



#==================================================
#mach_from_presure_ratio
#added fluid object
#==================================================    
def mach_from_pressure_ratio(p1: float, p2: float, M1: float, ds=0, gas=air) -> float:
    """Return the Mach number given a Mach number and the local pressures
    
    Notes
    -----
    Given the local pressure in two regions and the Mach number in one,
    return the Mach number in the second region. Default arguments
    are for air and isentropic flow.

    Parameters
    ----------
    p1 : `float`
        Pressure in region 1 \n
    p2 : `float`
        Pressure in region 2 \n
    M1 : `float`
        Mach number in region 1 \n
    ds : `float`
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
    >>> M2 = gd.mach_from_pressure_ratio(p1=10, p2=2, M1=1) 
    >>> M2
    2.1220079294384067
    >>>
    """

    gamma, R = gas.gamma, gas.R
    M2 = (((p1/p2 * np.exp(ds/R))**((gamma-1)/gamma) * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
    return M2



#==================================================
#mach_from_temperature_ratio
#added fluid class
#==================================================
def mach_from_temperature_ratio(T1: float, T2: float, M1: float, gas=air) -> float:
    """Return the Mach number given a Mach number and two local temperatures
    
    Notes
    -----
    Given the local temperatures in two regions and the mach number in one,
    return the Mach number in the second region. Default fluid is air.

    Parameters
    ----------
    T1 : `float`
        Temperature in region 1 \n
    T2 : `float`
        Temperature in region 2 \n
    M1 : `float`
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
    >>> M2 = gd.mach_from_temperature_ratio(T1=300, T2=150, M1=1)
    >>> M2
    2.6457513110645907
    >>>
    """

    gamma = gas.gamma
    M2 = (( T1/T2 * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
    return M2



#==================================================
#pressure_from_mach_ratio
#added fluid class
#==================================================
def pressure_from_mach_ratio(M1: float, M2: float, p1: float, ds=0, gas=air) -> float:
    """Return the pressure given a pressure in one region and the two Mach numbers
    
    Notes
    -----
    Given the Mach numbers in two regions and the pressure in one,
    return the missing pressure from the second region. Default arguments
    are for air and isentropic flow.
    
    Parameters
    ----------
    M1 : `float`
        Mach number in region 1 \n
    M2 : `float`
        Mach number in region 2 \n
    p1 : `float`
        Pressure in region 1 \n
    ds : `float`
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
    >>> p2 = gd.pressure_from_mach_ratio(M1=1, M2=2, p1=10)  
    >>> p2
    2.4192491286747444
    >>>
    """

    gamma, R = gas.gamma, gas.R
    p2 = p1 * ((1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2))**(gamma/(gamma-1)) * np.exp(-ds/R)
    return p2



#==================================================
#temperature_from_mach_ratio
#added fluid class
#==================================================
def temperature_from_mach_ratio(M1: float, M2: float, T1: float, gas=air) -> float:
    """Return the temperature given a temperature in one region and the two Mach numbers
    
    Notes
    -----
    Given the local Mach number in two regions and the temperature in one,
    return the missing temperature from the second region. Default fluid
    is air.

    Parameters
    ----------
    M1 : `float`
        Mach number in region 1 \n
    M2 : `float`
        Mach number in region 2 \n
    T1 : `float`
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
    >>> T2 = gd.temperature_from_mach_ratio(M1=1, M2=2, T1=297.15) 
    >>> T2
    198.10000000000002
    >>>
    """

    gamma = gas.gamma
    T2 = T1 * (1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2)
    return T2



#==================================================
#entropy_produced
#added fluid class
#add examples
#==================================================
def entropy_produced(pt1: float, pt2: float, gas=air) -> float:
    """Return the change in specific entropy from the stagnation pressure ratio
    
    Notes
    -----
    Given two stagnation pressures and the fluid, determine the entropy
    produced per unit mass.

    Parameters
    ----------
    pt1 : `float`
        Stagnation pressure in region 1 \n
    pt2 : `float`
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
    >>> gd.entropy_produced(pt1=10, pt2=9, gas='air')
    30.238467993796142 #J / kg K
    >>> 
    """
    
    R = gas.R   
    ds = -R*np.log(pt2/pt1)
    return ds



#==================================================
#mach_area_ratio_choked
#added fluid class 
#==================================================
def mach_area_star_ratio(M: float, gas=air) -> float:
    """Returns the ratio of A / A* given the Mach number.
    
    Notes
    -----
    Given the Mach number and the ratio of specific heats, return the area
    ratio of the Mach number given to the area where Mach number is equal 
    to 1. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
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
    >>> A_Astar =gd.mach_area_ratio_choked(M=3)
    >>> A_Astar
    4.23456790123457
    >>>
    """
    if M == 0:
        return float('inf')

    gamma = gas.gamma
    a_star_ratio = 1/M * ((1 + (gamma-1)/2 * M**2) / ((gamma+1)/2)) ** ((gamma+1)/(2*(gamma-1)))
    return a_star_ratio



#==================================================
#mach_area_ratio
#added fluid class
#==================================================
def mach_area_ratio(M1: float, M2: float, gas=air, ds=0) -> float:
    """Return the area ratio given the two Mach numbers

    Notes
    -----
    Given two mach numbers, return the area ratio required to accelerate
    or deaccelerate the flow accordingly. Default fluid is air.

    Parameters
    ----------
    M1 : `float` 
        Mach number in region 1 \n
    M2 : `float`
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
    >>> A2_A1 = gd.mach_area_ratio(M1=1.5, M2=2.5)
    >>> A2_A1
    2.241789331255894   #area ratio
    >>>
    """

    gamma, R = gas.gamma, gas.R
    A2_A1 = M1/M2 * ((1 + (gamma-1)/2 * M2**2 )/(1 + (gamma-1)/2 * M1**2 ))**((gamma+1)/(2*(gamma-1))) * np.exp(ds/R)
    return A2_A1



#==================================================
#mach_from_area_ratio

#==================================================    
def  mach_from_area_ratio(a_ratio: float, gas=air) ->list:
    """Return the possible mach numbers given a choked area ratio A / A*
    
    Notes
    -----
    Given a ratio of area over an area where Mach = 1, return the subsonic and supersonic
    Mach numbers for the change area.

    Parameters
    ----------
    a_ratio : `float`
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

    def zero(M, gas):
        return mach_area_star_ratio(M=M, gas=gas) - a_ratio

    subsonic = fsolve(zero, args=(gas), x0=0.5)
    supersonic = fsolve(zero, args=(gas), x0=5)
    sols = [subsonic[0], supersonic[0]]
    return sols



#==================================================
# mass_flux_max
# added fluid class
#==================================================
def mass_flux_max(pt: float, Tt: float, gas=air, metric=True) -> float:
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
    pt : `float`
        The stagnation pressure.\n
    Tt : `float`
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
    >>> from gas_dynamics.extra import air_us  
    >>> air_us.units
    'Btu / lbm-R'
    >>> flux = gd.mass_flux_max( pt=500, Tt=500, gas=air_us)  
    >>> flux
    2.097208828890205
    >>>
    """

    gamma, R = gas.gamma, gas.R
    if metric == False:
        gc = 32.174 #lbm -ft / lbf -s^2
    else:
        gc = 1
    
    mdot_a_star = (((gc*gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * pt/(Tt**.5))
    return mdot_a_star



#==================================================
#mass_flux
#added fluid class
#==================================================
def mass_flux(M: float, pt: float, Tt: float, gas = air, metric = True) -> float:
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
    M : `float`
        The mach number. Should not exceed 1 \n
    pt : `float`
        The stagnation pressure \n
    Tt : `float`
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
    >>> flux = gd.mass_flux(M=.8, pt=1e6, Tt=500) 
    >>> flux
    1741.3113452036841
    >>> #us standard, input units are psi and Rankine, output units are lbm/s/in^2
    >>> from gas_dynamics.extra import air_us  
    >>> air_us.units
    'Btu / lbm-R'
    >>> flux = gd.mass_flux(M=.8, pt=500, Tt=500, gas=air_us) 
    >>> flux
    2.01998480961849
    >>>
    """

    gamma, R = gas.gamma, gas.R
    if metric == False:
        gc = 32.174 #lbm-ft / lbf s^2
    else:
        gc = 1
    
    term1 = M * (1 + M**2 * (gamma-1)/2 )**(-(gamma+1)/(2*(gamma-1)))
    mass_flux = term1 * (gamma*gc/R)**.5 * pt/(Tt**.5)
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
        t_list = [stagnation_temperature_ratio(M=i, gas=f) for i in mach_nums]
        p_list = [stagnation_pressure_ratio(M=i, gas=f) for i in mach_nums]
        a_list = [mach_area_star_ratio(M=i, gas=f) for i in mach_nums]
        rho_list = [stagnation_density_ratio(M=i, gas=f) for i in mach_nums]
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