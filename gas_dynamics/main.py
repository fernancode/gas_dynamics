"""Equatons, plots, and tables for solving compresible flow problems.
Included are functions to solve problems relating to compressible flow,
from stagnation relations to determining the mach number from changes in
local properties. Tables can be made for any gas and its respective 
ratio of specific heats, as well as plots and charts for shock and
mach number relationships.
  Typical usage example:
  determine the mach number after a normal shock, given a starting
  mach number of 2
  >>> import gas_dynamics as gd
  >>> mach_before_shock = 2
  >>> mach_after_shock = gd.shock_mach(M1=mach_before_shock)
  >>> mach_after_shock
  0.5773502691896257
  >>>
  Generate isentropic flow stagnation relations for methane
  >>> import gas_dynamics as gd
  >>> gd.stagnation_ratios(range=[0,2], inc=.2, gas='methane')
  M: 0.000   |   P/Pt: 1.000    |    T/Tt: 1.000    |    A/A*: inf    |   rho/rho_t: 1.000
  M: 0.200   |   P/Pt: 0.974    |    T/Tt: 0.994    |    A/A*: 2.988    |   rho/rho_t: 0.980
  M: 0.400   |   P/Pt: 0.901    |    T/Tt: 0.975    |    A/A*: 1.600    |   rho/rho_t: 0.924
  M: 0.600   |   P/Pt: 0.794    |    T/Tt: 0.946    |    A/A*: 1.192    |   rho/rho_t: 0.839
  M: 0.800   |   P/Pt: 0.669    |    T/Tt: 0.907    |    A/A*: 1.039    |   rho/rho_t: 0.737
  M: 1.000   |   P/Pt: 0.542    |    T/Tt: 0.862    |    A/A*: 1.000    |   rho/rho_t: 0.629
  M: 1.200   |   P/Pt: 0.425    |    T/Tt: 0.813    |    A/A*: 1.032    |   rho/rho_t: 0.523
  M: 1.400   |   P/Pt: 0.325    |    T/Tt: 0.761    |    A/A*: 1.121    |   rho/rho_t: 0.426
  M: 1.600   |   P/Pt: 0.243    |    T/Tt: 0.709    |    A/A*: 1.267    |   rho/rho_t: 0.342
  M: 1.800   |   P/Pt: 0.179    |    T/Tt: 0.659    |    A/A*: 1.474    |   rho/rho_t: 0.271
  M: 2.000   |   P/Pt: 0.130    |    T/Tt: 0.610    |    A/A*: 1.754    |   rho/rho_t: 0.213
Copyright 2020 by Fernando A de la Fuente
All rights reserved
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

#GOOD
def fluid(fluid=[], metric=True, R=[], gamma=[]):
    """Return the ratio of specific heats and gas constant for the fluid
    Description
    -----------
    Return the ratio of specific heats and gas constant for the fluid
    Parameters
    ----------
    fluid : `string`
        The fluid\n
    metric : `bool`
        Use Metric or US Standard \n
    R : `float`
        Gas constant, if known \n
    gamma : `gamma`
        Ratio of specific heats, if known \n
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> R, gamma = gd.fluid(fluid='Nitrogen') 
    >>> R, gamma
    (1.4, 296)
    >>>       
    """

    if R or gamma != []:
        return gamma, R
    
    if fluid == []:
        gamma, R = 1.4, 287
        return gamma, R

    fluid = fluid.lower()
    if metric == True:
        #N m / kg K
        if fluid == 'air':
            gamma, R = 1.4, 287
        elif fluid == 'argon':
            gamma, R = 1.67, 208
        elif fluid == 'carbon dioxide':
            gamma, R = 1.29, 189
        elif fluid == 'carbon monoxide':
            gamma, R = 1.4, 297
        elif fluid == 'helium':
            gamma, R = 1.67, 2080
        elif fluid == 'hydrogen':
            gamma, R = 1.41, 4120
        elif fluid == 'methane':
            gamma, R = 1.32, 519
        elif fluid == 'nitrogen':
            gamma, R = 1.4, 296
        elif fluid == 'oxygen':
            gamma, R = 1.4, 260
        elif fluid == 'water' or fluid == 'steam':
            gamma, R = 1.33, 461
        elif fluid == 'custom':
            gamma = float(input('gamma: '))
            R = float(input('R (J/kg K): '))
        else:
            gamma, R = 1.4, 287
        return gamma, R

    if metric == False:
        # ft lbf / lbm R
        if fluid == 'air':
            gamma, R = 1.4, 53.3
        elif fluid == 'argon':
            gamma, R = 1.67, 38.7
        elif fluid == 'carbon dioxide':
            gamma, R = 1.29, 35.1
        elif fluid == 'carbon monoxide':
            gamma, R = 1.4, 55.2
        elif fluid == 'helium':
            gamma, R = 1.67, 386
        elif fluid == 'hydrogen':
            gamma, R = 1.41, 766
        elif fluid == 'methane':
            gamma, R = 1.32, 96.4
        elif fluid == 'nitrogen':
            gamma, R = 1.4, 55.1
        elif fluid == 'oxygen':
            gamma, R = 1.4, 48.3
        elif fluid == 'water' or fluid == 'steam':
            gamma, R = 133, 85.7
        else:
            gamma, R = 1.4, 55.3
        return gamma, R

#GOOD
def plot_stagnation_ratios(range=[.1,5],inc=.01, gasses=['air','methane','argon']):
    """Plot the isentropic stagnation relationships for different gasses
    
    Description
    -----------
    Plots Mach number vs T/T, P/Pt, A/A*, rho/rho_t, for a list of
    specific heat ratios.
    
    Parameters
    ----------
    range : `list`
        The starting and ending Mach # in a list, ex: [.01,5] \n
    inc : `float`
        The increment between min and max \n
    gasses : `list`
        A list of the gasses to be plotted, ex ['air','methane','argon'] \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.plot_stagnation_ratios()
    """

    plt.style.use('dark_background')
    fig, axs = plt.subplots(2,2)
    title = "Isentropic Stagnation Relations"
    fig.suptitle(title)

    Mach_min = range[0]
    Mach_max = range[1]

    for f in gasses:
        mach_nums = [i for i in np.arange(Mach_min,Mach_max+inc,inc)]
        t_list = [stagnation_temperature_ratio(M=i, gas=f) for i in mach_nums]
        p_list = [stagnation_pressure_ratio(M=i, gas=f) for i in mach_nums]
        a_list = [mach_area_choked_ratio(M=i, gas=f) for i in mach_nums]
        rho_list = [stagnation_density_ratio(M=i, gas=f) for i in mach_nums]
        labl = '\u03B3 ' + f #str(g)

    #T/Tt
        axs[0,0].plot(mach_nums,t_list,label=labl)
        axs[0,0].set_xlabel('Mach Number')
        axs[0,0].set_ylabel('T / Tt')
        axs[0,0].grid(b=True, which='major', color = 'w', linestyle = '--', alpha = .1)
        axs[0,0].grid(b=True, which='minor', color = 'w', linestyle = '--', alpha = .1)
        axs[0,0].minorticks_on()
        axs[0,0].legend()

    #P/Pt
        axs[0,1].plot(mach_nums,p_list,label=labl)
        axs[0,1].set_xlabel('Mach Number')
        axs[0,1].set_ylabel('P / Pt')
        axs[0,1].grid(b=True, which='major', color = 'w', linestyle = '--', alpha = .1)
        axs[0,1].grid(b=True, which='minor', color = 'w', linestyle = '--', alpha = .1)
        axs[0,1].minorticks_on()
        axs[0,1].legend()
    #A/A*
        axs[1,0].plot(mach_nums,a_list,label=labl)
        axs[1,0].set_xlabel('Mach Number')
        axs[1,0].set_ylabel('A / A*')
        axs[1,0].grid(b=True, which='major', color = 'w', linestyle = '--', alpha = .1)
        axs[1,0].grid(b=True, which='minor', color = 'w', linestyle = '--', alpha = .1)
        axs[1,0].minorticks_on()
        axs[1,0].legend()

    #rho/rho_t
        axs[1,1].plot(mach_nums,rho_list,label=labl)
        axs[1,1].set_xlabel('Mach Number')
        axs[1,1].set_ylabel('rho / rho_t')
        axs[1,1].grid(b=True, which='major', color = 'w', linestyle = '--', alpha = .1)
        axs[1,1].grid(b=True, which='minor', color = 'w', linestyle = '--', alpha = .1)
        axs[1,1].minorticks_on()
        axs[1,1].legend()
    plt.show()

#GOOD
def stagnation_ratios(range=[0,5], inc=.1, gas='air') -> str:
    """Returns the isentropic flow tables in the given range.
    
    Description
    -----------
    Given a ratio of specific heats, print out the stagnation
    temperature ratio, stagnation pressure ratio, the area to 
    choked area ratio, and the stagnation density ratio for every
    incremental Mach number.
    
    Parameters
    ----------
    range : `list`
        The starting and ending Mach numbers in a list, ex: [.01,5] \n
    inc : `float`
        The increment between min and max \n
    gas : `str`
        The fluid \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.stagnation_ratios(range=[0,2], inc=.2, gas='nitrogen')  
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

    gamma, R = fluid(gas)

    Mach_min = range[0]
    Mach_max = range[1]

    mach_nums = [i for i in np.arange(Mach_min,Mach_max+inc,inc)]
    t_list = [stagnation_temperature_ratio(M=i, gas=gas)for i in mach_nums]
    p_list = [stagnation_pressure_ratio(M=i, gas=gas) for i in mach_nums]
    a_list = [mach_area_choked_ratio(M=i, gas=gas) for i in mach_nums]
    rho_list = [stagnation_density_ratio(M=i, gas=gas) for i in mach_nums]

    labl = '\u03B3 = ' + str(gamma)
    print("Isentropic Flow Parameters for " + labl)
    for index, num in enumerate(mach_nums):
        print('M: %0.3f' % num, '  |   P/Pt: %0.3f' % p_list[index], '   |    T/Tt: %0.3f' % t_list[index],  '   |    A/A*: %0.3f' % a_list[index],  '   |   rho/rho_t: %0.3f ' % rho_list[index])
    print("\n \n \n")

#GOOD
def sonic_velocity(gas='air' ,metric=True, T=273.15) -> float:
    """Returns the local speed of sound.
    
    Description
    -----------
    Given a ratio of specific heats, gas constant, and temperature,
    this function returns the locoal speed of sound. Default fluid is air.
    
    Parameters
    ----------
    gas: `str`
        The fluid. \n
    T : `float`
        The temperature \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.sonic_velocity('air',T=500)
    448.2186966202994
    >>> gd.sonic_velocity(gas='argon', T=300)
    """

    gamma, R = fluid(gas, metric)
    a = (gamma*R*T)**.5
    return a

#GOOD
def entropy_produced(pt1: float, pt2: float, gas='air', metric=True) -> float:
    """Return the change in specific entropy from the stagnation pressure ratio
    
    Description
    -----------
    Given two stagnation pressures and the fluid, determine the entropy
    produced per unit mass
    
    Parameters
    ----------
    pt1 : `float`
        Stagnation pressure in region 1 \n
    pt2 : `float`
        Stagnation pressure in region 2 \n
    gas : `float`
        The fluid \n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.entropy_produced(pt1=10, pt2=9, gas='air')
    30.238467993796142 #J / kg K
    >>> 
    """
    
    gamma, R = fluid(gas, metric)   
    ds = -R*np.log(pt2/pt1)
    return ds

#GOOD
def pressure_from_mach_ratio(M1: float, M2: float, p1: float, ds=0, gas='air', metric=True) -> float:
    """Return the pressure given a starting pressure and the two Mach numbers
    
    Description
    -----------
    Given the Mach numbers in two regions and the pressure in one,
    return the missing pressure from the second region. Default arguments
    are for air and isentropic flow.
    
    Parameters:
    ----------
    M1 : `float`
        Mach number in region 1 \n
    M2 : `float`
        Mach number in region 2 \n
    p1 : `float`
        Pressure in region 1 \n
    ds : `float`
        Change in entropy, if any \n
    gas : `str`
        The fluid \n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> p2 = gd.pressure_from_mach_ratio(M1=1, M2=2, p1=10)  
    >>> p2
    2.4192491286747444
    >>>
    """

    gamma, R = fluid(gas, metric)
    p2 = p1 * ((1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2))**(gamma/(gamma-1)) * np.exp(-ds/R)
    return p2

#GOOD
def mach_from_pressure_ratio(p1: float, p2: float, M1: float, ds=0, gas='air', metric=True) -> float:
    """Return the Mach number given a Mach number and the local pressures
    
    Description
    -----------
    Given the local pressure in two regions and the Mach number in one,
    return the missing Machnumber# in the second region. Default arguments
    are for air and isentropic flow.
    
    Parameters:
    ----------
    p1 : `float`
        Pressure in region 1 \n
    p2 : `float`
        Pressure in region 2 \n
    M1 : `float`
        Mach number in region 1 \n
    ds : `float`
        Change in entropy, if any \n
    gas : `str`
        The fluid \n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M2 = gd.mach_from_pressure_ratio(p1=10, p2=2, M1=1) 
    >>> M2
    2.1220079294384067
    >>>
    """

    gamma, R = fluid(gas, metric)
    M2 = (((p1/p2 * np.exp(ds/R))**((gamma-1)/gamma) * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
    return M2

#GOOD
def temperature_from_mach_ratio(M1: float, M2: float, T1: float, gas='air', metric=True) -> float:
    """Return the temperature given a temperature and the two Mach numbers
    
    Description
    -----------
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
    gas : `str` 
        The fluid \n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> T2 = gd.temperature_from_mach_ratio(M1=1, M2=2, T1=297.15) 
    >>> T2
    198.10000000000002
    >>>
    """

    gamma, R = fluid(gas, metric)
    T2 = T1 * (1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2)
    return T2

#GOOD
def mach_from_temperature_ratio(T1: float, T2: float, M1: float, gas='air', metric=True) -> float:
    """Return the Mach number given a Mach number and two local temperatures
    
    Description
    -----------
    Given the local temperatures in two regions and the mach number in one,
    return the missing mach number from the second region. Default fluid
    is air.
    
    Parameters
    ----------
    T1 : `float`
        Temperature in region 1 \n
    T2 : `float`
        Temperature in region 2 \n
    M1 : `float`
        Mach number in region 1 \n
    gas : `float`
        The fluid \n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M2 = gd.mach_from_temperature_ratio(T1=300, T2=150, M1=1)
    >>> M2
    2.6457513110645907
    >>>
    """

    gamma, R = fluid(gas, metric)
    M2 = (( T1/T2 * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
    return M2

#GOOD
def mach_area_ratio(M1: float, M2: float, gas='air', ds=0, metric=True) -> float:
    """Return the area ratio given the two Mach numberss
    Description
    -----------
    Given two mach numbers, return the area ratio required to accelerate
    or deaccelerate the flow accordingly. Default fluid is air.
    
    Parameters
    ----------
    M1 : `float` 
        Mach number in region 1 \n
    M2 : `float`
        Mach number in region 2 \n
    gas : `str`
        The Fluid\n
    ds : `float`
        Entropy produced, if any \n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> A2_A1 = gd.mach_area_ratio(M1=1.5, M2=2.5)
    >>> A2_A1
    2.241789331255894   #area ratio
    >>>
    """

    gamma, R = fluid(gas, metric)
    A2_A1 = M1/M2 * ((1 + (gamma-1)/2 * M2**2 )/(1 + (gamma-1)/2 * M1**2 ))**((gamma+1)/(2*(gamma-1))) * np.exp(ds/R)
    return A2_A1

#GOOD
def mach_area_choked_ratio(M: float, gas='air', metric=True) -> float:
    """Returns the ratio of A / A* given the Mach number.
    
    Description
    -----------
    Given the Mach number and the ratio of specific heats, return the area
    ratio of the Mach number given to Mach number equal to 1. Default fluid
    is air.
    
    Parameters
    ----------
    M : `float`
        Mach Number \n
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> A_Astar =gd.mach_area_choked_ratio(M=3)
    >>> A_Astar
    4.23456790123457
    >>>
    """

    gamma, R = fluid(gas, metric)
    a_star_ratio = 1/M * ((1 + (gamma-1)/2 * M**2) / ((gamma+1)/2)) ** ((gamma+1)/(2*(gamma-1)))
    return a_star_ratio

#GOOD
def stagnation_pressure(pt=None, M=None, p=None, gas='air', metric=True) -> float:
    """Returns the stagnation pressure given pressure and Mach number.
    Description
    -----------
    Given a pressure, Mach number, and a ratio of specific heats, return
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
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
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

    gamma, R = fluid(gas, metric)
    if pt == []:
        pt = p* ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        return pt
   
    if M == []:
        M = (((pt/p)**((gamma-1)/gamma) -1 ) * 2/(gamma-1) ) ** .5
        return M

    if p == []:
        p = pt / ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        return p

#GOOD
def stagnation_pressure_ratio(M: float, gas='air', metric=True) -> float:
    """Returns the pressure ratio of p / p_t
    
    Description
    -----------
    Given a Mach number and ratio of specific heats, return the relation of
    p / p_t. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number\n
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> p_pt = gd.stagnation_pressure_ratio(M=3)
    >>> p_pt
    0.027223683703862817
    >>>
    """

    gamma, R = fluid(gas, metric)
    denom = 1 + (gamma-1)/2 * M**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio

#GOOD
def stagnation_temperature(T=None, Tt=None , M=None, gas='air', metric=True) -> float :
    """Returns the stagnation temperature given temperature and Mach number.
    
    Description
    -----------
    Given a temperature, Mach number, and a ratio of specific heats, 
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
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
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

    gamma, R = fluid(gas, metric)
    if Tt == None:
        Tt = T * ( 1 + (gamma-1)/2 * M**2)
        return Tt

    if M == None:
        M = ((Tt /T - 1) * 2/(gamma-1))**.5
        return M
        
    if T == None:
        T = Tt/( 1 + (gamma-1)/2 * M**2)
        return T

#TODO: ADD EXAMPLES
def stagnation_temperature_ratio(M: float, gas='air', metric=True) -> float:
    """Returns the temperature ratio of T / T_t
    
    Description
    -----------
    Given a Mach number and ratio of specific heats, return the relation of
    T / T_t.  Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number\n
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    """

    gamma, R = fluid(gas, metric)
    Tt_ratio = 1 / (1+(gamma-1)/2 *M**2)
    return Tt_ratio

#TODO: add examples
def stagnation_density_ratio(M: float, gas='air', metric=True) -> float:
    """Returns the density ratio of rho/ rho_t
    
    Description
    -----------
    Given a Mach number and ratio of specific heats, return the relation
    of rho / rho_t. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach # \n
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    Examples
    --------
    
    """

    gamma, R = fluid(gas, metric)
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M**2 )) ** (1 / (gamma-1))
    return rho_t_ratio

#GOOD
def choked_mdot(pt: float, Tt: float, gas='air', metric=True) -> float:
    """Returns the maximum flow rate per unit choked area
    
    Description
    -----------
    Given stagnation pressure, stagnation temperature, and the fluid, 
    return the flow rate per unit choked area. Default fluid is air. 
    
    Check your units! metric units need to be in Pa \n
    #TODO: figure out what the std units output is
    
    Parameters
    ----------
    pt : `float`
        The stagnation pressure.\n
    Tt : `float`
        The stagnation temperature.\n
    gas : `str`
        The Fluid.\n
    metric : `bool`
        Use metric or US standard.\n
    
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
    """

    gamma, R = fluid(gas, metric)
    if metric == True:
        mdot_a_star = (((gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * pt/(Tt**.5))
        return mdot_a_star

    if metric == False:
        gc = 32.174 #lbm -ft / lbf -s^2
        mdot_a_star = (((gc*gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * pt/(Tt**.5))
        return mdot_a_star

#GOOD:
def shock_tables(range=[1,5], inc=.01, gas='air', metric=True) -> str:
    """Returns shock tables for a range of Mach numberss.
    
    Description
    -----------
    Given a range of Mach numberss and a ratio of specific heats, generate
    the standing normal shock tables for every incremental Mach number
    in between.
    
    Parameters
    ----------
    range : `list`
        The starting and ending Mach # in a list, ie: [1,5]. \n
    inc : `float`
        The step size for the tables. \n
    gas : `str`
        The fluid \n
    metric : `bool`
        Use metric or US standard \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.shock_tables(range=[1,2], inc=.1)
    Normal Shock Parameters for γ = 1.4
    M: 1.00   |   M2: 1.0000   |    p2/p1: 1.0000   |    T2/T1: 1.0000   |   dV/a: 0.0000   |   pt2/pt1: 1.000000
    M: 1.10   |   M2: 0.9118   |    p2/p1: 1.2450   |    T2/T1: 1.0649   |   dV/a: 0.1591   |   pt2/pt1: 0.998928
    M: 1.20   |   M2: 0.8422   |    p2/p1: 1.5133   |    T2/T1: 1.1280   |   dV/a: 0.3056   |   pt2/pt1: 0.992798
    M: 1.30   |   M2: 0.7860   |    p2/p1: 1.8050   |    T2/T1: 1.1909   |   dV/a: 0.4423   |   pt2/pt1: 0.979374
    M: 1.40   |   M2: 0.7397   |    p2/p1: 2.1200   |    T2/T1: 1.2547   |   dV/a: 0.5714   |   pt2/pt1: 0.958194
    M: 1.50   |   M2: 0.7011   |    p2/p1: 2.4583   |    T2/T1: 1.3202   |   dV/a: 0.6944   |   pt2/pt1: 0.929787
    M: 1.60   |   M2: 0.6684   |    p2/p1: 2.8200   |    T2/T1: 1.3880   |   dV/a: 0.8125   |   pt2/pt1: 0.895200
    M: 1.70   |   M2: 0.6405   |    p2/p1: 3.2050   |    T2/T1: 1.4583   |   dV/a: 0.9265   |   pt2/pt1: 0.855721   
    M: 1.80   |   M2: 0.6165   |    p2/p1: 3.6133   |    T2/T1: 1.5316   |   dV/a: 1.0370   |   pt2/pt1: 0.812684
    M: 1.90   |   M2: 0.5956   |    p2/p1: 4.0450   |    T2/T1: 1.6079   |   dV/a: 1.1447   |   pt2/pt1: 0.767357
    M: 2.00   |   M2: 0.5774   |    p2/p1: 4.5000   |    T2/T1: 1.6875   |   dV/a: 1.2500   |   pt2/pt1: 0.720874
    >>> 
    """

    Mach_min = range[0]
    Mach_max = range[1]

    mach_nums = [i for i in np.arange(Mach_min,Mach_max+inc,inc)]
    M2 = [shock_mach(M1=i, gas=gas)for i in mach_nums]
    p2_p1 = [shock_pressure_ratio(M=i, gas=gas)for i in mach_nums]
    T2_T1 = [shock_temperature_ratio(M=i, gas=gas) for i in mach_nums]
    dv_a = [shock_dv_a(M=i, gas=gas) for i in mach_nums]
    pt2_pt1 = [shock_stagnation_ratio(M=i, gas=gas) for i in mach_nums]
    
    gamma, R = fluid(gas, metric)

    labl = '\u03B3 = ' + str(gamma)
    print("Normal Shock Parameters for " + labl)
    for index, num in enumerate(mach_nums):
        print("M: " + f"{num:.2f}" + "   |"+"   M2: " + f"{M2[index]:.4f}" + "   | " + "   p2/p1: " + f"{p2_p1[index]:.4f}" + "   | "+"   T2/T1: " + f"{T2_T1[index]:.4f}" + "   |"+"   dV/a: " + f"{dv_a[index]:.4f}"+ "   |"+"   pt2/pt1: " + f"{pt2_pt1[index]:.6f}" )
    print("\n \n \n")

#TODO: fix runtime double scalars warning
def shock_mach(M1=None, M2=None, gas='air', metric=True) -> float:
    """Returns the Mach number after a standing normal shock or the Mach number prior a standing normal shock
    
    Description
    -----------
    Given a starting Mach number M1 and the ratio of specific heats,
    return the Mach number M2 that immediately follows the shock. If M1 is
    not specified and M2 is, returns the Mach number prior to the shock.
    Default fluid is air.
    
    Parameters
    ----------
    M1 : `float`
        The Mach # before the shock\n
    M2 : `float`
        The Mach # after the shock\n
    gas : `str`
        The Fluid\n
    metric : `bool
        Use metric or US standard\n
    
    Examples
    --------
    >>> M2 = gd.shock_mach(M1=1.5) 
    >>> M2
    0.7010887416930995
    >>> M1 = gd.shock_mach(M2)
    >>> M1
    1.4999999999999998
    >>>
    """

    gamma, R = fluid(gas, metric)
    if M2 == None:
        #TODO: keep getting runtime double scalars warning what is this
        M2 = ((M1**2 + 2/(gamma-1)) / ((2*gamma / (gamma-1)) * M1**2 - 1))**.5
        return M2

    if M1 == None:
        M1 = ((-2/(gamma-1) -M2**2 ) / (1- ((2*gamma)/(gamma-1))*M2**2))**.5
        return M1

#GOOD #FIXME: maybe split these kinds of functions into two
def shock_pressure_ratio(M=None ,p2_p1=None, gas='air', metric=True) -> float:
    """Returns the pressure ratio after a standing normal shock for a given Mach number
    
    Description
    -----------
    Given a starting Mach number and a ratio of specific heats, this
    function returns the ratio of p2 / p1 across a standing normal shock.
    If Mach number is not specified and p2_p1 is, function returns Mach
    number. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    p2_p1 : `float`
        The pressure ratio\n
    gas : `str`
        The ratio of specific heats \n
    metric : `bool`
        Use metric or US standard \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> p2_p1 = gd.shock_pressure_ratio(M=1.5) 
    >>> p2_p1   
    2.4583333333333335
    >>>
    """

    gamma, R = fluid(gas, metric)

    if p2_p1 == None:
        p2_p1 = 2*gamma / (gamma+1) * M**2 - (gamma-1)/(gamma+1)
        return p2_p1

    if M == None:
        M = ((gamma+1)/(2*gamma) * (p2_p1 + (gamma-1)/(gamma+1)) )**.5
        return M

#GOOD
def shock_temperature_ratio(M: float, gas='air', metric=True) -> float:
    """Returns the temperature ratio after a standing normal shock for a given Mach number
    
    Description
    -----------
    Given a starting Mach number and a ratio of specific heats, this function
    returns the ratio of T2 / T1 across a standing normal shock. Default fluid
    is air.
    
    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> T2_T1 = gd.shock_temperature_ratio(M=1.5)
    >>> T2_T1
    1.320216049382716
    >>>
    """

    gamma, R = fluid(gas, metric)
    term1 = (1 + (gamma-1)/2 * M**2)
    term2 = (2*gamma)/(gamma-1) * M**2 -1
    term3 = (gamma+1)**2 / (2*(gamma-1)) * M**2
    t2_t1 = (term1 * term2) / term3
    return t2_t1

#TODO: add examples
def shock_dv_a(M: float, gas='air', metric=True) -> float:
    """Returns change in velocity over the local speed of sound after a normal shock.
    
    Description
    ----------
    Given a starting Mach # and a ratio of specific heats, this function
    returns the velocity change across a standing normal shock divided by the
    local speed of sound. Default fluid is air
    
    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    """

    gamma, R = fluid(gas, metric)
    dv_a = 2/(gamma+1) * (M**2 -1)/ M
    return dv_a

#TODO: add examples
def shock_stagnation_ratio(M: float, gas='air', metric=True) -> float:
    """Returns stagnation pressure ratio after a normal shock.
    
    Description
    -----------
    Given a starting Mach # and a ratio of specific heats, this function
    returns the ratio of pt2/pt1 across a standing normal shock. Default
    fluid is air.
    
    Parameters
    ----------
    M : `float`
        The starting Mach # \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    """

    gamma, R = fluid(gas, metric)
    term1 = (gamma+1)/2*M**2
    term2 = 1 + (gamma-1)/2 * M**2
    term3 = (term1/term2) ** (gamma/(gamma-1))
    term4 = (2*gamma / (gamma+1) * M**2 - ((gamma-1)/(gamma+1)))**(1/(1-gamma))
    return term3 * term4

#GOOD
def shock_flow_deflection(M: float, theta: float, gas='air', metric=True) -> float:
    """Returns flow deflection angle from Mach number and Oblique shock angle
    Description
    -----------
    Given the Mach number prior to the oblique shock, the angle of the oblique
    shock in degrees, and the ratio of specific heats, this function returns
    the angle that the flow is turned. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number before the shock \n
    theta : `float` 
        The shock angle in degrees \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> deflect = gd.shock_flow_deflection(M=2, theta = 22.5)
    >>> deflect
    -10.856560004139958
    >>>        
    """

    gamma, R = fluid(gas, metric)
    theta = theta * np.pi / 180 #degrees to radians
    dirac = np.arctan( 2 * 1/np.tan(theta) * (M**2 * np.sin(theta)**2 - 1 ) / (M**2 * (gamma + np.cos(2*theta)) + 2 ))
    dirac = dirac * 180 / np.pi
    return dirac

#GOOD
def shock_angle(M: float, dirac: float, gas='air', metric=True) -> float:
    """Return the shock angle given the Mach number prior to the shock and the deflection angle
    
    Description
    -----------
    Given the Mach number prior to the oblique shock, the angle of the flow
    deflection, and the ratio of specific heats, this functions returns the
    angle that is formed by the shock. Default ratio of specific heats is
    for air
    
    Parameters
    ----------
    M : `float`
        The Mach number before the shock \n
    dirac : `float`
        The flow deflection angle in degrees\n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> shocks = gd.shock_angle(M=2, dirac = -10) 
    >>> shocks
    [23.014012220565785, 96.29991962425305]
    >>> 
    """

    dirac = dirac * np.pi / 180
    gamma, R = fluid(gas, metric)
    def func(theta, M=M, dirac=dirac, gamma=gamma):
        zero = 2 * 1/np.tan(theta) * (M**2 * np.sin(theta)**2 - 1 ) / (M**2 * (gamma + np.cos(2*theta)) + 2 ) - np.tan(dirac)
        return zero

    weak = fsolve(func, x0=0.001, args=(M, dirac, gamma))
    strong = fsolve(func, x0=np.pi/2, args=(M, dirac, gamma))
    shock_angles = [weak[0] * 180/np.pi , strong[0] * 180/np.pi]
    return shock_angles

#TODO: edit docstring
def shock_mach_given_angles(theta: float, dirac: float, gas='air', metric=True) -> float:
    """Return the Mach number given the shock angle and flow deflection
    
    Description
    -----------
    #TODO: edit meee
    
    Parameters
    ----------
    theta : `float`
        The shock angle, in degrees \n
    dirac : `float`
        The flow deflection angle, in degrees \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = gd.shock_mach_given_angles(theta=22.5, dirac=10) 
    >>> M
    3.9293486839798955
    >>>
    """

    gamma, R = fluid(gas, metric)
    theta = theta * np.pi / 180
    dirac = dirac * np.pi / 180

    def func(M, theta=theta, dirac=dirac, gamma=gamma):
        '''
        Zero function for solving for the mach #
        '''
        zero = 2 * 1/np.tan(theta) * (M**2 * np.sin(theta)**2 - 1 ) / (M**2 * (gamma + np.cos(2*theta)) + 2 ) - np.tan(dirac)
        return zero

    sol = fsolve(func, x0=0.001, args=(theta, dirac, gamma))
    return sol[0]

#TODO: add examples
def prandtl_meyer_turn(M: float, gas='air', metric=True) -> float:
    """Returns the angle through which a flow has turned to reach a Mach number
    
    Description
    -----------
    Given a Mach number and ratio of specific heats, calculate angle through
    which a flow has turned to reach the Mach number given. Also known as
    the Prandtl-Meyer function. Default fluid is air
    
    Parameters
    ----------
    M : `float`
        The Mach number \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    #TODO:
    """

    gamma, R = fluid(gas, metric)
    nu = ((gamma+1)/(gamma-1))**.5 * np.arctan(((M**2-1)*(gamma-1)/(gamma+1))**.5) - np.arctan((M**2-1)**.5)
    return nu

#TODO: make this function
def prandtl_meyer_mach(gas='air', metric=True):
    """Returns the Mach number given an angle through which the flow has turned
    
    Description
    -----------
    Parameters:
    ----------
    Examples
    --------
    >>> 
    >>> 
    """
    gamma, R = fluid(gas, metric)
    
#GOOD for now
def shock_oblique_charts(Mach_max=6, gas='air', metric=True, lite=True):
    """Generate 2-D Oblique Shock Charts
    
    Description
    -----------
    Displays two plots,
    1) Mach number versus oblique shock wave angle and the corresponding
    deflection angles for each case. 2) Mach number versus resulting Mach
    number after an oblique shock and the corresponding deflection angles
    for each case. Default ratio of specific heats is for air.
    
    Parameters
    ----------
    Mach_max : `float`
        The upper limit Mach # for the chart \n
    gas : `str`
        The fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> gd.shock_oblique_charts()
    >>>
    """

    n = 1000
    gamma, R = fluid(gas, metric)

    theta = np.linspace(.001,90, n)
    mach = np.linspace(1,Mach_max,n)
    MACH,THETA = np.meshgrid(mach,theta)
    dirac = shock_flow_deflection(M=MACH, theta=THETA, gas=gas)
    plt.style.use('dark_background')
    fig, (ax1,ax2) = plt.subplots(1,2)
    levels=[0, 5,10,15,20,25,30,35,40]
    h = ax1.contour(mach,theta,dirac,levels=levels,cmap='tab10')
    ax1.clabel(h, inline =1, fontsize=10)
    minor_ticks_mach = np.arange(1,Mach_max+.1,.1)
    minor_ticks_theta = np.arange(0,91,1)
    ax1.set_yticks(minor_ticks_theta, minor=True)
    ax1.set_xticks(minor_ticks_mach, minor=True)
    ax1.grid(which='major',color='w',linestyle = '--', alpha=.5,)
    ax1.grid(which='minor',color='w',linestyle = '--', alpha=.1)
    ax1.set(xlabel = 'Mach #')
    ax1.set(ylabel = 'Oblique Shock Wave Angle')
    ax1_1 = ax1.twinx()
    ax1_1.set(ylabel = 'Flow Deflection Angle')
    ax1_1.get_xaxis().set_visible(False)
    ax1_1.get_yaxis().set_visible(False)

    if lite == True:
        n = 200
    else:
        n = 750

    total, counter = n**2, 0 
    mach_before = np.linspace(1,Mach_max, n)
    mach_after = np.linspace(.001,Mach_max, n)
    dirac = np.zeros((n,n))
    for row, m1 in enumerate(mach_before):
        percent =counter / total * 100
        print(' %0.1f%% complete' %percent)
        for col, m2 in enumerate(mach_after):
            dirac[col][row] = dirac_from_machs(M1=m1, M2=m2, gas=gas)
            counter += 1


    h = ax2.contour(mach_before, mach_after, dirac , levels=levels, cmap='tab10')
    ax2.clabel(h, inline = 1, fontsize=10)
    minor_ticks_mach_before = np.arange(0,Mach_max+.1,.1)
    minor_ticks_mach_after = np.arange(0,Mach_max,.1)
    ax2.set_yticks(minor_ticks_mach_before, minor=True)
    ax2.set_xticks(minor_ticks_mach_after, minor=True)
    ax2.grid(which='major',color='w',linestyle = '--', alpha=.5,)
    ax2.grid(which='minor',color='w',linestyle = '--', alpha=.1)
    ax2.set(xlabel = 'Mach # Before')
    ax2.set(ylabel = 'Mach # After')

    #plot title stuff
    gas = gas.lower()
    string = 'Oblique Shock Charts for ' + gas[0].upper() + gas[1:]
    fig.suptitle(string)
    fig.tight_layout(pad=2.0)
    plt.show()

#TODO: make dosctring
def dirac_from_machs(M1=[], M2=[], gas='air'):
    """Return the flow deflection angle and the shock angle required to go from one Mach # to a second Mach #
    Description
    -----------
    Parameters
    ----------
    Examples
    --------
    """   
    #get rid of the out of domain answers, if m2 >m1, throw it out.
    if M2 >= M1:
        return 0
    if shock_mach(M1) > M2:
        return 0

    def zero(theta, M1=[], M2=[], gas=gas):
        """Local function for testing different shock angles to see if they work.
        
        """
        theta_deg = degrees(theta)
        M1n = M1 * np.sin(theta)
        M2n = shock_mach(M1n)
        dirac = shock_flow_deflection(M=M1, theta=theta_deg, gas=gas, metric=True)
        M2_prime = M2n / np.sin(radians(theta_deg-dirac))
        zero = M2_prime - M2
        return zero

    thetas = np.linspace(0, np.pi/2, 90)
    for num, theta in enumerate(thetas[:-1]):
        zero1 = zero(thetas[num], M1=M1, M2=M2, gas=gas)
        zero2 = zero(thetas[num+1], M1=M1, M2=M2, gas=gas)
        if zero2 < 0:
            theta = lin_interpolate(0, zero1, zero2, thetas[num], thetas[num+1])
            return shock_flow_deflection(M=M1, theta=degrees(theta), gas=gas)


#TODO: make docstrings
def degrees(theta: float) -> float:
    """Convert from radians to degrees
    Parameters
    ----------
    theta : `float`
        The angle in radians \n
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> rad = gd.degrees(theta=3.14159)
    >>> rad
    179.99984796050427
    >>> 
    """

    return theta * 180/np.pi

#TODO: make docstrings
def radians(theta: float) -> float:
    """Convert from degrees to radians
    Parameters
    ----------
    theta : `float`
        The angle in degrees\n
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> deg = gd.radians(theta=180)
    >>> deg
    3.141592653589793
    >>> 
    """

    return theta * np.pi/180


def sind(theta: float) -> float:
    """Sine given degrees
    Parameters
    ----------
    theta : `float`
        The angle in degrees
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> sine = gd.sind(theta=45)
    >>> sine
    0.7071067811865476
    >>>
    """
    
    theta_radians = theta * np.pi / 180
    return np.sin(theta_radians)
    

def arcsind(sin: float) -> float:
    """Return the arcsine in degrees
    Parameters
    ----------
    sin : `float`
        The sine of theta
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> theta = gd.arcsind(sin = .7071)
    >>> theta
    44.99945053347443
    >>>
    """

    theta_rad = np.arcsin(sin)
    theta = degrees(theta_rad)
    return theta


def cosd(theta: float) -> float:
    """Cosine given degrees
    Parameters
    ----------
    theta : `float`
        The angle in degrees
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> cosine = gd.cosd(theta=45)
    >>> cosine
    0.7071067811865476
    >>>
    """

    theta_radians = theta * np.pi / 180
    return np.cos(theta_radians)


def arccosd(cos: float) -> float:
    """Return the arccosine in degrees
    Parameters
    ----------
    cos : `float`
        The cosine of theta
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> theta = gd.arccosd(.7071) 
    >>> theta
    45.00054946652557
    >>>
    """

    theta_rad = np.arccos(cos)
    theta = degrees(theta_rad)
    return theta


def tand(theta: float) -> float:
    """Tangent given degrees
    Parameters
    ----------
    theta : `float`
        The angle in degrees
    Examples
    --------
    >>> tan = gd.tand(theta=45)
    >>> tan
    0.9999999999999999
    >>>  
    """

    theta_rad = theta * np.pi / 180
    return np.tan(theta_rad)


def arctand(tan: float) -> float:
    """Return the arctangent in degrees
    Parameters
    ----------
    tan : `float`
        The tangent of theta
    Examples
    --------
    >>> theta = gd.arctand(1)
    >>> theta
    45.0
    >>>
    """
    
    theta_rad = np.arctan(tan)
    theta = degrees(theta_rad)
    return theta

#TODO: make docstrings
def area_dia(dia=[], area=[]):
    """
    Description
    -----------
    Given area or diameter, return the unknown
    Parameters
    ----------
    dia: Diameter \n
    area: Area \n
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> area = gd.area_dia(dia=1)    
    >>> area
    0.7853981633974483
    >>> dia = gd.area_dia(area=area)
    >>> dia 
    1.0
    >>>
    """

    if area == []:
        area = np.pi *dia**2 / 4
        return area
    if dia == []:
        dia = (area * 4 / np.pi )**.5
        return dia


#GOOD
def lin_interpolate(x, x0, x1, y0, y1):
    """Linear interpolation formula
    
    Description
    -----------
    Given two x values and their corresponding y values, interpolate for
    an unknown y value at x
    
    Parameters
    ----------
    x : `float
        The x value at which the y value is unknown \n
    x0 : `float`
        The x value of the first known pair \n
    x1 : `float
        The x value of the second known pair \n
    y0 : `float`
        The y value of the first known pair \n
    y1 : `float`
        The y value of the second known pair \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> y = gd.lin_interpolate(x=5, x0=4, x1=6, y0=25, y1=33)
    >>> y
    29.0
    >>>
    """

    y = y0 + (x-x0) * (y1-y0)/(x1-x0)
    return y