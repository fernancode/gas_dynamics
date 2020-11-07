#!usr/bin/env
"""Equatons, plots, and tables for working with shocks..
Included are functions to solve problems relating to properties across
normal and oblique shocks, from stagnation relations and changes in static
conditions to flow deflections. Tables can be made for any gas and its respective 
ratio of specific heats, as well as plots and charts of relationships.

  Typical usage example:
  Generate a shock tabgle for argon (two columns omitted for readability)
  >>> gd.shock_tables(range=[1,2], inc=.1, gas='argon') 
  Normal Shock Parameters for γ = 1.67
  M: 1.00   |   M2: 1.0000   |    p2/p1: 1.0000   |    T2/T1: 1.0000   |   
  M: 1.10   |   M2: 0.9131   |    p2/p1: 1.2627   |    T2/T1: 1.0985   |  
  M: 1.20   |   M2: 0.8463   |    p2/p1: 1.5504   |    T2/T1: 1.1956   |  
  M: 1.30   |   M2: 0.7935   |    p2/p1: 1.8631   |    T2/T1: 1.2933   |   
  M: 1.40   |   M2: 0.7509   |    p2/p1: 2.2009   |    T2/T1: 1.3934   |  
  M: 1.50   |   M2: 0.7158   |    p2/p1: 2.5637   |    T2/T1: 1.4968   |   
  M: 1.60   |   M2: 0.6866   |    p2/p1: 2.9515   |    T2/T1: 1.6042   |  
  M: 1.70   |   M2: 0.6620   |    p2/p1: 3.3643   |    T2/T1: 1.7162   |  
  M: 1.80   |   M2: 0.6410   |    p2/p1: 3.8021   |    T2/T1: 1.8331   |   
  M: 1.90   |   M2: 0.6229   |    p2/p1: 4.2649   |    T2/T1: 1.9552   |   
  M: 2.00   |   M2: 0.6073   |    p2/p1: 4.7528   |    T2/T1: 2.0827   |  

Copyright 2020 by Fernando A de la Fuente
All rights reserved
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from gas_dynamics.extra import ( fluid, air ,radians, degrees, sind, arcsind, cosd, arccosd, tand, arctand, lin_interpolate )



#==================================================
#shock_mach

#fix runtime double scalars warning
#==================================================
def shock_mach(M: float, gas='air', metric=True) -> float:
    """Returns the Mach number after a standing normal shock
    
    Description
    -----------
    Given a starting Mach number M1 and the ratio of specific heats,
    return the Mach number M2 that immediately follows the shock.
    Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number before the shock\n
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> M2 = gd.shock_mach(M1=1.5) 
    >>> M2
    0.7010887416930995
    >>> 
    """

    gamma, R = fluid(gas, metric)
    #TODO: keep getting runtime double scalars because of div/0
    M2 = ((M**2 + 2/(gamma-1)) / ((2*gamma / (gamma-1)) * M**2 - 1))**.5
    return M2



#==================================================
#shock_mach_before

#fix runtime double scalars warning
#==================================================
def shock_mach_before(M: float, gas='air', metric=True) -> float:
    """Returns the Mach number before a standing normal shock
    
    Description
    -----------
    Given a starting Mach number M1 and the ratio of specific heats,
    return the Mach number M2 that immediately follows the shock.
    Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number after the shock\n
    gas : `str`
        The Fluid\n
    metric : `bool`
        Use metric or US standard\n
    
    Examples
    --------
    >>> import gas_dynamics as gd 
    >>> M2 = 0.7010887416930995
    >>> M1 = gd.shock_mach(M2)
    >>> M1
    1.4999999999999998
    >>>
    """

    gamma, R = fluid(gas, metric)
    M1 = ((-2/(gamma-1) -M**2 ) / (1- ((2*gamma)/(gamma-1))*M**2))**.5
    return M1



#==================================================
#prandtl_meyer_mach
#need to add examples, parameters, descriptions
#maybe split these kinds of functions into two
#TODO: SPLIT INTO TWO FUNCS
#==================================================
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



#==================================================
#shock_temperature_ratio
#need to add examples, parameters, descriptions

#==================================================
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



#==================================================
# shock_dv_a
#need to add examples, parameters, descriptions, maybe change name?
#TODO: ADD EXAMPLES
#==================================================
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



#==================================================
#shock_stagnation ratio
# examples!
#TODO: ADD EXMAPLES
#==================================================
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



#==================================================
#shock_tables
# examples!

#==================================================
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
    M2 = [shock_mach(M=i, gas=gas)for i in mach_nums]
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



#==================================================
#shock_flow_deflection
# good!

#==================================================
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
    dirac = arctand( 2 * 1/tand(theta) * (M**2 * sind(theta)**2 - 1 ) / (M**2 * (gamma + cosd(2*theta)) + 2 ))    
    return dirac



#==================================================
#shock_angle
# good!

#==================================================
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

    gamma, R = fluid(gas, metric)
    def func(theta, M=M, dirac=dirac, gamma=gamma):
        zero = 2 * 1/tand(theta) * (M**2 * sind(theta)**2 - 1 ) / (M**2 * (gamma + cosd(2*theta)) + 2 ) - tand(dirac)
        return zero

    weak = fsolve(func, x0=0.001, args=(M, dirac, gamma))
    strong = fsolve(func, x0=90, args=(M, dirac, gamma))
    shock_angles = [weak[0], strong[0]]
    return shock_angles



#==================================================
#shock_mach_given angles
#need to add descriptions

#==================================================
def shock_mach_given_angles(theta: float, dirac: float, gas='air', metric=True) -> float:
    """Return the Mach number given the shock angle and flow deflection
    
    Description
    -----------
    
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

    def func(M, theta=theta, dirac=dirac, gamma=gamma):
        '''
        Zero function for solving for the mach #
        '''
        zero = 2 * 1/tand(theta) * (M**2 * sind(theta)**2 - 1 ) / (M**2 * (gamma + cosd(2*theta)) + 2 ) - tand(dirac)
        return zero

    sol = fsolve(func, x0=0.001, args=(theta, dirac, gamma))
    return sol[0]



#==================================================
#prandtl_meyer_turn
#need to add examples

#==================================================
def prandtl_meyer_turn(M: float, gas='air', metric=True) -> float:
    """Returns the angle through which a flow has turned to reach a Mach number
    
    Description
    -----------
    Given a Mach number and ratio of specific heats, calculate angle through
    which a flow has turned to reach the Mach number given from a starting
    Mach number of 1. Also known as the Prandtl-Meyer function. Default fluid
    is air
    
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

    """

    gamma, R = fluid(gas, metric)
    nu = ((gamma+1)/(gamma-1))**.5 * arctand(((M**2-1)*(gamma-1)/(gamma+1))**.5) - arctand((M**2-1)**.5)
    return nu



#==================================================
#prandtl_meyer_mach
#need to add examples, parameters, descriptions

#==================================================
def prandtl_meyer_mach(nu: float, gas='air', metric=True) -> float:
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
    
    def get_mach(M: float, nu=nu, gas=gas) -> float:
        return prandtl_meyer_turn(M, gas=gas) - nu
    
    sol = fsolve(get_mach, x0=1.5, args=(nu, gas))
    return sol[0]



#==================================================
#shock_oblique_charts
#need to add examples, parameters, descriptions

#==================================================
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



#==================================================
#dirac from machs
#need to add examples, parameters, descriptions

#==================================================
def dirac_from_machs(M1=[], M2=[], gas='air') -> float:
    """Return the flow deflection angle and the shock angle
    required to go from one Mach # to a second Mach #
    
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
        M1n = M1 * sind(theta)
        M2n = shock_mach(M1n)
        dirac = shock_flow_deflection(M=M1, theta=theta, gas=gas, metric=True)
        M2_prime = M2n / sind(theta-dirac)
        zero = M2_prime - M2
        return zero

    thetas = np.linspace(0, 90, 90)
    for num, theta in enumerate(thetas[:-1]):
        zero1 = zero(thetas[num], M1=M1, M2=M2, gas=gas)
        zero2 = zero(thetas[num+1], M1=M1, M2=M2, gas=gas)
        if zero2 < 0:
            theta = lin_interpolate(0, zero1, zero2, thetas[num], thetas[num+1])
            return shock_flow_deflection(M=M1, theta=theta, gas=gas)