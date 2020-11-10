#!usr/bin/env
#Equations, plots, and tables for working with shocks.
#Included are functions for solving problems related to properties across
#normal and oblique shocks, from stagnation relations and changes to
#conditions to flow deflections. Tables can be made for any gas and its respective 
#ratio of specific heats, as well as plots and charts of relationships.
#
#  Typical usage example:
#  Generate a shock tabgle for argon (two columns omitted for readability)
#  >>> gd.shock_tables(range=[1,2], step=.1, gas=argon) 
#  Normal Shock Parameters for γ = 1.67
#  M: 1.00   |   M2: 1.0000   |    p2/p1: 1.0000   |    T2/T1: 1.0000   |   
#  M: 1.10   |   M2: 0.9131   |    p2/p1: 1.2627   |    T2/T1: 1.0985   |  
#  M: 1.20   |   M2: 0.8463   |    p2/p1: 1.5504   |    T2/T1: 1.1956   |  
#  M: 1.30   |   M2: 0.7935   |    p2/p1: 1.8631   |    T2/T1: 1.2933   |   
#  M: 1.40   |   M2: 0.7509   |    p2/p1: 2.2009   |    T2/T1: 1.3934   |  
#  M: 1.50   |   M2: 0.7158   |    p2/p1: 2.5637   |    T2/T1: 1.4968   |   
#  M: 1.60   |   M2: 0.6866   |    p2/p1: 2.9515   |    T2/T1: 1.6042   |  
#  M: 1.70   |   M2: 0.6620   |    p2/p1: 3.3643   |    T2/T1: 1.7162   |  
#  M: 1.80   |   M2: 0.6410   |    p2/p1: 3.8021   |    T2/T1: 1.8331   |   
#  M: 1.90   |   M2: 0.6229   |    p2/p1: 4.2649   |    T2/T1: 1.9552   |   
#  M: 2.00   |   M2: 0.6073   |    p2/p1: 4.7528   |    T2/T1: 2.0827   |  
#
#Copyright 2020 by Fernando A de la Fuente
#All rights reserved

#TODO: everything needs a 
#returns
#   `type`


import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from gas_dynamics.extra import ( radians, degrees, sind, arcsind, cosd, arccosd, tand, arctand, lin_interpolate )
from gas_dynamics.fluids import fluid, air



#==================================================
#shock_mach
#fix runtime double scalars warning
#TODO: maybe add a "returns" for all functions
#FIXME: the desriptions sound wordy but i think it makes sense to be like that:/
#==================================================
def shock_mach(M: float, gas=air) -> float:
    """Returns the Mach number after a standing normal shock


    Notes
    -----
    Given a starting Mach number M1 and the ratio of specific heats,
    return the Mach number M2 that immediately follows the shock.
    The default fluid is air.


    Parameters
    ----------
    M : `float`
        The Mach number before the shock\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The Mach number after the shock\n

    Examples
    --------
    >>> M2 = gd.shock_mach(M1=1.5) 
    >>> M2
    0.7010887416930995
    >>> 
    """
    
    gamma = gas.gamma

    #TODO: still getting runtime double scalars, filter out for this better
    M2 = ((M**2 + 2/(gamma-1)) / ((2*gamma / (gamma-1)) * M**2 - 1))**.5
    return M2



#==================================================
#shock_mach_before
#==================================================
def shock_mach_before(M: float, gas=air) -> float:
    """Returns the Mach number before a standing normal shock
    
    Notes
    -----
    Given the Mach number after the shock and the ratio of specific heats,
    return the Mach number that immediately precedes the shock.
    Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number after the shock\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The mach number before the shock\n
    
    Examples
    --------
    >>> import gas_dynamics as gd 
    >>> M2 = 0.7010887416930995
    >>> M1 = gd.shock_mach(M2)
    >>> M1
    1.4999999999999998
    >>>
    """

    gamma = gas.gamma
    M1 = ((-2/(gamma-1) -M**2 ) / (1- ((2*gamma)/(gamma-1))*M**2))**.5
    return M1



#==================================================
#prandtl_meyer_mach
#added fluid class, split into two funcs
#==================================================
def shock_pressure_ratio(M: float, gas=air) -> float:
    """Returns the pressure ratio after a standing normal shock for a given Mach number
    
    Notes
    -----
    Given a starting Mach number and a ratio of specific heats, this
    function returns the ratio of pressure 2 over pressure 1 across a
    standing normal shock. 

    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The pressure ratio p2/p1\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> p2_p1 = gd.shock_pressure_ratio(M=1.5) 
    >>> p2_p1   
    2.4583333333333335
    >>>
    """

    gamma = gas.gamma
    p2_p1 = 2*gamma / (gamma+1) * M**2 - (gamma-1)/(gamma+1)
    return p2_p1



#==================================================
#shock mach from pressure ratio
#==================================================
def shock_mach_from_pressure_ratio(p2_p1: float, gas=air) -> float:
    """Returns the mach number across a standing normal shock given a pressure ratio
    
    Notes
    -----
    Given the ratio of pressure behind the shock over pressure before
    the shock and the ratio of specific heats, this function returns the 
    Mach number before the shock. Default fluid is air.
    
    Parameters
    ----------
    p2_p1 : `float`
        The pressure ratio\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The mach number prior the shock\n
    
    Examples
    --------
    >>> import gas_dynamics as gd 
    >>> p2_p1 = 6
    >>> M1 = gd.shock_mach_from_pressure_ratio(p2_p1)
    >>> M1
    2.29906813420444
    >>> 
    """

    gamma = gas.gamma
    M = ((gamma+1)/(2*gamma) * (p2_p1 + (gamma-1)/(gamma+1)) )**.5
    return M



#==================================================
#shock_temperature_ratio
#==================================================
def shock_temperature_ratio(M: float, gas=air) -> float:
    """Returns the temperature ratio after a standing normal shock for a given Mach number
    
    Notes
    -----
    Given a starting Mach number and a ratio of specific heats, this function
    returns the ratio of temperature 2 over temperature 1 across a standing normal
    shock. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The temperature ratio T2/T1\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> T2_T1 = gd.shock_temperature_ratio(M=1.5)
    >>> T2_T1
    1.320216049382716
    >>>
    """

    gamma = gas.gamma
    term1 = (1 + (gamma-1)/2 * M**2)
    term2 = (2*gamma)/(gamma-1) * M**2 -1
    term3 = (gamma+1)**2 / (2*(gamma-1)) * M**2
    t2_t1 = (term1 * term2) / term3
    return t2_t1



#==================================================
# shock_dv_a
#need to add examples, maybe change name?
#TODO: ADD EXAMPLES
#==================================================
def shock_dv_a(M: float, gas=air) -> float:
    """Returns change in velocity over the local speed of sound after a normal shock.
    
    Notes
    -----
    Given a starting Mach number and a ratio of specific heats, this function
    returns the velocity change across a standing normal shock divided by the
    local speed of sound. Default fluid is air
    
    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The change in Mach number\n

    Examples
    --------

    """

    gamma = gas.gamma
    dv_a = 2/(gamma+1) * (M**2 -1)/ M
    return dv_a



#==================================================
#shock_stagnation ratio
#==================================================
def shock_stagnation_ratio(M: float, gas=air) -> float:
    """Returns stagnation pressure ratio after a normal shock.
    
    Notes
    -----
    Given a starting Mach number and a ratio of specific heats, this function
    returns the ratio of stagnation presure 2 over stagnation pressure 1
    across a standing normal shock. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The starting Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n    

    Returns
    -------
    float
        The stagnation pressure ratio pt2/pt1\n

    Examples
    --------
    >>> pt2_pt1 = gd.shock_stagnation_ratio(M=2)
    >>> pt2_pt1
    0.7208738614847454
    >>>
    """

    gamma = gas.gamma
    term1 = (gamma+1)/2*M**2
    term2 = 1 + (gamma-1)/2 * M**2
    term3 = (term1/term2) ** (gamma/(gamma-1))
    term4 = (2*gamma / (gamma+1) * M**2 - ((gamma-1)/(gamma+1)))**(1/(1-gamma))
    return term3 * term4



#==================================================
#shock_tables
# examples!

#==================================================
def shock_tables(range=[1,5], step=.01, gas=air) -> str:
    """Returns shock tables for a range of Mach numbers.
    
    Notes
    -----
    Given a range of Mach numbers and a ratio of specific heats, generate
    the standing normal shock tables for every incremental Mach number
    in between.
    
    Parameters
    ----------
    range : `list`
        The starting and ending Mach # in a list, ie: [1,5]. \n
    step : `float`
        The step size for the tables. \n
    gas : `fluid`
        A user defined fluid object. Default is air \n    
    
    Returns
    -------
    str
        The shock table\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> gd.shock_tables(range=[1,2], step=.1)
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
    if Mach_min < 1:
        Mach_min = 1

    Mach_max = range[1]

    mach_nums = [i for i in np.arange(Mach_min,Mach_max+step,step)]
    M2 = [shock_mach(M=i, gas=gas)for i in mach_nums]
    p2_p1 = [shock_pressure_ratio(M=i, gas=gas)for i in mach_nums]
    T2_T1 = [shock_temperature_ratio(M=i, gas=gas) for i in mach_nums]
    dv_a = [shock_dv_a(M=i, gas=gas) for i in mach_nums]
    pt2_pt1 = [shock_stagnation_ratio(M=i, gas=gas) for i in mach_nums]
    gamma = gas.gamma

    labl = '\u03B3 = ' + str(gamma)
    print("Normal Shock Parameters for " + gas.name + ", " + labl)
    for index, num in enumerate(mach_nums):
        print("M: " + f"{num:.2f}" + "   |"+"   M2: " + f"{M2[index]:.4f}" + "   | " + "   p2/p1: " + f"{p2_p1[index]:.4f}" + "   | "+"   T2/T1: " + f"{T2_T1[index]:.4f}" + "   |"+"   dV/a: " + f"{dv_a[index]:.4f}"+ "   |"+"   pt2/pt1: " + f"{pt2_pt1[index]:.6f}" )
    print("\n \n \n")



#==================================================
#shock_flow_deflection
# good! added fluid class
#==================================================
def shock_flow_deflection(M: float, theta: float, gas=air) -> float:
    """Returns flow deflection angle from Mach number and oblique shock angle
    
    Notes
    -----
    Given the Mach number prior to the oblique shock, the angle of the oblique
    shock in degrees, and the ratio of specific heats, this function returns
    the angle that the flow is turned. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number before the shock \n
    theta : `float`
        The shock angle in degrees \n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The deflection angle\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> deflect = gd.shock_flow_deflection(M=2, theta = 22.5)
    >>> deflect
    -10.856560004139958
    >>>        
    """

    gamma = gas.gamma
    dirac = arctand( 2 * 1/tand(theta) * (M**2 * sind(theta)**2 - 1 ) / (M**2 * (gamma + cosd(2*theta)) + 2 ))    
    return dirac



#==================================================
#shock_angle
# good! added fluid class
#==================================================
def shock_angle(M: float, dirac: float, gas=air) -> float:
    """Return the shock angle given the Mach number prior to the shock and the deflection angle
    
    Notes
    -----
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
    gas : `fluid`
        A user defined fluid object. Default is air \n    
    
    Returns
    -------
    float
        The shock angle\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> shocks = gd.shock_angle(M=2, dirac = -10) 
    >>> shocks
    [23.014012220565785, 96.29991962425305]
    >>> 
    """

    gamma = gas.gamma
    def func(theta, M=M, dirac=dirac, gamma=gamma):
        zero = 2 * 1/tand(theta) * (M**2 * sind(theta)**2 - 1 ) / (M**2 * (gamma + cosd(2*theta)) + 2 ) - tand(dirac)
        return zero

    weak = fsolve(func, x0=0.001, args=(M, dirac, gamma))
    strong = fsolve(func, x0=90, args=(M, dirac, gamma))
    shock_angles = [weak[0], strong[0]]
    return shock_angles



#==================================================
#shock_mach_given angles
#==================================================
def shock_mach_given_angles(theta: float, dirac: float, gas=air) -> float:
    """Return the Mach number given the shock angle and flow deflection
    
    Notes
    -----
    Given the angle of the shock and the angle that the flow has turned,
    return the mach number that preceded the shock.

    Parameters
    ----------
    theta : `float`
        The shock angle in degrees \n
    dirac : `float`
        The flow deflection angle in degrees \n
    gas : `fluid`
        A user defined fluid object. Default is air \n    
    
    Returns
    -------
    float
        The mach number\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = gd.shock_mach_given_angles(theta=22.5, dirac=10) 
    >>> M
    3.9293486839798955
    >>>
    """

    gamma = gas.gamma
    def func(M, theta=theta, dirac=dirac, gamma=gamma):
        '''
        Zero function for solving for the mach number
        '''
        zero = 2 * 1/tand(theta) * (M**2 * sind(theta)**2 - 1 ) / (M**2 * (gamma + cosd(2*theta)) + 2 ) - tand(dirac)
        return zero

    sol = fsolve(func, x0=0.001, args=(theta, dirac, gamma))
    return sol[0]



#==================================================
#prandtl_meyer_turn
#need to add examples
#==================================================
def prandtl_meyer_turn(M: float, gas=air) -> float:
    """Returns the angle through which a flow has turned to reach a Mach number
    
    Notes
    -----
    Given a Mach number and ratio of specific heats, calculate the angle of a turn
    through which a flow has traversed to reach the Mach number given, from a  Mach number
    of 1. Also known as the Prandtl-Meyer function. Default fluid is air.
    
    Parameters
    ----------
    M : `float`
        The Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n    
    
    Returns
    -------
    float
        The angle through which the flow has turned\n

    Examples
    --------

    """

    gamma = gas.R
    nu = ((gamma+1)/(gamma-1))**.5 * arctand(((M**2-1)*(gamma-1)/(gamma+1))**.5) - arctand((M**2-1)**.5)
    return nu



#==================================================
#prandtl_meyer_mach
#need to add examples
#==================================================
def prandtl_meyer_mach(nu: float, gas=air) -> float:
    """Returns the Mach number given an angle through which the flow has turned
    
    Notes
    -----
    Given a smooth turn through which a flow has turned and the ratio of specific
    heats, return the Mach number after the turn.

    Parameters
    ----------
    nu : `float`
        The turn angle \n    
    gas : `fluid`
        A user defined fluid object. Default is air \n    

    Returns
    -------
    float
        The mach number\n

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
def shock_oblique_charts(Mach_max=6, gas=air, lite=True, dark=True):
    """Generate 2-D Oblique Shock Charts
    
    Notes
    -----
    Displays two plots,
    1) Mach number versus oblique shock wave angle and the corresponding
    deflection angles for each case. 2) Mach number versus resulting Mach
    number after an oblique shock and the corresponding deflection angles
    for each case. Default ratio of specific heats is for air.
    
    Parameters
    ----------
    Mach_max : `float`
        The upper limit Mach number for the chart \n
    gas : `fluid`
        A user defined fluid object. Default is air \n    
    lite : `bool`
        Calculate 40,000 pts in the mesh vs 562,500
    dark : `bool`
        Dark mode for the plots. Default is true.
    
    Examples
    --------
    >>> gd.shock_oblique_charts()
    >>>
    """

    n = 1000
    gamma = gas.gamma

    if dark == True:
        plt.style.use('dark_background')
        grid_color = 'w'
    else: 
        plt.style.use('default')
        grid_color = 'k'


    theta = np.linspace(.001,90, n)
    mach = np.linspace(1,Mach_max,n)
    MACH,THETA = np.meshgrid(mach,theta)
    dirac = shock_flow_deflection(M=MACH, theta=THETA, gas=gas)
    fig, (ax1,ax2) = plt.subplots(1,2)
    levels=[0, 5,10,15,20,25,30,35,40]
    h = ax1.contour(mach,theta,dirac,levels=levels,cmap='tab10')
    ax1.clabel(h, inline =1, fontsize=10)
    minor_ticks_mach = np.arange(1,Mach_max+.1,.1)
    minor_ticks_theta = np.arange(0,91,1)
    ax1.set_yticks(minor_ticks_theta, minor=True)
    ax1.set_xticks(minor_ticks_mach, minor=True)
    ax1.grid(which='major', color=grid_color, linestyle = '--', alpha=.5,)
    ax1.grid(which='minor', color=grid_color, linestyle = '--', alpha=.1)
    ax1.set(xlabel = 'Mach #')
    ax1.set(ylabel = 'Oblique Shock Wave Angle')

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
    ax2.grid(which='major',color=grid_color,linestyle = '--', alpha=.5,)
    ax2.grid(which='minor',color=grid_color,linestyle = '--', alpha=.1)
    ax2.set(xlabel = 'Mach # Before')
    ax2.set(ylabel = 'Mach # After')

    #plot title stuff
    string = 'Oblique Shock and Flow Deflection Angles for ' + gas.name
    fig.suptitle(string)
    fig.tight_layout(pad=2.0)
    plt.show()



#==================================================
#dirac from machs
##TODO: better numerical method for this
##TODO: examples
#==================================================
def dirac_from_machs(M1: float, M2: float, gas=air) -> float:
    """Return the flow deflection angle given the mach number before
    and after the oblique shock
    
    Notes
    -----
    Given two mach numbers, iteratively solve for the shock angle and
    flow deflection to satisfy the system.

    Parameters
    ----------
    M1 : `float`
        The initial mach number \n
    M2 : `float`
        The mach number after the event \n
    gas : `fluid`
        A user defined fluid object. Default is air \n 

    Returns
    -------
    float
        The flow deflection angle\n

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
        dirac = shock_flow_deflection(M=M1, theta=theta, gas=gas)
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