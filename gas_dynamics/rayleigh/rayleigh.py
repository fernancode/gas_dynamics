#!usr/bin/env
#Equations, plots, and tables for working with constant area flow
#where heat transfer plays a major role and frictional effects are
#neglected, otherwise known as Rayleigh flow
#
#
#Copyright 2020 by Fernando A de la Fuente
#All rights reserved
from gas_dynamics.fluids import fluid, air
from scipy.optimize import fsolve


#==================================================
#rayleigh pressure ratio
#==================================================
def rayleigh_pressure_ratio(mach_initial: float, mach_final:float, gas=air) -> float:
    """Return the pressure ratio pressure_final / pressure_initial given the two Mach numbers

    Notes
    -----
    Given two Mach numbers, determine the resulting pressure ratio pressure two
    over pressure one in the non-adiabatic constant area frictionless flow. 
    Default fluid is air.

    Parameters
    ----------
    mach_initial : `flaot`
        The Mach number at region 1\n
    mach_final : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The rayleigh pressure ratio pressure_final / pressure_initial \n
    
    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, mach_final = 3, 2
    >>> p2_p1 = gd.rayleigh_pressure_ratio(mach_initial, mach_final)
    >>> p2_p1
    2.0606060606060606
    >>>
    """
    gamma = gas.gamma
    p2_p1 = (1 + gamma * mach_initial**2)/(1 + gamma *  mach_final**2)
    return p2_p1



#==================================================
#rayleigh temperature ratio
#==================================================
def rayleigh_temperature_ratio(mach_initial: float, mach_final: float, gas=air) -> float:
    """Return the temperature ratio T2/T1 given the two Mach numbers

    Notes
    -----
    Given two Mach numbers, determine the resulting Temperature ratio
    temperature two over temperature one in the non-adiabatic constant 
    area frictionless flow. Default fluid is air.    

    Parameters
    ----------
    mach_initial : `flaot`
        The Mach number at region 1\n
    mach_final : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The rayleigh temperature ratio T2 / T1 \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, mach_final = .8, .3
    >>> T2_T1 = gd.rayleigh_temperature_ratio(mach_initial, mach_final)
    >>> T2_T1
    0.39871485855083627
    >>>
    """

    gamma = gas.gamma
    T2_T1 = ((1 + gamma * mach_initial**2)/(1 + gamma *  mach_final**2))**2 * (mach_final**2)/(mach_initial**2)
    return T2_T1



#==================================================
#rayleigh density ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_density_ratio(mach_initial: float, mach_final: float, gas=air) -> float:
    """Return the density ratio rho2/rho1 given the two Mach numbers

    Notes
    -----
    Given two Mach numbers, determine the resulting density ratio
    density two over density one in the non-adiabatic constant 
    area frictionless flow. Default fluid is air.    

    Parameters
    ----------
    mach_initial : `flaot`
        The Mach number at region 1\n
    mach_final : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The rayleigh density ratio rho2 / rho1 \n

    Examples
    --------

    """

    gamma = gas.gamma
    rho2_rho1 = (1 + gamma * mach_final**2)/(1 + gamma *  mach_initial**2) * (mach_initial**2)/(mach_final**2)
    return rho2_rho1



#==================================================
#rayleigh stagnation temperature ratio
#==================================================
def rayleigh_stagnation_temperature_ratio(mach_initial: float, mach_final: float, gas=air) -> float:
    """Return the stagnation temperature ratio Tt2/Tt1 given the two Mach numbers

    Notes
    -----
    Given two Mach numbers, determine the resulting stagnation temperature ratio
    stagnation temperature two over stagnation temperature one in the non-adiabatic 
    constant area frictionless flow. Default fluid is air.    

    Parameters
    ----------
    mach_initial : `flaot`
        The Mach number at region 1\n
    mach_final : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The rayleigh stagnation temperature ratio pt2 / pt1 \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, mach_final = .8, .3
    >>> Tt2_Tt1 = gd.rayleigh_stagnation_temperature_ratio(mach_initial, mach_final)
    >>> Tt2_Tt1
    0.35983309042974404
    >>>
    """
    gamma = gas.gamma
    Tt2_Tt1 = ((1+gamma*mach_initial**2)/(1+gamma*mach_final**2))**2 * (mach_final**2)/(mach_initial**2) * ((1+(gamma-1)/2*mach_final**2)/(1+(gamma-1)/2*mach_initial**2))
    return Tt2_Tt1



#==================================================
#rayleigh stagnation pressure ratio
#==================================================
def rayleigh_stagnation_pressure_ratio(mach_initial: float, mach_final: float, gas=air) -> float:
    """Return the stagnation pressure ratio pt2/pt1 given the two Mach numbers

    Notes
    -----
    Given two Mach numbers, determine the resulting stagnation pressure ratio
    stagnation pressure two over stagnation pressure one in the non-adiabatic 
    constant area frictionless flow. Default fluid is air.    

    Parameters
    ----------
    mach_initial : `flaot`
        The Mach number at region 1\n
    mach_final : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The rayleigh stagnation pressure ratio pt2 / pt1 \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, mach_final = .8, .3
    >>> pt2_pt1 = gd.rayleigh_stagnation_pressure_ratio(mach_initial, mach_final)
    >>> pt2_pt1
    1.1758050380938454
    >>>
    """

    gamma = gas.gamma
    pt2_pt1 = (1+gamma*mach_initial**2)/(1+gamma*mach_final**2) * ((1+(gamma-1)/2*mach_final**2)/(1+(gamma-1)/2*mach_initial**2))**(gamma/(gamma-1))
    return pt2_pt1



#==================================================
#rayleigh mach from pressure ratio
#==================================================
def rayleigh_mach_from_pressure_ratio(mach_initial: float, pressure_initial: float, pressure_final: float, gas=air) -> float:
    """Return the mach number given the Mach number and two pressures

    Notes
    -----
    Given the initial Mach number, initial pressure, and final pressure, determine
    the resulting Mach number in the non-adiabatic constant area frictionless flow with 
    heat transfer. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The Mach number\n
    pressure_initial : `flaot`
        pressure 1\n
    pressure_final : `float`
        pressure 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The resulting Mach number \n

    Examples
    --------
    >>> import gas_dynamics as gd                                        
    >>> mach_final = gd.rayleigh_mach_from_pressure_ratio(mach_initial=.8, pressure_initial=1.5, pressure_final=2.5) 
    >>> mach_final
    0.31350552512789054
    >>>
    """

    gamma = gas.gamma
    mach_final = (((pressure_initial*(1+gamma*mach_initial**2))/pressure_final -1)/gamma)**.5
    return mach_final



#==================================================
#rayleigh mach from temperature ratio
#==================================================
def rayleigh_mach_from_temperature_ratio(mach_initial: float, temperature_initial: float, temperature_final: float, gas=air) -> float:
    """Return the Mach number given the Mach number and two temperatures

    Notes
    -----
    Given the initial Mach number, initial temperature, and final temperature, determine
    the resulting Mach number in the non-adiabatic constant area frictionless flow with 
    heat transfer. Default fluid is air.

    Parameters
    ----------
    mach : `float`
        The Mach number\n
    temperature_initial : `float`
        Temperature 1\n
    temperature_final : `float`
        Temperature 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The resulting Mach number \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, T1, T2 = .8, 250, 100
    >>> mach_final = gd.rayleigh_mach_from_temperature_ratio(mach_initial, T1, T2)
    >>> mach_final
    0.3006228581671002
    >>>
    """

    gamma = gas.gamma
    T2_T1 = temperature_final/temperature_initial
    def zero(mach_final, mach_initial, T2_T1,gas):
        return rayleigh_temperature_ratio(mach_initial=mach_initial, mach_final=mach_final, gas=gas) - T2_T1
    
    if mach_initial < 1:
        x0 = .5
    elif mach_initial > 1:
        x0 = 1.5
    else:
        x0 = 1

    sol = fsolve(zero, args=(mach_initial, T2_T1, gas), x0=x0)
    return sol[0]



#==================================================
#rayleigh mach from stagnation temperature ratio
#==================================================
def rayleigh_mach_from_stagnation_temperature_ratio(mach_initial: float, stagnation_temperature_initial: float, stagnation_temperature_final: float, gas=air) -> float:
    """Return the Mach number given the Mach number and two stagnation temperatures

    Notes
    -----
    Given the initial Mach number, initial stagnation temperature, and final stagnation
    temperature, determine the resulting Mach number in the non-adiabatic constant area
    frictionless flow with heat transfer. Default fluid is air.

    Parameters
    ----------
    mach_initial : `float`
        The Mach number\n
    stagnation_temperature_initial : `float`
        Stagnation temperature 1\n
    stagnation_temperature_final : `float`
        Stagnation temperature 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The resulting Mach number \n


    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, Tt1, Tt2 = .8, 250, 105
    >>> mach_final = gd.rayleigh_mach_from_stagnation_temperature_ratio(mach_initial, Tt1, Tt2)
    >>> mach_final
    0.33147520792270446
    >>>
    """

    gamma = gas.gamma
    Tt2_Tt1 = stagnation_temperature_final/stagnation_temperature_initial
    def zero(mach_final, mach_initial, Tt2_Tt1, gas):
        return rayleigh_stagnation_temperature_ratio(mach_initial=mach_initial, mach_final=mach_final, gas=gas) - Tt2_Tt1
    
    if mach_initial < 1:
        x0 = .5
    elif mach_initial > 1:
        x0 = 1.5
    else:
        x0 = 1

    sol = fsolve(zero, args=(mach_initial, Tt2_Tt1, gas), x0=x0)
    return sol[0]



#==================================================
#rayleigh mach from stagnation pressure ratio
#==================================================
def rayleigh_mach_from_stagnation_pressure_ratio(mach_initial: float, stagnation_pressure_initial: float, stagnation_pressure_final: float, gas=air) -> float:
    """Return the Mach number given the Mach number and two stagnation pressures

    Notes
    -----
    Given the initial Mach number, initial stagnation pressure, and final stagnation
    pressure, determine the resulting Mach number in the non-adiabatic constant area
    frictionless flow with heat transfer. Default fluid is air.

    Parameters
    ----------
    M : `float`
        The Mach number\n
    stagnation_pressure_initial : `float`
        Stagnation pressure 1\n
    stagnation_pressure_final : `float`
        Stagnation pressure 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The resulting Mach number\n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> mach_initial, pt1, pt2 = .8, 2.3, 2.6
    >>> mach_final = gd.rayleigh_mach_from_stagnation_pressure_ratio(mach_initial, pt1, pt2)
    >>> mach_final
    0.40992385119887326
    >>>
    """

    gamma = gas.gamma
    pt2_pt1 = stagnation_pressure_final/stagnation_pressure_initial
    def zero(mach_final, mach_initial, pt2_pt1, gas):
        return rayleigh_stagnation_pressure_ratio(mach_initial=mach_initial, mach_final=mach_final, gas=gas) - pt2_pt1

    if mach_initial < 1:
        x0 = .5
    elif mach_initial > 1:
        x0 = 1.5
    else:
        x0 = 1

    sol = fsolve(zero, args=(mach_initial, pt2_pt1, gas), x0=x0)
    return sol[0]



#==================================================
#rayleigh pressure star ratio
#==================================================
def rayleigh_pressure_star_ratio(mach: float, gas=air) -> float:
    """Return the ratio of pressure over pressure where Mach is equal to one

    Notes
    -----
    Given a mach number, return the ratio of the pressure of the mach number
    over the pressure over where the mach number is equal to one for a non-adiabatic
    fricionless constant area system. Default fluid is air.

    Parameters
    ----------
    M : `flaot`
        The Mach number \n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The ratio of p / p* \n

    Examples
    --------
    >>> gas_dynamics as gd
    >>> M = 2
    >>> p_pstar = gd.rayleigh_pressure_star_ratio(M)
    >>> p_pstar
    0.36363636363636365
    >>>
    """

    gamma = gas.gamma
    p_pstar = (gamma+1)/(1+gamma*mach**2)
    return p_pstar



#==================================================
#rayleigh temperature star ratio
#==================================================
def rayleigh_temperature_star_ratio(mach: float, gas=air) -> float:
    """Return the ratio of temperature over temperature where Mach is equal to one

    Notes
    -----
    Given a mach number, return the ratio of the temperature of the mach number
    over the temperature over where mach number is equal to one for a non-adiabatic
    frictionless constant area system. Default fluid is air.

    Parameters
    ----------
    mach : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The ratio of T / T* \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 2
    >>> T_Tstar = gd.rayleigh_temperature_star_ratio(M)
    >>> T_Tstar
    0.5289256198347108
    >>>
    """

    gamma = gas.gamma
    T_Tstar = (mach**2 * (1+gamma)**2) / (1+gamma*mach**2)**2
    return T_Tstar



#==================================================
#rayleigh density star ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_density_star_ratio(mach: float, gas=air) -> float:
    """Return the ratio of density over density where Mach is equal to one

    Notes
    -----
    Given a mach number, return the ratio of the density of the mach number
    over the density over where mach number is equal to one for a non-adiabatic
    constant area frictionless system. Default fluid is air.

    Parameters
    ----------
    mach : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The ratio of rho / rho* \n

    Examples
    --------

    """

    gamma = gas.gamma
    rho_rhostar = (1+gamma*mach**2)/((1+gamma)*mach**2)
    return rho_rhostar



#==================================================
#rayleigh stagnation pressure star ratio
#==================================================
def rayleigh_stagnation_pressure_star_ratio(mach: float, gas=air) -> float:
    """Return the ratio of stagnation pressure over stagnation pressure where 
    Mach is equal to one

    Notes
    -----
    Given a mach number, return the ratio of the stagnation pressure of the mach number
    over the stagnation pressure where mach number is equal to one for a non-adiabatic
    frictionless system. Default fluid is air.

    Parameters
    ----------
    mach : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The ratio of pt / pt* \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 2
    >>> pt_ptstar = gd.rayleigh_stagnation_pressure_star_ratio(M)
    >>> pt_ptstar
    1.5030959785260414
    >>>
    """

    gamma = gas.gamma
    pt_ptstar = (1+gamma) / (1 + gamma*mach**2) * ((1 + (gamma-1)/2 * mach**2)/((gamma+1)/2))**(gamma/(gamma-1))
    return pt_ptstar



#==================================================
#rayleigh stagnation temperature star ratio
#==================================================
def rayleigh_stagnation_temperature_star_ratio(mach: float, gas=air) -> float:
    """Return the ratio of stagnation temperature over stagnation temperature where 
    Mach is equal to one

    Notes
    -----
    Given a mach number, return the ratio of the stagnation temperature of the mach number
    over the stagnation temperature where mach number is equal to one for a non-adiabatic
    frictionless constant area system. Default fluid is air.

    Parameters
    ----------
    mach : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Returns
    -------
    float
        The ratio of Tt / Tt* \n

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> M = 2
    >>> Tt_Ttstar = gd.rayleigh_stagnation_temperature_star_ratio(M)
    >>> Tt_Ttstar
    0.793388429752066
    >>>
    """

    gamma = gas.gamma
    Tt_Ttstar = (2*(gamma+1)*mach**2)/((1+gamma*mach**2)**2) * (1 + (gamma-1)/2 * mach**2)
    return Tt_Ttstar


#==================================================
#rayleigh heat flux
#==================================================
def rayleigh_heat_flux(stagnation_temperature_initial: float, stagnation_temperature_final: float, gas=air) -> float:
    """Return the heat per unit mass in our out given the two stagnation temperatures and the fluid

    Notes
    -----
    Given the initial and final stagnation temperatures, determine the specific heat
    flux to satisfy the constant area frictionless flow system. Default fluid is air.

    Parameters
    ----------
    Tt1 : `float`
        The stagnation temperature at region 1\n
    Tt2 : `float`
        The stagnation temperature at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Returns
    -------
    float
        The change in specific heat. \n
    
    Examples
    --------
    >>> from gas_dynamics.fluids import air 
    >>> air.cp = 1000
    >>> Tt1, Tt2 = 280, 107.9
    >>> gd.rayleigh_heat_flux(Tt1=Tt1, Tt2=Tt2, ,air) 
    -172100.0
    >>>
    """

    cp = gas.cp
    q = cp*(stagnation_temperature_final-stagnation_temperature_initial)
    return q