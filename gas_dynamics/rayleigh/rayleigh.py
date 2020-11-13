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
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_pressure_ratio(M1: float, M2:float, gas=air) -> float:
    """Return the pressure ratio p2 / p1 given the two Mach numbers

    Notes
    -----

    Parameters
    ----------
    M1 : `flaot`
        The Mach number at region 1\n
    M2 : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """
    gamma = gas.gamma
    p2_p1 = (1 + gamma * M1**2)/(1 + gamma *  M2**2)
    return p2_p1



#==================================================
#rayleigh temperature ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_temperature_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the temperature ratio T2/T1 given the two Mach numbers

    Notes
    -----

    Parameters
    ----------
    M1 : `flaot`
        The Mach number at region 1\n
    M2 : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """

    gamma = gas.gamma
    T2_T1 = ((1 + gamma * M1**2)/(1 + gamma *  M2**2))**2 * (M2**2)/(M1**2)
    return T2_T1



#==================================================
#rayleigh density ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_density_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the density ratio rho2/rho1 given the two Mach numbers

    Notes
    -----

    Parameters
    ----------
    M1 : `flaot`
        The Mach number at region 1\n
    M2 : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """

    gamma = gas.gamma
    rho2_rho1 = (1 + gamma * M2**2)/(1 + gamma *  M1**2) * (M1**2)/(M2**2)
    return rho2_rho1



#==================================================
#rayleigh stagnation temperature ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_stagnation_temperature_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the stagnation temperature ratio Tt2/Tt1 given the two Mach numbers

    Notes
    -----

    Parameters
    ----------
    M1 : `flaot`
        The Mach number at region 1\n
    M2 : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """
    gamma = gas.gamma
    Tt2_Tt1 = ((1+gamma*M1**2)/(1+gamma*M2**2))**2 * (M2**2)/(M1**2) * ((1+(gamma-1)/2*M2**2)/(1+(gamma-1)/2*M1**2))
    return Tt2_Tt1



#==================================================
#rayleigh stagnation pressure ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_stagnation_pressure_ratio(M1: float, M2: float, gas=air) -> float:
    """Return the stagnation pressure ratio pt2/pt1 given the two Mach numbers

    Notes
    -----

    Parameters
    ----------
    M1 : `flaot`
        The Mach number at region 1\n
    M2 : `float`
        The Mach number at region 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """

    gamma = gas.gamma
    pt2_pt1 = (1+gamma*M1**2)/(1+gamma*M2**2) * ((1+(gamma-1)/2*M2**2)/(1+(gamma-1)/2*M1**2))**(gamma/(gamma-1))
    return pt2_pt1



#==================================================
#rayleigh mach from pressure ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_mach_from_pressure_ratio(M: float, p1: float, p2: float, gas=air) -> float:
    """Return the mach number given the Mach number and two pressures

    Notes
    -----

    Parameters
    ----------
    M : `float`
        The Mach number\n
    p1 : `flaot`
        pressure 1\n
    p2 : `float`
        pressure 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """

    gamma = gas.gamma
    M2 = (((p2*(1+gamma*M**2))/p2 -1)/gamma)**.5
    return M2



#==================================================
#rayleigh mach from temperature ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_mach_from_temperature_ratio(M: float, T1: float, T2: float, gas=air) -> float:
    """Return the Mach number given the Mach number and two temperatures

    Notes
    -----

    Parameters
    ----------
    M : `float`
        The Mach number\n
    T1 : `flaot`
        Temperature 1\n
    T2 : `float`
        Temperature 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------
    
    """

    gamma = gas.gamma
    T2_T1 = T2/T1
    def zero(M2, M1, T2_T1,gas):
        return rayleigh_temperature_ratio(M1=M1, M2=M2, gas=gas) - T2_T1
    
    sol = fsolve(zero, args=(M, T2_T1, gas), x0=1)
    return sol[0]



#==================================================
#rayleigh mach from stagnation temperature ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_mach_from_stagnation_temperature_ratio(M: float, Tt1: float, Tt2: float, gas=air) -> float:
    """Return the Mach number given the Mach number and two stagnation temperatures

    Notes
    -----

    Parameters
    ----------
    M : `float`
        The Mach number\n
    Tt1 : `float`
        Stagnation temperature 1\n
    Tt2 : `float`
        Stagnation temperature 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------
    
    """

    gamma = gas.gamma
    Tt2_Tt1 = Tt2/Tt1
    def zero(M2, M1, Tt2_Tt1, gas):
        return rayleigh_stagnation_temperature_ratio(M1=M1, M2=M2, gas=gas) - Tt2_Tt1
    
    sol = fsolve(zero, args=(M, Tt2_Tt1, gas), x0=1)
    return sol[0]



#==================================================
#rayleigh mach from stagnation pressure ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_mach_from_stagnation_pressure_ratio(M: float, pt1: float, pt2: float, gas=air) -> float:
    """Return the Mach number given the Mach number and two stagnation pressures

    Notes
    -----

    Parameters
    ----------
    M : `float`
        The Mach number\n
    pt1 : `float`
        Stagnation pressure 1\n
    pt2 : `float`
        Stagnation pressure 2\n
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------
    
    """

    gamma = gas.gamma
    pt2_pt1 = pt2/pt1
    def zero(M2, M1, pt2_pt1, gas):
        return rayleigh_stagnation_pressure_ratio(M1=M1, M2=M2, gas=gas) - pt2_pt1

    sol = fsolve(zero, args=(M, pt2_pt1, gas), x0=1)
    return sol[0]




#==================================================
#rayleigh pressure star ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_pressure_star_ratio(M: float, gas=air) -> float:
    """

    Notes
    -----

    Parameters
    ----------
    M : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------


    """

    gamma = gas.gamma
    p_pstar = (gamma+1)/(1+gamma*M**2)
    return p_pstar



#==================================================
#rayleigh temperature star ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_temperature_star_ratio(M: float, gas=air) -> float:
    """

    Notes
    -----

    Parameters
    ----------
    M : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------

    """

    gamma = gas.gamma
    T_Tstar = (M**2 * (1+gamma)**2) / (1+gamma*M**2)**2
    return T_Tstar



#==================================================
#rayleigh density star ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_density_star_ratio(M: float, gas=air) -> float:
    """

    Notes
    -----

    Parameters
    ----------
    M : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------
    """

    gamma = gas.gamma
    rho_rhostar = (1+gamma*M**2)/((1+gamma)*M**2)
    return rho_rhostar



#==================================================
#rayleigh stagnation pressure star ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_stagnation_pressure_star_ratio(M: float, gas=air) -> float:
    """

    Notes
    -----

    Parameters
    ----------
    M : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------

    """

    gamma = gas.gamma
    pt_ptstar = (1+gamma) / (1 + gamma*M**2) * ((1 + (gamma-1)/2 * M**2)/((gamma+1)/2))**(gamma/(gamma-1))
    return pt_ptstar



#==================================================
#rayleigh stagnation temperature star ratio
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_stagnation_temperature_star_ratio(M: float, gas=air) -> float:
    """

    Notes
    -----

    Parameters
    ----------
    M : `flaot`
        The Mach number\n
    gas : `fluid`
        A user defined fluid object. Default is air \n

    Examples
    --------


    """

    gamma = gas.gamma
    Tt_Ttstar = (2*(gamma+1)*M**2)/((1+gamma*M**2)**2) * (1 + (gamma-1)/2 * M**2)
    return Tt_Ttstar


#==================================================
#rayleigh heat flux
#TODO: docstring and examples
#TODO: verify
#==================================================
def rayleigh_heat_flux(Tt1: float, Tt2: float, gas=air) -> float:
    """Return the heat per unit mass in our out given the two stagnation temperatures and the fluid

    Notes
    -----

    Parameters
    ----------
    Tt1 : `flaot`
        The stagnation temperature at region 1
    Tt2 : `float`
        The stagnation temperature at region 2
    gas : `fluid`
        A user defined fluid object. Default is air \n
    
    Examples
    --------

    """
    cp = gas.cp
    q = cp*(Tt2-Tt1)
    return q