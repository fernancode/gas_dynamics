###This script is a large list of functions from the oscar biblarz Gas Dynamics book
# I wrote this to be able to import the following functions into a script and to further my understanding of gas dynamics
# and how the equations are derived. Maybe also some layperson definitions of what is happening to help undersanding. Its also just for fun!

#TODO: add equations for perfect gases in terms of mach number and how to get them
#TODO: add stagnation equations and how to get them and what they mean ch3
#TODO: add stagnation equations from mach ch4




###The following three equations are for the general behavior of an ideal and arbitrary fluid with no losses
###They are used to understand the relation of increasing mach number on pressure, density, Area, and Velocity
# These equations can be derived from the following:
# 1. Energy equation w/assumptions steady 1D flow, adiabatic (dq = 0), no losses (dsi = 0), neglect potential energy (dz = 0), no shaft work (dws = 0)
# 2. Property relation Tds = dh - dp/rho
# 3. Continuity (mdot = constant)
# 4. Sonic velocity, a^2 = gc * dp/drho at constant entropy
# 5. Mach number M^2 = V^2/a^2
 
def dp_from_Mach(rho,V,M,da,A,metric=True):
    """
    Takes density, velocity, mach number, change in area, area to calculate change in pressre for an arbitrary fluid

    :param rho: Density
    :param V: Velocity
    :param M: Mach number
    :param dA: Change in Area
    :param A: Area    
    """

    if metric == True:
        dp = rho*V^2 * (1/(1-M^2)) * (dA/A)
        return dp
    else:
        gc = 32.174 #formatted like this so gravitational constant usage for imperial is clear.
        dp = (rho*V^2)/gc * (1/(1-M^2))*(dA/A)
        return dp

def drho_from_Mach(rho,M,dA,A):
    """ 
    Takes density, mach number, change in area, and area to calculate change in density for an arbitrary fluid

    :param rho: density
    :param M: Mach number
    :param dA: change in area
    :param A: Area    
    """

    drho = rho * M^2/(1-M^2) * dA/A
    return drho

def dv_from_Mach(V,M,dA,A):
    """Takes velocity, mach number, change in area, and area to calculate dv for an arbitrary fluid
    
    :param V: velocity
    :param M: Mach number
    :param dA: change in area
    :param A: area
    """

    dV = -V * (1/(1-M^2)) * (dA/A)

#TODO: perfect gas with losses
#TODO: * ref concept explainer


def Temp_stgn_ratio(M,gamma = 1.4):
    """ Given Mach number and gamma, returns the relation of T / Tt

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    
    Tt_ratio = 1 / (1+(gamma-1)/2*M**2)
    return Tt_ratio

def Pressure_stgn_ratio(M,gamma = 1.4):
    """ given Mach number and gamma, returns the relation of P / Pt
    
    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    denom = 1 + (gamma-1)/2 * M**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio

def Area_stgn_ratio(M,gamma = 1.4):
    """given Mach number and gamma, returns the realtion of A / A*

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    A_star_ratio = 1/M * ((1 + (gamma-1)/2 * M**2) / ((gamma+1)/2)) ** ((gamma+1)/(2*(gamma-1)))
    return A_star_ratio

def Density_stgn_ratio(M,gamma = 1.4):
    """given Mach number and gamma, returns the realtion of rho / rho_t

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M**2 )) ** (1 / (gamma-1))
    return rho_t_ratio


def func(x,y):
    """func descrioption Energy Equation
    -
    Lots of words
    -
    :param x: this is x
    :param y: this is y
    ]"""
    #do lots of stuff and calculate stuff
    x = x #bla
    y = y #bla
    z = x+y #bla

    return z