###This script is a large list of functions from the oscar biblarz Gas Dynamics book
# I wrote this to be able to import the following functions into a script and to further my understanding of gas dynamics
# and how the equations are derived. Maybe also some layperson definitions of what is happening to help undersanding. Its also just for fun!
import numpy as np
import matplotlib.pyplot as plt


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
    return dV


def T_stgn(M,gamma = 1.4):
    """ Given Mach number and gamma, returns the relation of T / Tt

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    
    Tt_ratio = 1 / (1+(gamma-1)/2 *M**2)
    return Tt_ratio


def P_stgn(M,gamma = 1.4):
    """ given Mach number and gamma, returns the relation of P / Pt
    
    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    denom = 1 + (gamma-1)/2 * M**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio


def Area_stgn(M,gamma = 1.4):
    """given Mach number and gamma, returns the relation of A / A*

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    A_star_ratio = 1/M * ((1 + (gamma-1)/2 * M**2) / ((gamma+1)/2)) ** ((gamma+1)/(2*(gamma-1)))
    return A_star_ratio


def rho_stgn(M,gamma = 1.4):
    """given Mach number and gamma, returns the relation of rho / rho_t

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M**2 )) ** (1 / (gamma-1))
    return rho_t_ratio


def plot_stgn_ratios(min=.1,max=10,increment=.01,gamma=[1.4]):
    """Plot all Mach number vs T/T, P/Pt, A/A*, rho/rho_t, for a given gamma

    :param min: min mach number
    :param max: max mach number
    :param increment: increment between min and max
    :param gamma: default gamma
    """
    fig, axs = plt.subplots(2,2)
    title = "Isentropic Stagnation Relations"
    fig.suptitle(title)

    for g in gamma:
        mach_nums = [i for i in np.arange(min,max,increment)]
        t_list = [T_stgn(i,g)for i in mach_nums]
        p_list = [P_stgn(i,g) for i in mach_nums]
        a_list = [Area_stgn(i,g) for i in mach_nums]
        rho_list = [rho_stgn(i,g) for i in mach_nums]
        labl = '\u03B3 = ' + str(g)

    #T/Tt
        axs[0,0].plot(mach_nums,t_list,label=labl)
        axs[0,0].set_xlabel('Mach Number')
        axs[0,0].set_ylabel('T / Tt')
        axs[0,0].grid(b=True, which='major')
        axs[0,0].grid(b=True, which='minor')
        axs[0,0].minorticks_on()
        axs[0,0].legend()

    #P/Pt
        axs[0,1].plot(mach_nums,p_list,label=labl)
        axs[0,1].set_xlabel('Mach Number')
        axs[0,1].set_ylabel('P / Pt')
        axs[0,1].grid(b=True, which='major')
        axs[0,1].grid(b=True, which='minor')
        axs[0,1].minorticks_on()
        axs[0,1].legend()
    #A/A*
        axs[1,0].plot(mach_nums,a_list,label=labl)
        axs[1,0].set_xlabel('Mach Number')
        axs[1,0].set_ylabel('A / A*')
        axs[1,0].grid(b=True, which='major')
        axs[1,0].grid(b=True, which='minor')
        axs[1,0].minorticks_on()
        axs[1,0].legend()

    #rho/rho_t
        axs[1,1].plot(mach_nums,rho_list,label=labl)
        axs[1,1].set_xlabel('Mach Number')
        axs[1,1].set_ylabel('rho / rho_t')
        axs[1,1].grid(b=True, which='major')
        axs[1,1].grid(b=True, which='minor')
        axs[1,1].minorticks_on()
        axs[1,1].legend()
    plt.show()


def print_stgn_ratios(min=.1,max=10,increment=.1,gamma = [1.4]):
    """Print all Mach number vs T/T, P/Pt, A/A*, rho/rho_t, for a given gamma or list of gammas 

    ex: to print several gammas,  
    :param min: min mach number
    :param max: max mach number
    :param increment: increment between min and max
    :param gamma: default gamma
    """
    for g in gamma:
        mach_nums = [i for i in np.arange(min,max,increment)]
        t_list = [T_stgn(i,g)for i in mach_nums]
        p_list = [P_stgn(i,g) for i in mach_nums]
        a_list = [Area_stgn(i,g) for i in mach_nums]
        rho_list = [rho_stgn(i,g) for i in mach_nums]

        for index, num in enumerate(mach_nums):
            print(" M: " + f"{num:.2f}" + "   |"+"   P/Pt: " + f"{p_list[index]:.3f}" + "   | " + "   T/Tt: " + f"{t_list[index]:.3f}" + "   | "+"   A/A*: " + f"{a_list[index]:.3f}" + "   |"+"   rho/rho_t: " + f"{rho_list[index]:.3f}")

        print("\n \n \n")

def func(x,y):
    """func descrioption Energy Equation
    -
    Lots of words
    -
    :param x: this is x
    :param y: this is y
    ]"""
    #do lots of stuff and calculate stuff
    print(x,y)

    return z