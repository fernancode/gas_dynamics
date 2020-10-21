import numpy as np
import matplotlib.pyplot as plt

#TODO: add error handling, automatic switching for get='T1' to get='T2'
#TODO: add shock tables generator and shock functions
#TODO: add docstrings for functions
#TODO: add a bunch of values for R and gamma




#R values for various gases
R_air = 287



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


def T_stgn_ratio(M,gamma = 1.4):
    """ Given Mach number and gamma, returns the relation of T / Tt

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    
    Tt_ratio = 1 / (1+(gamma-1)/2 *M**2)
    return Tt_ratio


def p_stgn_ratio(M,gamma = 1.4):
    """ given Mach number and gamma, returns the relation of P / Pt
    
    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    denom = 1 + (gamma-1)/2 * M**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio


def a_star_ratio(M=[], A_ratio=[], gamma = 1.4):
    """given Mach number and gamma, returns the relation of A / A*

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    if A_ratio==[]:
        A_star_ratio = 1/M * ((1 + (gamma-1)/2 * M**2) / ((gamma+1)/2)) ** ((gamma+1)/(2*(gamma-1)))
        return A_star_ratio
    
    if M==[]:
        t1 = (gamma+1) / (2*(gamma-1))
        t2 = -1*(gamma+1)/(2*(gamma-1))
        M = ( t1 * A_ratio**t2) ** ((gamma-1)/(gamma+1))
        return M


def rho_stgn_ratio(M,gamma = 1.4):
    """given Mach number and gamma, returns the relation of rho / rho_t

    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    """
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M**2 )) ** (1 / (gamma-1))
    return rho_t_ratio


def plot_stgn_ratios(Mach_min=.01,Mach_max=5,increment=.01,gamma=[1.4]):
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
        mach_nums = [i for i in np.arange(Mach_min,Mach_max+increment,increment)]
        t_list = [T_stgn_ratio(i,g)for i in mach_nums]
        p_list = [p_stgn_ratio(i,g) for i in mach_nums]
        a_list = [a_star_ratio(M=i, gamma=g) for i in mach_nums]
        rho_list = [rho_stgn_ratio(i,g) for i in mach_nums]
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


def print_stgn_ratios(Mach_min=0,Mach_max=5,increment=.1,gamma = [1.4]):
    """Print all Mach number vs T/T, P/Pt, A/A*, rho/rho_t, for a given gamma or list of gammas 

    ex: to print several gammas,  
    :param min: min mach number
    :param max: max mach number
    :param increment: increment between min and max
    :param gamma: default gamma
    """
    for g in gamma:
        mach_nums = [i for i in np.arange(Mach_min,Mach_max+increment,increment)]
        t_list = [T_stgn_ratio(i,g)for i in mach_nums]
        p_list = [p_stgn_ratio(i,g) for i in mach_nums]
        a_list = [a_star_ratio(M=i,gamma=g) for i in mach_nums]
        rho_list = [rho_stgn_ratio(i,g) for i in mach_nums]

        labl = '\u03B3 = ' + str(g)
        print("Isentropic Flow Parameters for " + labl)
        for index, num in enumerate(mach_nums):
            print("M: " + f"{num:.2f}" + "   |"+"   P/Pt: " + f"{p_list[index]:.3f}" + "   | " + "   T/Tt: " + f"{t_list[index]:.3f}" + "   | "+"   A/A*: " + f"{a_list[index]:.3f}" + "   |"+"   rho/rho_t: " + f"{rho_list[index]:.3f}")
        print("\n \n \n")


def pressure_mach_ratio(p1=[],p2=[],M1=[],M2=[],get='p2',gamma=1.4,R=286.9,ds=0):
    """ Specify whether you need Mach number or Pressure, and provide the three knowns ex. get = 'P2', M1, M2, P1 will return the missing pressure. default arguments are gamma = 1.4, R = 286 , ds = 0

    """
    if get == 'p2':
        p2 = p1 * ((1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2))**(gamma/(gamma-1)) * np.exp(-ds/R)
        return p2

    elif get =='M2':
        M2 = (((p1/p2 * np.exp(ds/R))**((gamma-1)/gamma) * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
        return M2

    elif get == 'ds':
        ds = R * np.log( p1/p2 * ((1 + (gamma-1)/2 * M1**2 )/(1 + (gamma-1)/2 * M2**2))**(gamma/(gamma-1)))
        return ds
    else:
        print('Incorrect argument')


def entropy_produced(pt1=[], pt2=[], R=[]):
    ds = -R*np.log(pt2/pt1)
    return ds


def temperature_mach_ratio(T1=[],T2=[],M1=[],M2=[],get='T2',gamma=1.4):
    """Specify whether you need Mach number or Temperature, and provide the three knowns ex. get = 'T2', M1, M2, T1 will return the missing temperature. Default arguments are gamma = 1.4

    """
    if get == 'T2':
        T2 = T1 * (1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2)
        return T2

    elif get == 'M2':
        M2 = (( T1/T2 * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
        return M2

    else:
        print('Incorrect argument')


def area_mach_ratio(M1,M2,gamma=1.4,R=286.9,ds=0):
    """Specify whether you need Mach number or Area, and provide the three knowns ex. get = 'A2', M1, M2, A1 will return the missing Area. Default arguments are gamma = 1.4, R = 286.9, ds = 0. 

    """
    A2_A1 = M1/M2 * ((1 + (gamma-1)/2 * M2**2 )/(1 + (gamma-1)/2 * M1**2 ))**((gamma+1)/(2*(gamma-1))) * np.exp(ds/R)
    return A2_A1


def stgn_pressure(p=[], M=[], p_t=[], gamma=1.4, get='p_t'):
    if get == 'p_t':
        p_t = p* ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        return p_t
   
    elif get == 'M':
        M = (((p_t/p)**((gamma-1)/gamma) -1 ) * 2/(gamma-1) ) ** .5
        return M

    elif get == 'p':
        p = p_t / ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        return p

    else:
        print('Incorrect argument')


def sonic_velocity(gamma=1.4,R=287,T=273.15):
    a = (gamma*R*T)**.5
    return a


def stgn_temperature(T_t =[], T=[] , M=[], get = 'T_t', gamma = 1.4):
    if get == 'T_t':
        T_t = T * ( 1 + (gamma-1)/2 * M**2)
        return T_t

    elif get == 'M':
        M = ((T_t /T - 1) * 2/(gamma-1))**.5
        return M
        
    elif get =='T':
        T = T_t/( 1 + (gamma-1)/2 * M**2)
        return T

    else:
        print('Incorrect argument')


def mdot_a_star(p_t=[], T_t=[], R=286.9, gamma=1.4):
    """ Returns the maximum flow rate over a_star (m_dot/a_star)
    Check your units for pressure! Need to be in Pa!
    
    """
    mdot_a_star = (((gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * p_t/(T_t**.5))
    return mdot_a_star


def shock_tables(Mach_min=1, Mach_max=5, increment=.05, gamma=1.4):
    '''
    Generate shock tables for a range of Mach numbers for a given gamma
    '''
    mach_nums = [i for i in np.arange(Mach_min,Mach_max+increment,increment)]
    M2 = [shock_mach(i,gamma)for i in mach_nums]
    p2_p1 = [shock_pressure_ratio(i,gamma)for i in mach_nums]
    T2_T1 = [shock_temperature_ratio(i,gamma) for i in mach_nums]
    dv_a = [shock_dv_a(i,gamma) for i in mach_nums]
    pt2_pt1 = [shock_stagnation_ratio(i,gamma) for i in mach_nums]
    
    labl = '\u03B3 = ' + str(gamma)
    print("Normal Shock Parameters for " + labl)
    for index, num in enumerate(mach_nums):
        print("M: " + f"{num:.2f}" + "   |"+"   M2: " + f"{M2[index]:.4f}" + "   | " + "   p2/p1: " + f"{p2_p1[index]:.4f}" + "   | "+"   T2/T1: " + f"{T2_T1[index]:.4f}" + "   |"+"   dV/a: " + f"{dv_a[index]:.4f}"+ "   |"+"   pt2/pt1: " + f"{pt2_pt1[index]:.6f}" )
    print("\n \n \n")



def shock_mach(M1=[], M2=[], gamma=1.4):
    if M2==[]:
        M2 = ((M1**2 + 2/(gamma-1)) / ((2*gamma / (gamma-1)) * M1**2 - 1))**.5
        return M2

    if M1==[]:
        M1 = ((-2/(gamma-1) -M2**2 ) / (1- ((2*gamma)/(gamma-1))*M2**2))**.5
        return M1

def shock_pressure_ratio(M=[], p2_p1=[], gamma=1.4):
    '''
    Gives p2/p1 after a standing normal shock for a given Mach number
    '''
    if p2_p1 == []:
        p2_p1 = 2*gamma / (gamma+1) * M**2 - (gamma-1)/(gamma+1)
        return p2_p1

    if M == []:
        M = ((gamma+1)/(2*gamma) * (p2_p1 + (gamma-1)/(gamma+1)) )**.5
        return M


def shock_temperature_ratio(M1=[], gamma=1.4):
    '''
    Gives T2/T1 after a standing normal shock for a given Mach number
    '''
    term1 = (1 + (gamma-1)/2 * M1**2)
    term2 = (2*gamma)/(gamma-1) * M1**2 -1
    term3 = (gamma+1)**2 / (2*(gamma-1)) * M1**2
    t2_t1 = (term1 * term2) / term3
    return t2_t1

def shock_dv_a(M1=[], gamma=1.4):
    '''
    Gives velocity change across a standing normal shock for a given Mach number
    '''
    dv_a = 2/(gamma+1) * (M1**2 -1)/ M1
    return dv_a


def shock_stagnation_ratio(M1=[], gamma=1.4):
    '''
    Gives pt2/pt1 across a standing normal shock for a given Mach number
    '''
    term1 = (gamma+1)/2*M1**2
    term2 = 1 + (gamma-1)/2 * M1**2
    term3 = (term1/term2) ** (gamma/(gamma-1))
    term4 = (2*gamma / (gamma+1) * M1**2 - ((gamma-1)/(gamma+1)))**(1/(1-gamma))
    return term3 * term4


def shock_oblique_charts(gamma=1.4):
    
    #generate values and plot them
    theta = np.linspace(0,np.pi/2, 100)
    mach = np.linspace(1,6,100)
    MACH,THETA = np.meshgrid(mach,theta) 
    dirac = np.arctan( 2 * 1/np.tan(THETA) * (MACH**2 * np.sin(THETA)**2 - 1 ) / (MACH**2 * (gamma + np.cos(2*THETA)) + 2 ))
    
    h = plt.contourf(mach,theta,dirac)
    plt.show()


def func(x,y):
    """func descrioption Energy Equation
    -
    Lots of words
    -
    :param x: this is x
    :param y: this is y
    ]
    """
    #do lots of stuff and calculate stuff
    print(x,y)

    return z