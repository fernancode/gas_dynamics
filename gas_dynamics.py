import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

#TODO: add error handling, automatic switching for get='T1' to get='T2'
#TODO: add a bunch of values (in a class? object?) for R and gamma


#R values for various gases
R_air = 287
#TODO: function where i gave any function gas='air' (metric=true) default, and that function in the other functions grabs the R and gamma values from an object that has all of the values for both metric and std.


def plot_stagnation_ratios(range=[.1,5],inc=.01,gamma=[1.2, 1.3, 1.4]):
    '''Plots Mach # vs T/T, P/Pt, A/A*, rho/rho_t, for a list of specific heat ratios.

    Parameters:
    :param range: The starting and ending Mach # in a list, ex: [.01,5]
    :param inc: The increment between min and max
    :param gamma: A list of the ratio of specific heats to be plotted, ex [1.2, 1.3, 1.4]
    '''
    fig, axs = plt.subplots(2,2)
    title = "Isentropic Stagnation Relations"
    fig.suptitle(title)

    Mach_min = range[0]
    Mach_max = range[1]

    for g in gamma:
        mach_nums = [i for i in np.arange(Mach_min,Mach_max+inc,inc)]
        t_list = [stagnation_temperature_ratio(M=i, gamma=g)for i in mach_nums]
        p_list = [stagnation_pressure_ratio(M=i, gamma=g) for i in mach_nums]
        a_list = [mach_area_choked_ratio(M=i, gamma=g) for i in mach_nums]
        rho_list = [stagnation_density_ratio(M=i, gamma=g) for i in mach_nums]
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


def stagnation_ratios(range=[0,5], inc=.1, gamma=1.4):
    '''Returns the isentropic flow relation tables for all Mach #'s in the given range.
    
    Given a ratio of specific heats, print out the stagnation temperature ratio, stagnation pressure ratio, the area to choked area ratio, and the stagnation density ratio for every incremental Mach #.

    Parameters:
    :param range: The starting and ending Mach # in a list, ex: [.01,5]
    :param inc: The increment between min and max
    :param gamma: The ratio of specific heats
    '''
    Mach_min = range[0]
    Mach_max = range[1]


    mach_nums = [i for i in np.arange(Mach_min,Mach_max+inc,inc)]
    t_list = [stagnation_temperature_ratio(M=i, gamma=gamma)for i in mach_nums]
    p_list = [stagnation_pressure_ratio(M=i, gamma=gamma) for i in mach_nums]
    a_list = [mach_area_choked_ratio(M=i, gamma=gamma) for i in mach_nums]
    rho_list = [stagnation_density_ratio(M=i, gamma=gamma) for i in mach_nums]

    labl = '\u03B3 = ' + str(gamma)
    print("Isentropic Flow Parameters for " + labl)
    for index, num in enumerate(mach_nums):
        print('M: %0.3f' % num, '  |   P/Pt: %0.3f' % p_list[index], '   |    T/Tt: %0.3f' % t_list[index],  '   |    A/A*: %0.3f' % a_list[index],  '   |   rho/rho_t: %0.3f ' % rho_list[index])
    print("\n \n \n")


def sonic_velocity(gamma=1.4,R=287,T=273.15):
    '''Returns the local speed of sound.

    Given a ratio of specific heats, gas constant, and temperature, this function returns the locoal speed of sound.
    
    :param gamma: The ratio of specific heats.
    :param R: The gas constant.
    :param T: The temperature
    '''
    a = (gamma*R*T)**.5
    return a


def entropy_produced(pt1=[], pt2=[], R=[]):
    ds = -R*np.log(pt2/pt1)
    return ds


def mach_pressure_ratio(p1=[],p2=[],M1=[],M2=[],get='p2',gamma=1.4,R=286.9,ds=0):
    '''
    #TODO: edit this docstring
     Specify whether you need Mach number or Pressure, and provide the three knowns ex. get = 'P2', M1, M2, P1 will return the missing pressure. default arguments are gamma = 1.4, R = 286 , ds = 0

    '''
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


def mach_temperature_ratio(T1=[],T2=[],M1=[],M2=[],get='T2',gamma=1.4):
    '''
    #TODO:edit this docstring
    Specify whether you need Mach number or Temperature, and provide the three knowns ex. get = 'T2', M1, M2, T1 will return the missing temperature. Default arguments are gamma = 1.4

    '''
    if get == 'T2':
        T2 = T1 * (1 + ((gamma-1)/2) *M1**2)/(1 + ((gamma-1)/2) *M2**2)
        return T2

    elif get == 'M2':
        M2 = (( T1/T2 * (1 + (gamma-1)/2 * M1**2) - 1) * 2/(gamma-1))**0.5
        return M2

    else:
        print('Incorrect argument')


def mach_area_ratio(M1,M2,gamma=1.4,R=286.9,ds=0):
    '''
    #TODO:edit this docstring

    '''
    A2_A1 = M1/M2 * ((1 + (gamma-1)/2 * M2**2 )/(1 + (gamma-1)/2 * M1**2 ))**((gamma+1)/(2*(gamma-1))) * np.exp(ds/R)
    return A2_A1


def mach_area_choked_ratio(M=[], a_ratio=[], gamma = 1.4):
    '''Returns the ratio of A / A* given the Mach #.
    
    Given the Mach # and the ratio of specific heats, return the ratio of the area at the current Mach # to the area at which Mach # = 1. Default ratio of specific heats is for air.

    Parameters:
    :param M: Mach Number
    :param gamma: ratio of specific heats (default 1.4)
    '''
    a_star_ratio = 1/M * ((1 + (gamma-1)/2 * M**2) / ((gamma+1)/2)) ** ((gamma+1)/(2*(gamma-1)))
    return a_star_ratio


def stagnation_pressure(p=[], M=[], pt=[], gamma=1.4):
    '''Returns the stagnation pressure given pressure and Mach #.

    Given a pressure, Mach #, and a ratio of specific heats, return the stagnation pressure. Alternatively, provided two arguments the function will return the missing one. Default ratio of specific heats is for air.

    Parameters:
    :param pt: The stagnation pressure.
    :param p: The pressure.
    :param M: The Mach #
    :param gamma: The ideal gas constant 
    '''
    if pt == []:
        pt = p* ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        return pt
   
    if M == []:
        M = (((pt/p)**((gamma-1)/gamma) -1 ) * 2/(gamma-1) ) ** .5
        return M

    if p == []:
        p = pt / ( 1 + (gamma-1)/2 * M**2)** (gamma/(gamma-1))
        return p


def stagnation_pressure_ratio(M, gamma = 1.4):
    """Returns the pressure ratio of p / p_t

    Given a Mach # and ratio of specific heats, return the relation of p / p_t. Default ratio of specific heats is for air.

    Parameters:
    :param M: The Mach #
    :param gamma: The ratio of specific heats
    """
    denom = 1 + (gamma-1)/2 * M**2
    Pt_ratio = (1 / denom ) ** (gamma/(gamma-1))
    return Pt_ratio


def stagnation_temperature(Tt =[], T=[] , M=[], gamma = 1.4):
    '''Returns the stagnation temperature given temperature and Mach #.

    Given a temperature, Mach #, and a ratio of specific heats, this function returns the stagnation temperature. Alternatively, provided two arguments the function will return the missing one. Default ratio of specific heats is for air.

    Parameters:
    :param Tt: The stagnation temperature.
    :param T: The temperature.
    :param M: The Mach #
    :param gamma: The ratio of specific heats
    '''

    if Tt == []:
        Tt = T * ( 1 + (gamma-1)/2 * M**2)
        return Tt

    if M == []:
        M = ((Tt /T - 1) * 2/(gamma-1))**.5
        return M
        
    if T == []:
        T = Tt/( 1 + (gamma-1)/2 * M**2)
        return T


def stagnation_temperature_ratio(M, gamma = 1.4):
    """Returns the temperature ratio of T / T_t

    Given a Mach # and ratio of specific heats, return the relation of T / T_t.  Default ratio of specific heats is for air.

    Parameters:
    :param M: The Mach #
    :param gamma: The ratio of specific heats
    """
    
    Tt_ratio = 1 / (1+(gamma-1)/2 *M**2)
    return Tt_ratio


def stagnation_density_ratio(M, gamma = 1.4):
    """Returns the density ratio of rho/ rho_t

    Given a Mach # and ratio of specific heats, return the relation of rho / rho_t. Default ratio of specific heats is for air.

    Parameters:
    :param M: The Mach #
    :param gamma: The ratio of specific heats
    """
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M**2 )) ** (1 / (gamma-1))
    return rho_t_ratio


def choked_mdot(pt=[], Tt=[], R=286.9, gamma=1.4, metric=True):
    """Returns the maximum flow rate over A* (m_dot/a_star)

    Given stagnation pressure, stagnation temperature, the gas constant, and ratio of specific heats, return the given flow rate per unit A*, or flow rate for a given choked flow area. Default gas constant and ratio of specific heats are for air. 
    
    Check your units! metric units need to be in Pa
    #TODO: figure out what the std units output are
    
    Parameters:
    :param pt: The stagnation pressure.
    :param Tt: The stagnation temperature.
    :param R: The gas constant J / kg K
    :param gamma: The ratio of specific heats
    """
    if metric == True:
        mdot_a_star = (((gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * pt/(Tt**.5))
        return mdot_a_star

    if metric == False:
        gc = 32.174 #lbm -ft / lbf -s^2
        mdot_a_star = (((gc*gamma/(R))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5 * pt/(Tt**.5))
        return mdot_a_star


def shock_tables(range=[1,5], inc=.01, gamma=1.4):
    '''Returns shock tables for a range of Mach #'s.

    Given a range of Mach #'s and a ratio of specific heats, generate the standing normal shock tables for every incremental Mach # in-between.

    Parameters:
    :param range: The starting and ending Mach # in a list, ie: [1,5].
    :param inc: The step size for the tables.
    :param gamma: The ratio of specific heats
    '''
    Mach_min = range[0]
    Mach_max = range[1]

    mach_nums = [i for i in np.arange(Mach_min,Mach_max+inc,inc)]
    M2 = [shock_mach(M1=i, gamma=gamma)for i in mach_nums]
    p2_p1 = [shock_pressure_ratio(M=i, gamma=gamma)for i in mach_nums]
    T2_T1 = [shock_temperature_ratio(M=i, gamma=gamma) for i in mach_nums]
    dv_a = [shock_dv_a(M=i, gamma=gamma) for i in mach_nums]
    pt2_pt1 = [shock_stagnation_ratio(M=i, gamma=gamma) for i in mach_nums]
    
    labl = '\u03B3 = ' + str(gamma)
    print("Normal Shock Parameters for " + labl)
    for index, num in enumerate(mach_nums):
        print("M: " + f"{num:.2f}" + "   |"+"   M2: " + f"{M2[index]:.4f}" + "   | " + "   p2/p1: " + f"{p2_p1[index]:.4f}" + "   | "+"   T2/T1: " + f"{T2_T1[index]:.4f}" + "   |"+"   dV/a: " + f"{dv_a[index]:.4f}"+ "   |"+"   pt2/pt1: " + f"{pt2_pt1[index]:.6f}" )
    print("\n \n \n")


def shock_mach(M1=[], M2=[], gamma=1.4):
    '''Returns the Mach # after a standing normal shock.

    Given a starting Mach # M1 and the ratio of specific heats, return the Mach # M2 that immediately follows the shock. If M1 is not specified and M2 is, returns the Mach # prior to the shock. Default ratio of specific heats is for air.

    Parameters:
    :param M1: The Mach # before the shock
    :param M2: The Mach # after the shock
    :param gamma: The ratio of specific heats
    '''
    if M2 == []:
        M2 = ((M1**2 + 2/(gamma-1)) / ((2*gamma / (gamma-1)) * M1**2 - 1))**.5
        return M2

    if M1 == []:
        M1 = ((-2/(gamma-1) -M2**2 ) / (1- ((2*gamma)/(gamma-1))*M2**2))**.5
        return M1

def shock_pressure_ratio(M=[], p2_p1=[], gamma=1.4):
    '''Returns the pressure ratio after a standing normal shock for a given Mach #
    
    Given a starting Mach # and a ratio of specific heats, this function returns the ratio of p2 / p1 across a standing normal shock. If Mach # is not specified and p2_p1 is, function returns Mach #. Default ratio of specific heats is for air.
    
    Parameters:
    :param M: The starting Mach #
    :param gamma: The ratio of specific heats
    '''
    if p2_p1 == []:
        p2_p1 = 2*gamma / (gamma+1) * M**2 - (gamma-1)/(gamma+1)
        return p2_p1

    if M == []:
        M = ((gamma+1)/(2*gamma) * (p2_p1 + (gamma-1)/(gamma+1)) )**.5
        return M


def shock_temperature_ratio(M=[], gamma=1.4):
    '''Returns the temperature ratio after a standing normal shock for a given Mach number
    
    Given a starting Mach # and a ratio of specific heats, this function returns the ratio of T2 / T1 across a standing normal shock. Default ratio of specific heats is for air.
    
    Parameters:
    :param M: The starting Mach #
    :param gamma: The ratio of specific heats
    '''
    term1 = (1 + (gamma-1)/2 * M**2)
    term2 = (2*gamma)/(gamma-1) * M**2 -1
    term3 = (gamma+1)**2 / (2*(gamma-1)) * M**2
    t2_t1 = (term1 * term2) / term3
    return t2_t1

def shock_dv_a(M=[], gamma=1.4):
    '''Returns change in velocity over the local speed of sound after a normal shock.
    
    Given a starting Mach # and a ratio of specific heats, this function returns the velocity change across a standing normal shock divided by the local speed of sound. Default ratio of specific heats is for air.
    
    Parameters:
    :param M: The starting Mach #
    :param gamma: The ratio of specific heats
    '''
    dv_a = 2/(gamma+1) * (M**2 -1)/ M
    return dv_a


def shock_stagnation_ratio(M=[], gamma=1.4):
    '''Returns stagnation pressure ratio after a normal shock.
    
    Given a starting Mach # and a ratio of specific heats, this function returns the ratio of pt2/pt1 across a standing normal shock. Default ratio of specific heats is for air.
    
    
    Parameters:
    :param M: The starting Mach #
    :param gamma: The ratio of specific heats
    '''
    term1 = (gamma+1)/2*M**2
    term2 = 1 + (gamma-1)/2 * M**2
    term3 = (term1/term2) ** (gamma/(gamma-1))
    term4 = (2*gamma / (gamma+1) * M**2 - ((gamma-1)/(gamma+1)))**(1/(1-gamma))
    return term3 * term4


def shock_oblique_charts(Mach_max=6,gamma=1.4):
    '''Generate 2-D Oblique Shock Charts

    Displays two plots,
    1) Mach # versus oblique shock wave angle and the corresponding deflection angles for each case.
    2) Mach # versus resulting Mach # after an oblique shock and the corresponding deflection angles for each case. 
    Default ratio of specific heats is for air.
    
    Parameters:
    :param Mach_max: The upper limit Mach # for the chart
    :param gamma: The ratio of specific heats
    '''

    #generate values and plot them
    theta = np.linspace(.001,np.pi/2, 1000)
    mach = np.linspace(1,Mach_max,1000)
    MACH,THETA = np.meshgrid(mach,theta) 
    dirac = np.arctan( 2 * 1/np.tan(THETA) * (MACH**2 * np.sin(THETA)**2 - 1 ) / (MACH**2 * (gamma + np.cos(2*THETA)) + 2 ))
    
    theta_deg = []
    [theta_deg.append(n * 180/np.pi) for n in theta]
    dirac_deg=[]
    [dirac_deg.append(n*180/np.pi) for n in dirac]
    
    fig, (ax1,ax2) = plt.subplots(1,2)
    levels=[5,10,15,20,25,30,35,40]
    h = ax1.contour(mach,theta_deg,dirac_deg,levels=levels,cmap='tab10')
    ax1.clabel(h, inline =1, fontsize=10)
    minor_ticks_mach = np.arange(1,Mach_max+.1,.1)
    minor_ticks_theta = np.arange(0,91,1)

    ax1.set_yticks(minor_ticks_theta, minor=True)
    ax1.set_xticks(minor_ticks_mach, minor=True)
    ax1.grid(which='major',color='k',linestyle = '--', alpha=.5,)
    ax1.grid(which='minor',color='k',linestyle = '--', alpha=.1)
    ax1.set(xlabel = 'Mach #')
    ax1.set(ylabel = 'Oblique Shock Wave Angle')
    ax1_1 = ax1.twinx()
    ax1_1.set(ylabel = 'Flow Deflection Angle')
    ax1_1.get_xaxis().set_visible(False)
    ax1_1.get_yaxis().set_visible(False)

    #TODO: the same table but y axis is Mach number after an oblique shock.

    fig.tight_layout(pad=2.0)
    plt.show()