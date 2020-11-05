#!usr/bin/env
"""
Extra functions for random stuff
"""

import numpy as np



#==================================================
#fluid
#remove custom feature
# add doc for '$1.4,286'
#==================================================
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
        elif fluid[0] == '$':
            new_fluid = fluid[1:]
            gamma = float(new_fluid.split(',')[0])
            R = float(new_fluid.split(',')[1])            
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



#==================================================
#degrees
#==================================================
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



#==================================================
#radians
#==================================================
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



#==================================================
#sine degrees
#==================================================
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
    


#==================================================
#arc sine degrees
#==================================================
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



#==================================================
#cosine degrees
#==================================================
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



#==================================================
#arc cosine degrees
#==================================================
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



#==================================================
#tangent degrees
#==================================================
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



#==================================================
#arc tangent degrees
#==================================================
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



#==================================================
#are to diameter
#TODO: deprecate and split into two
#==================================================
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



#==================================================
#linear interpolate
#==================================================
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