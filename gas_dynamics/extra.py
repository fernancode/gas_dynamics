#!usr/bin/env
"""
Extra functions for random stuff
"""

import numpy as np


#==================================================
#fluid
#==================================================
class fluid:
    """A class to represent the fluid and its properties

    Attributes
    ----------
    gamma : `float`
        The ratio of specific heats \n
    R : `float`
        The gas constant for the fluid \n
    units : `str`
        The unit system defining the gas constant \n

    Methods
    -------
    No methods at this time

    Examples
    --------
    >>> import gas_dynamics as gd
    >>> methane = gd.fluid('methane', 1.3, 518.2)
    >>> methane.gamma
    1.3
    >>> methane.R
    518.2
    >>> methane.units
    'metric'
    
    Conversely we can set the units
    >>> methane = fluid('methane', 1.3, 0.1238, units = 'btu / lbm-R')
    >>> methane.gamma
    1.3
    >>> methane.R
    0.1238
    >>> methane.units
    'btu / lbm-R'

    """

    def __init__(self, name: str, gamma: float, R: float, units='metric'):
        """Construct the necessary attributes for the fluid object

        Parameters
        ----------
        name : `str`
            The name of the fluid\n
        gamma : `float`
            The ratio of specific heats \n
        R : `flaot`
            The gas constant for the fluid \n
        units : `str`
            The units being used for the gas constant. Default is metric \n

        """
        self.name = name
        self.gamma = gamma
        self.R = R
        self.units = units



#Initialize some fluids in the metric system
#==================================================
air = fluid(name='Air', gamma=1.4, R=286.9, units='J / kg-K')

argon = fluid(name='Argon', gamma=1.67, R=208, units='J / kg-K')

CO2 = fluid(name='Carbon Dioxide', gamma=1.29, R=189, units='J / kg-K')

CO = fluid(name='Carbon Monoxide', gamma=1.4, R=297, units='J / kg-K')

hydrogen = fluid(name='Hydrogen', gamma=1.41, R=4120, units='J / kg-K')

helium = fluid(name='Helium', gamma=1.67, R=2080, units='J / kg-K')

methane = fluid(name='Methane', gamma=1.32, R=519, units='J / kg-K')

nitrogen = fluid(name='Nitrogen', gamma=1.4, R=296, units='J / kg-K')

O2 = fluid(name='Oxygen', gamma=1.4, R=260, units='J / kg-K')

water = fluid(name='water', gamma=1.33, R=461, units='J / kg-K')


#and in the british standard (Btu / lbm-R)
#============================================================
air_us = fluid(name='Air', gamma=1.4, R=53.3, units='Btu / lbm-R')

argon_us = fluid(name='Argon', gamma=1.67, R=38.7, units='Btu / lbm-R')

CO2_us = fluid(name='Carbon Dioxide', gamma=1.29, R=35.1, units='Btu / lbm-R')

CO_us = fluid(name='Carbon Monoxide', gamma=1.4, R=55.2, units='Btu / lbm-R')

hydrogen_us = fluid(name='Hydrogen', gamma=1.41, R=766, units='Btu / lbm-R')

helium_us = fluid(name='Helium', gamma=1.67, R=386, units='Btu / lbm-R')

methane_us = fluid(name='Methane', gamma=1.32, R=96.4, units='Btu / lbm-R')

nitrogen_us = fluid(name='Nitrogen', gamma=1.4, R=55.1, units='Btu / lbm-R')

O2_us = fluid(name='Oxygen', gamma=1.4, R=48.3, units='Btu / lbm-R')

water_us = fluid(name='water', gamma=1.33, R=85.7, units='Btu / lbm-R')



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