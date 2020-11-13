#!usr/bin/env
###
#The fluid class and some common fluids and their properties in metric and standard
###

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
        self.cp = None
        self.cv = None
        self.R = R
        self.gc = 1 #the gravitational constant. default 1 for metric
        self.units = units
        self.rho = None
        self.temperature = None 
        self.velocity = None 
        self.a = None #speed of sound
        self.mach = None
        self.mass_veolcity = None #rho * velocity
        self.a = None #speed of sound

    def set_a(self):
        """Set the local speed of sound for the fluid given its state

        """
        
        self.a = (self.temperature * self.R * self.gamma * self.gc)**.5

    def set_mach(self):
        """Set the mach number for the fluid given its velocity and speed of sound

        """
        self.set_a()
        self.mach = self.velocity / self.a
        
    def set_mass_velocity(self):
        """Set the mass velocity for the fluid given its state

        """

        self.mass_velocity = self.velocity * self.rho

    

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