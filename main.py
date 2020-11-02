import gas_dynamics as gd
from scipy.optimize import fsolve

gd.shock_oblique_charts(Mach_max=5, gas='Argon',lite=False)

#gd.dirac_from_machs(M1=2,M2=1.5)

