import gas_dynamics as gd

#passed
gd.stagnation_ratios()
gd.stagnation_ratios(range=[1,2], inc=.1, gamma=1.5)

#passed
print(gd.sonic_velocity())
print(gd.sonic_velocity(gamma=1.5, R=288, T=300))
print('\n')

#need to check and fix docstrings
#gd.mach_pressure_ratio()
#gd.mach_temperature_ratio()
#gd.mach_area_ratio()

#passed
print(gd.mach_area_choked_ratio(M=1))
print('\n')

#passed
pt = gd.stagnation_pressure(p=100, M=1, gamma=1.5)
p = gd.stagnation_pressure(pt=pt, M=1, gamma=1.5)
M = gd.stagnation_pressure(p=p, pt=pt, gamma=1.5)
print(pt,p,M)
print('\n')

#passed
print(gd.stagnation_pressure_ratio(1))
print(gd.stagnation_pressure_ratio(1, gamma=1.5))
print('\n')

#passed
Tt = gd.stagnation_temperature(T=273, M=1, gamma=1.5)
T = gd.stagnation_temperature(Tt=Tt, M=1, gamma=1.5)
M = gd.stagnation_temperature(Tt=Tt, T=T , gamma=1.5)
print(Tt,T,M)
print('\n')

#passed
print(gd.stagnation_temperature_ratio(1))
print(gd.stagnation_temperature_ratio(1, gamma=1.5))
print('\n')

#passed
print(gd.stagnation_density_ratio(1))
print(gd.stagnation_density_ratio(1, gamma=1.5))
print('\n')

#passed
print(gd.choked_mdot(pt=100000, Tt=273))

#passed
gd.shock_tables()
gd.shock_tables(range=[1,2], inc=.1, gamma=1.5)

#passed
gd.shock_oblique_charts(gamma=1.5)

gd.plot_stagnation_ratios()

