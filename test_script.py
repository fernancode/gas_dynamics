import gas_dynamics as gd

#passed
gd.stagnation_ratios()
gd.stagnation_ratios(range=[1,2], inc=.1, gas='argon')

#passed
print(gd.sonic_velocity())
print(gd.sonic_velocity(gas='argon', T=300))
print('\n')

#need to check and fix docstrings
#gd.mach_pressure_ratio()
#gd.mach_temperature_ratio()
#gd.mach_area_ratio()

#passed
print(gd.mach_area_choked_ratio(M=1))
print('\n')

#passed
pt = gd.stagnation_pressure(p=100, M=1, gas='argon')
p = gd.stagnation_pressure(pt=pt, M=1, gas='argon')
M = gd.stagnation_pressure(p=p, pt=pt, gas='argon')
print(pt,p,M)
print('\n')

#passed
print(gd.stagnation_pressure_ratio(1))
print(gd.stagnation_pressure_ratio(1, gas='argon'))
print('\n')

#passed
Tt = gd.stagnation_temperature(T=273, M=1, gas='argon')
T = gd.stagnation_temperature(Tt=Tt, M=1, gas='argon')
M = gd.stagnation_temperature(Tt=Tt, T=T , gas='argon')
print(Tt,T,M)
print('\n')

#passed
print(gd.stagnation_temperature_ratio(1))
print(gd.stagnation_temperature_ratio(1, gas='argon'))
print('\n')

#passed
print(gd.stagnation_density_ratio(1))
print(gd.stagnation_density_ratio(1, gas='argon'))
print('\n')

#passed
print(gd.choked_mdot(pt=100000, Tt=273))

#passed
gd.shock_tables()
gd.shock_tables(range=[1,2], inc=.1, gas='argon')

#passed
gd.shock_oblique_charts(gas='arGon')

gd.plot_stagnation_ratios()

gd.sonic_velocity()