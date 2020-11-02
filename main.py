import gas_dynamics as gd

gd.plot_stagnation_ratios()
#gd.shock_oblique_charts(Mach_max=5, gas='Argon',lite=False)

#add a bunch of crap
M2 = gd.shock_mach(M1=10, gas='argon')