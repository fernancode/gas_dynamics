import gas_dynamics as gd

#gd.plot_stgn_ratios(gamma = [1.2,1.4,1.6,1.8])

#gd.print_stgn_ratios(Mach_min=.1,Mach_max = 5,increment = .05,gamma=[1.33])

gd.shock_tables()

gd.shock_oblique_charts()

gd.shock_stagnation_ratio()
