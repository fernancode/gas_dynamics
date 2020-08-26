import gas_dynamics as gd

#gd.plot_stgn_ratios(.1,5,.1,gamma=[1.2,1.4,1.6,1.8])

gd.print_stgn_ratios(Mach_min=.1,Mach_max = 5,increment = .1,gamma=[1.2])