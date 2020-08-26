# Gas_Dynamics
Equations from the Oscar Biblarz Gas Dynamics book turned into python funcs with some explanations for solidifying my understanding. just for fun.

Getting Isentropic flow tables for a given gamma:
Arguments are min and max mach numbers, increment, and a list of gammas. Default valules are .01, 5, .1, and 1.4.

![Stagnation_relations](https://github.com/fernancode/gas_dynamics/blob/master/print_ratios.png)

Plotting Stagnation relations versus mach number for different gammas. Arguments are min and max mach numbers, increment, and a list of gammas. Default valules are .01, 5, .1, and 1.4.
'''
plot_stgn_ratios(gamma = [1.2,1.4,1.6,1.8])
'''
![Stagnation_plots](https://github.com/fernancode/gas_dynamics/blob/master/plot_ratios.png)
