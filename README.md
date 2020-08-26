# Gas_Dynamics
Equations from the Oscar Biblarz Gas Dynamics book turned into python funcs with some explanations for solidifying my understanding. just for fun.

Getting Isentropic flow tables for a given gamma:
Arguments are min and max mach numbers, increment, and a list of gammas. Default valules are .01, 5, .1, and 1.4.

```
gd.print_stgn_ratios(gamma = [1.2]);
Isentropic Flow Parameters for γ = 1.2
M: 0.10   |   P/Pt: 0.994   |    T/Tt: 0.999   |    A/A*: 5.953   |   rho/rho_t: 0.995
M: 0.20   |   P/Pt: 0.976   |    T/Tt: 0.996   |    A/A*: 3.026   |   rho/rho_t: 0.980
M: 0.30   |   P/Pt: 0.948   |    T/Tt: 0.991   |    A/A*: 2.073   |   rho/rho_t: 0.956
M: 0.40   |   P/Pt: 0.909   |    T/Tt: 0.984   |    A/A*: 1.615   |   rho/rho_t: 0.924
M: 0.50   |   P/Pt: 0.862   |    T/Tt: 0.976   |    A/A*: 1.356   |   rho/rho_t: 0.884
M: 0.60   |   P/Pt: 0.809   |    T/Tt: 0.965   |    A/A*: 1.199   |   rho/rho_t: 0.838
M: 0.70   |   P/Pt: 0.750   |    T/Tt: 0.953   |    A/A*: 1.100   |   rho/rho_t: 0.787
M: 0.80   |   P/Pt: 0.689   |    T/Tt: 0.940   |    A/A*: 1.041   |   rho/rho_t: 0.733
M: 0.90   |   P/Pt: 0.627   |    T/Tt: 0.925   |    A/A*: 1.010   |   rho/rho_t: 0.677
M: 1.00   |   P/Pt: 0.564   |    T/Tt: 0.909   |    A/A*: 1.000   |   rho/rho_t: 0.621
M: 1.10   |   P/Pt: 0.504   |    T/Tt: 0.892   |    A/A*: 1.009   |   rho/rho_t: 0.565
M: 1.20   |   P/Pt: 0.446   |    T/Tt: 0.874   |    A/A*: 1.034   |   rho/rho_t: 0.510
M: 1.30   |   P/Pt: 0.392   |    T/Tt: 0.855   |    A/A*: 1.075   |   rho/rho_t: 0.458
M: 1.40   |   P/Pt: 0.342   |    T/Tt: 0.836   |    A/A*: 1.132   |   rho/rho_t: 0.409
M: 1.50   |   P/Pt: 0.296   |    T/Tt: 0.816   |    A/A*: 1.205   |   rho/rho_t: 0.363
M: 1.60   |   P/Pt: 0.255   |    T/Tt: 0.796   |    A/A*: 1.296   |   rho/rho_t: 0.320
M: 1.70   |   P/Pt: 0.218   |    T/Tt: 0.776   |    A/A*: 1.407   |   rho/rho_t: 0.281
M: 1.80   |   P/Pt: 0.186   |    T/Tt: 0.755   |    A/A*: 1.540   |   rho/rho_t: 0.246
M: 1.90   |   P/Pt: 0.157   |    T/Tt: 0.735   |    A/A*: 1.697   |   rho/rho_t: 0.214
M: 2.00   |   P/Pt: 0.133   |    T/Tt: 0.714   |    A/A*: 1.884   |   rho/rho_t: 0.186
M: 2.10   |   P/Pt: 0.112   |    T/Tt: 0.694   |    A/A*: 2.103   |   rho/rho_t: 0.161
M: 2.20   |   P/Pt: 0.094   |    T/Tt: 0.674   |    A/A*: 2.359   |   rho/rho_t: 0.139
M: 2.30   |   P/Pt: 0.078   |    T/Tt: 0.654   |    A/A*: 2.660   |   rho/rho_t: 0.120
M: 2.40   |   P/Pt: 0.065   |    T/Tt: 0.635   |    A/A*: 3.011   |   rho/rho_t: 0.103
M: 2.50   |   P/Pt: 0.054   |    T/Tt: 0.615   |    A/A*: 3.421   |   rho/rho_t: 0.088
M: 2.60   |   P/Pt: 0.045   |    T/Tt: 0.597   |    A/A*: 3.898   |   rho/rho_t: 0.076
M: 2.70   |   P/Pt: 0.037   |    T/Tt: 0.578   |    A/A*: 4.455   |   rho/rho_t: 0.065
M: 2.80   |   P/Pt: 0.031   |    T/Tt: 0.561   |    A/A*: 5.103   |   rho/rho_t: 0.055
M: 2.90   |   P/Pt: 0.026   |    T/Tt: 0.543   |    A/A*: 5.858   |   rho/rho_t: 0.047
M: 3.00   |   P/Pt: 0.021   |    T/Tt: 0.526   |    A/A*: 6.735   |   rho/rho_t: 0.040
M: 3.10   |   P/Pt: 0.018   |    T/Tt: 0.510   |    A/A*: 7.755   |   rho/rho_t: 0.034
M: 3.20   |   P/Pt: 0.015   |    T/Tt: 0.494   |    A/A*: 8.940   |   rho/rho_t: 0.029
M: 3.30   |   P/Pt: 0.012   |    T/Tt: 0.479   |    A/A*: 10.315   |   rho/rho_t: 0.025
M: 3.40   |   P/Pt: 0.010   |    T/Tt: 0.464   |    A/A*: 11.910   |   rho/rho_t: 0.021
M: 3.50   |   P/Pt: 0.008   |    T/Tt: 0.449   |    A/A*: 13.759   |   rho/rho_t: 0.018
M: 3.60   |   P/Pt: 0.007   |    T/Tt: 0.436   |    A/A*: 15.899   |   rho/rho_t: 0.016
M: 3.70   |   P/Pt: 0.006   |    T/Tt: 0.422   |    A/A*: 18.376   |   rho/rho_t: 0.013
M: 3.80   |   P/Pt: 0.005   |    T/Tt: 0.409   |    A/A*: 21.238   |   rho/rho_t: 0.011
M: 3.90   |   P/Pt: 0.004   |    T/Tt: 0.397   |    A/A*: 24.543   |   rho/rho_t: 0.010
M: 4.00   |   P/Pt: 0.003   |    T/Tt: 0.385   |    A/A*: 28.355   |   rho/rho_t: 0.008
M: 4.10   |   P/Pt: 0.003   |    T/Tt: 0.373   |    A/A*: 32.748   |   rho/rho_t: 0.007
M: 4.20   |   P/Pt: 0.002   |    T/Tt: 0.362   |    A/A*: 37.805   |   rho/rho_t: 0.006
M: 4.30   |   P/Pt: 0.002   |    T/Tt: 0.351   |    A/A*: 43.619   |   rho/rho_t: 0.005
M: 4.40   |   P/Pt: 0.002   |    T/Tt: 0.341   |    A/A*: 50.297   |   rho/rho_t: 0.005
M: 4.50   |   P/Pt: 0.001   |    T/Tt: 0.331   |    A/A*: 57.959   |   rho/rho_t: 0.004
M: 4.60   |   P/Pt: 0.001   |    T/Tt: 0.321   |    A/A*: 66.737   |   rho/rho_t: 0.003
M: 4.70   |   P/Pt: 0.001   |    T/Tt: 0.312   |    A/A*: 76.785   |   rho/rho_t: 0.003
M: 4.80   |   P/Pt: 0.001   |    T/Tt: 0.303   |    A/A*: 88.271   |   rho/rho_t: 0.003
M: 4.90   |   P/Pt: 0.001   |    T/Tt: 0.294   |    A/A*: 101.387   |   rho/rho_t: 0.002
M: 5.00   |   P/Pt: 0.001   |    T/Tt: 0.286   |    A/A*: 116.344   |   rho/rho_t: 0.002
```
![Stagnation_relations](https://github.com/fernancode/gas_dynamics/blob/master/print_ratios.png)

Plotting Stagnation relations versus mach number for different gammas. Arguments are min and max mach numbers, increment, and a list of gammas. Default valules are .01, 5, .1, and 1.4.

```
plot_stgn_ratios(gamma = [1.2,1.4,1.6,1.8])
```
![Stagnation_plots](https://github.com/fernancode/gas_dynamics/blob/master/plot_ratios.png)
