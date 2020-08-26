# Gas_Dynamics
Equations from the Oscar Biblarz Gas Dynamics book turned into python funcs with some explanations for solidifying my understanding. just for fun.

Getting Isentropic flow tables for a given gamma:
Arguments are min and max mach numbers, increment, and a list of gammas. Default valules are .01, 5, .1, and 1.4.

```
>>> import gas_dynamics as gd
>>> gd.print_stgn_ratios(increment = .2, gamma = [1.2]);
Isentropic Flow Parameters for γ = 1.2
M: 0.00   |   P/Pt: 1.000   |    T/Tt: 1.000   |    A/A*: inf   |   rho/rho_t: 1.000
M: 0.20   |   P/Pt: 0.976   |    T/Tt: 0.996   |    A/A*: 3.026   |   rho/rho_t: 0.980
M: 0.40   |   P/Pt: 0.909   |    T/Tt: 0.984   |    A/A*: 1.615   |   rho/rho_t: 0.924
M: 0.60   |   P/Pt: 0.809   |    T/Tt: 0.965   |    A/A*: 1.199   |   rho/rho_t: 0.838
M: 0.80   |   P/Pt: 0.689   |    T/Tt: 0.940   |    A/A*: 1.041   |   rho/rho_t: 0.733
M: 1.00   |   P/Pt: 0.564   |    T/Tt: 0.909   |    A/A*: 1.000   |   rho/rho_t: 0.621
M: 1.20   |   P/Pt: 0.446   |    T/Tt: 0.874   |    A/A*: 1.034   |   rho/rho_t: 0.510
M: 1.40   |   P/Pt: 0.342   |    T/Tt: 0.836   |    A/A*: 1.132   |   rho/rho_t: 0.409
M: 1.60   |   P/Pt: 0.255   |    T/Tt: 0.796   |    A/A*: 1.296   |   rho/rho_t: 0.320
M: 1.80   |   P/Pt: 0.186   |    T/Tt: 0.755   |    A/A*: 1.540   |   rho/rho_t: 0.246
M: 2.00   |   P/Pt: 0.133   |    T/Tt: 0.714   |    A/A*: 1.884   |   rho/rho_t: 0.186
M: 2.20   |   P/Pt: 0.094   |    T/Tt: 0.674   |    A/A*: 2.359   |   rho/rho_t: 0.139
M: 2.40   |   P/Pt: 0.065   |    T/Tt: 0.635   |    A/A*: 3.011   |   rho/rho_t: 0.103
M: 2.60   |   P/Pt: 0.045   |    T/Tt: 0.597   |    A/A*: 3.898   |   rho/rho_t: 0.076
M: 2.80   |   P/Pt: 0.031   |    T/Tt: 0.561   |    A/A*: 5.103   |   rho/rho_t: 0.055
M: 3.00   |   P/Pt: 0.021   |    T/Tt: 0.526   |    A/A*: 6.735   |   rho/rho_t: 0.040
M: 3.20   |   P/Pt: 0.015   |    T/Tt: 0.494   |    A/A*: 8.940   |   rho/rho_t: 0.029
M: 3.40   |   P/Pt: 0.010   |    T/Tt: 0.464   |    A/A*: 11.910   |   rho/rho_t: 0.021
M: 3.60   |   P/Pt: 0.007   |    T/Tt: 0.436   |    A/A*: 15.899   |   rho/rho_t: 0.016
M: 3.80   |   P/Pt: 0.005   |    T/Tt: 0.409   |    A/A*: 21.238   |   rho/rho_t: 0.011
M: 4.00   |   P/Pt: 0.003   |    T/Tt: 0.385   |    A/A*: 28.355   |   rho/rho_t: 0.008
M: 4.20   |   P/Pt: 0.002   |    T/Tt: 0.362   |    A/A*: 37.805   |   rho/rho_t: 0.006
M: 4.40   |   P/Pt: 0.002   |    T/Tt: 0.341   |    A/A*: 50.297   |   rho/rho_t: 0.005
M: 4.60   |   P/Pt: 0.001   |    T/Tt: 0.321   |    A/A*: 66.737   |   rho/rho_t: 0.003
M: 4.80   |   P/Pt: 0.001   |    T/Tt: 0.303   |    A/A*: 88.271   |   rho/rho_t: 0.003
M: 5.00   |   P/Pt: 0.001   |    T/Tt: 0.286   |    A/A*: 116.344   |   rho/rho_t: 0.002
```
![Stagnation_relations](https://github.com/fernancode/gas_dynamics/blob/master/print_ratios.png)

Plotting Stagnation relations versus mach number for different gammas. Arguments are min and max mach numbers, increment, and a list of gammas. Default valules are .01, 5, .1, and 1.4.

```
plot_stgn_ratios(gamma = [1.2,1.4,1.6,1.8])
```
![Stagnation_plots](https://github.com/fernancode/gas_dynamics/blob/master/plot_ratios.png)
