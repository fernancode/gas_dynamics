import gas_dynamics as gd

#problem 1
a = gd.sonic_velocity(T=323.15+50)
M = stgn_pressure(p=2, p_t=2.2, get = 'M')
V = a*M
print('velocity is ', V)
