import gas_dynamics as gd

#problem 1a
a = gd.sonic_velocity(T=323.15)
M = gd.stgn_pressure(p=2, p_t=2.2, get = 'M')
V = a*M
print('1a) Compressible velocity is ',V)

#problem 1b
V = ( (2.2-2) * 2/2 * 286 * 323.15)**.5
print('1b Incompressible velocity is ',V)

#problem 1c
M = gd.stgn_pressure(p=2, p_t=5.5, get = 'M')
Vc = a*M
Vi = ( (5.5-2) * 2/2 * 286 * 323.15)**.5
print('1c) Compressible velocity is ,',Vc,', incompressible velocity is ',Vi)


#problem 2a
gamma = 1.32
p1 = 14
T1 = 500
V1 = 125
R = 518.3
M2 = 0.8
M1 = V1 / gd.sonic_velocity(gamma=gamma, R = R, T = T1)
T2 = gd.temperature_mach_ratio(T1=T1,M1=M1,M2=.8,get='T2',gamma=gamma)
a = gd.sonic_velocity(gamma=gamma,R = R, T = T2)
V2 = M2 * a
print('2a) Temperature 2 is', T2, 'K and velocity 2 is ', V2,'meters/second.')

#problem 2b
p2 = gd.pressure_mach_ratio(p1=p1, M1 = M1, M2 = M2, gamma=gamma, R = R)
print('2b) Pressure 2 assuming no friction losses is ', p2, ' bar')

#problem 2c
a_ratio = gd.area_mach_ratio(M1,M2,gamma=gamma, R = R)
print('2c) The area ratio is ', a_ratio)




#problem 3a
p1 = 7
T1 = 600
p2 = 4
T2 = 550
M2 = .9
gamma = 1.29
R = 184

M1 = gd.temperature_mach_ratio(T1=T2, T2=T1, M1=M2, get='M2',gamma=gamma)
V1 = M1 * gd.sonic_velocity(gamma=gamma,R=R,T=T1)
print('3a) The velocity is ', V1)

#problem 3b
ds = gd.pressure_mach_ratio(p1=p1, p2=p2, M1=M1, M2=M2, gamma=gamma, R=R, get='ds')
print(ds)

#problem 3c
a_ratio = gd.area_mach_ratio(M1,M2,gamma=gamma,R=R)
print('3c) Area ratio is ', a_ratio)


#problem 4a,b
M1 = 0.3
T1 = 450
p1 = 10
A1 = 0.1

a_ratio = gd.a_star_ratio(M=M1)
a_star = A1/a_ratio

p_t = gd.stgn_pressure(p=p1,M=M1)
T_t = gd.stgn_temperature(T=T1, M=M1)

m_dot = a_star * gd.mdot_a_star(p_t=p_t, T_t=T_t)
print('4a) Flow rate is ',m_dot,' kg/s')
print('4b) The throat area is ',a_star,'m^2')

#problem 4c
p_rec = 1.5
M2 = gd.pressure_mach_ratio(p1=p_t, p2=p_rec, M1=0, get='M2')
print(M2)
a2 = a_star * gd.a_star_ratio(M2)
print('4c) The exit area is ',a2,'m^2')