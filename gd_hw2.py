import gas_dynamics as gd

p_t = 20e5 / 1000               #stag p inlet
p_rec = 200
T_t = 40 + 273.15       #stag t inlet
a_throat = 10 * .0001   #cm^2 to meters
#we want 200kPa (p3) to equal p3/pt3 * pt3/pt1 * pt1
#from rearranging and getting pt3 we can get M2, knowing M2 we can get the area ratio for the third critical

M2 = gd.stagnation_pressure(pt=p_t, p=p_rec) #the mach number that would exist 
A2_A1 = gd.mach_area_ratio(M1=1, M2=M2)
A2 = a_throat * A2_A1

#for a shock to exist at the exit plane, the pressure just after the plane needs to be high enough such that a compression process occurs and causes a normal shock
#this presure ratio looks like p4 = p4/p3 * p3/pt3 * pt3/pt1 *pt1/p1
#can also be obtained diretly from the tables.
p2_p1_shock = gd.shock_pressure_ratio(M=M2)
p_for_shock = p2_p1_shock * p_rec

#receiver pressure for shock at throat would be such that we are at the first critical
#p2/p1 = 1, therefore p2 = p2/pt2 * pt2/pt1 * pt1/p1 * p2/p1
p2_p1 = gd.stagnation_pressure_ratio(M=1)
p2 = p2_p1 * p_t
#get stagnation values for throat and multiply by the pressure ratio for the expansion, check the tables for this value
pt2 = gd.stagnation_pressure(p=p2, M=1)
p_shock_throat = pt2 * 0.9315

print('1a) Exit area is %0.4f m^2' % A2)
print('1b) Pressure for shock at exit plane is %0.2f kPa' % p_for_shock)
print('1c) Pressure for shock at throat is %0.2f kPa' % p_shock_throat)
print('\n')

p_t = 1000      #kPa
T_t = 300       #kelvin
p2 = 280       #kpa prior to shock
M4 = .5
mdot = 5

a_star = mdot / gd.choked_mdot(pt=p_t*1000, Tt=T_t)
#use the pressure ratio to get the mach number and the area reqiured for that accelerated flow

M2 = gd.stagnation_pressure(p=p2, pt=p_t)
A2 = a_star * gd.mach_area_ratio(M1=1, M2=M2)

#get the stagnation pressure ratio across the shock and multiply by p_t to get pt4
pt4 = p_t * gd.shock_stagnation_ratio(M=M2)
#knowing pt4 and that M=.5, get p
prec = gd.stagnation_pressure(M=M4, pt=pt4)

#the mach number after the shock is a function of Mach only, use this to get A2/A1 for M2=.5
M3 = gd.shock_mach(M1=M2)
A4 = A2 * gd.mach_area_ratio(M1=M3, M2=M4)
M4 = 1.795

print('2a) Throat area: %0.5f m^2' % a_star)
print('2b) Area just prior to shock: %0.5f m^2' % A2)
print('2c) Receiver pressure is %0.3f bar' % prec)
print('2d) Final exit area is %0.5f m^2' % A4)
print('2e) Consulting the tables, the design mach number for the area ratio is %0.2f' % M4)
print('\n')

#Nitrogen!!!
M1 = 2.5    #initial mach
T1 = 150    #initial kelvin
p1 = 0.7    #bar
prec = 1   #receiver pressure
#an oblique shock forms to compensate the pressure difference. a pressure ratio forms of p_rec/p1
prec_p1 = prec/p1
#get the normal components across the shock
M1normal = gd.shock_pressure_ratio(p2_p1=prec_p1)
M2normal = gd.shock_mach(M1=M1normal)
theta = gd.np.arcsin(M1normal/M1)

dirac = gd.shock_flow_deflection(M=M1, theta=theta)
M2 = M2normal/gd.np.sin(theta-dirac)
T2 = T1 * gd.shock_temperature_ratio(M1normal)
dirac_degrees = dirac * 180 / gd.np.pi

#numerically solve for shock angle given dirac and mach
sol = gd.shock_angle(M=M2, dirac=dirac)
theta_2 = min(sol) #get the weak shock angle
M2n = M2*gd.np.sin(theta_2)
M3n = gd.shock_mach(M2n)
M3 = M3n / gd.np.sin(theta_2-dirac)
T3 = T2 * gd.shock_temperature_ratio(M2n)
p3 = prec * gd.shock_pressure_ratio(M2n)

print('3a) Resulting Mach number is %0.2f, ' %M2, 'Temperature is %0.2f K,' % T2 ,'and deflection angle is %0.2f degrees' %dirac_degrees)
print('3b) The flow deflection is %0.2f degrees to divert the flow back to the normal' %dirac_degrees)
print('3c) In region 3 the Mach number is %0.2f,' %M3, 'Temperature is %0.2f K,' %T3, 'Pressure is %0.2f.' %p3)
print('\n')

#problem 4
theta1= 40 * gd.np.pi/180
dirac = 15 * gd.np.pi/180
#get the mach 1 given the angles, get mach 2
M1 = gd.shock_mach_given_angles(theta=theta1, dirac=dirac)
M1n = M1 * gd.np.sin(theta1)
M2n = gd.shock_mach(M1n)
M2 = M2n / gd.np.sin(theta1-dirac)

#knowing mach 2 and flow deflection, get M3 and beta
shock_angles = gd.shock_angle(M=M2, dirac=dirac) 
theta2 = shock_angles[0]
theta2_deg = theta2 * 180/gd.np.pi
M2n_p = M2 * gd.np.sin(theta2)
M3n = gd.shock_mach(M2n_p)
#flow gets deflected back the other way so add dirac
M3 = M3n / gd.np.sin(theta2+dirac)

p3_p1 = gd.shock_pressure_ratio(M=M1n) * gd.shock_pressure_ratio(M=M2n_p)
pt3_pt1 = gd.shock_stagnation_ratio(M=M1n) * gd.shock_stagnation_ratio(M=M2n_p)
print('3a) Mach 3 is %0.2f' %M3,'and beta is %0.2f' %theta2_deg)
print('3b) The ratio p3/p1 is %0.3f' %p3_p1, 'and pt3/pt1 is %0.3f' %pt3_pt1)