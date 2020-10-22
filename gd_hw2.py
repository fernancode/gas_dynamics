import gas_dynamics as gd

p_t = 20e5 / 1000               #stag p inlet
p_rec = 200
T_t = 40 + 273.15       #stag t inlet
a_throat = 10 * .0001   #cm^2 to meters
#we want 200kPa (p3) to equal p3/pt3 * pt3/pt1 * pt1
#from rearranging and getting pt3 we can get M2, knowing M2 we can get the area ratio for the third critical

M2 = gd.stgn_pressure(p_t=p_t, p=p_rec, get='M')
A2_A1 = gd.area_mach_ratio(M1=1, M2=M2)
A2 = a_throat * A2_A1

#for a shock to exist at the exit plane, the pressure just after the plane needs to be high enough such that a compression process occurs and causes a normal shock
#this presure ratio looks like p4 = p4/p3 * p3/pt3 * pt3/pt1 *pt1/p1
#can also be obtained diretly from the tables.
p2_p1_shock = gd.shock_pressure_ratio(M=M2)
p_for_shock = p2_p1_shock * p_rec

#receiver pressure for shock at throat would be such that we are at the first critical
#p2/p1 = 1, therefore p2 = p2/pt2 * pt2/pt1 * pt1/p1 * p2/p1
p2_p1 = gd.p_stgn_ratio(M=1)
p2 = p2_p1 * p_t


print('1a) Exit area is %0.4f m^2' % A2)
print('1b) Pressure for shock at exit plane is %0.2f kPa' % p_for_shock)
print('1c) Pressure for shock at throat is %0.2f kPa' % p2)
print('\n')

p_t = 1000      #kPa
T_t = 300       #kelvin
p2 = 280       #kpa prior to shock
M4 = .5
mdot = 5
a_star = mdot / gd.mdot_a_star(p_t=p_t*1000, T_t=T_t)
#use the pressure ratio to get the mach number and the area reqiured for that accelerated flow
M2 = gd.stgn_pressure(p=p2, p_t=p_t, get='M')
A2 = a_star * gd.area_mach_ratio(M1=1, M2=M2)
#get the pressure ratio required for the shock and get p_rec
prec_p3 = gd.shock_pressure_ratio(M=M2)
prec = prec_p3 * p2/p_t * p_t / 1000
#the mach number after the shock is a function of Mach only, use this to get A2/A1 for M2=.5
M3 = gd.shock_mach(M1=M2)
A4 = A2 * gd.area_mach_ratio(M1=M3, M2=M4)
M4 = 1.795
#method 1, assume mach of .5 immediately following shock 
#M1 = gd.shock_mach(M2=M2)
#A2 = a_star * gd.area_mach_ratio(M1=1, M2=M1)
#prec_p3 = gd.shock_pressure_ratio(M1=M1)
#prec = prec_p3 * p3/p_t

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
#from the oblique shock chart grab dirac for theta and M1 of 2.5
dirac =  5 * gd.np.pi / 180
M2 = M2normal/gd.np.sin(theta-dirac)
T2 = T1 * gd.shock_temperature_ratio(M1normal)
dirac_degrees = dirac * 180 / gd.np.pi


print('3a) Resulting Mach number is %0.2f, ' %M2, 'Temperature is %0.2f K,' % T2 ,'and deflection angle is %0.0f degrees' %dirac_degrees)
print('3b) The flow deflection is 5 degrees to divert the flow back to the normal')

#TODO: this is wrong, fix
print('3c) By consulting the isentropc flow tables and seeing that a turn angle of 32.25 degrees could have resulted in a mach number pf 2.22 in region 2, we see that that for a turn of 32.25 + 5 degrees (to return the flow normal to the exit plane) this would result in a mach number of 2.43 for region 3. The resulting pressure is 0.04 bar and the resulting temperature is 68.7 Kelvin')
print('\n')


#problem 4
