import gas_dynamics as gd

#stagnation temps, pressures, gas properties 
P = 1.8616e7
T = 293.15
R = 296
gamma = 1.4

#my throat area
area_in = 1.571**2 * gd.np.pi / 4
area_m = area_in * 0.00064516

#mdot per choked unit area
mdot_astar = gd.choked_mdot(pt = P, Tt = T, R = R, gamma = gamma)
#meters squared to inches squared
mdot_astar_in2 = mdot_astar / 1550 



#need as much mdot as possible to choke the flow
known_mdot = 1
area2choke = known_mdot/mdot_astar

#given a flow rate and the stagnation conditions what area is necessary for choked flow
area2choke_in = area2choke / .00064516
dia_in = (area2choke_in*4/gd.np.pi)**.5

print( dia_in)