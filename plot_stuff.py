import matplotlib
import gas_dynamics as gd
import numpy as np

#m_nums = [i for i in np.arange(.01,10,.01)]

t_list = gd.Temp_stgn_ratio(.9)#for i in m_nums]
#p_list = [gd.Pressure_stgn_ratio(i) for i in m_nums]
#a_list = [gd.Area_stgn_ratio(i) for i in m_nums]
#rho_list = [gd.Density_stgn_ratio(i) for i in m_nums]

print(t_list)