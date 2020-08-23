module GasDynamics

"""
Takes density, velocity, mach number, change in area, area to calculate change in pressre for an arbitrary fluid
"""
function dp_from_mach(rho,V,M,da,A)
    dp = rho * V^2 * (1/(1-M^2)) * (dA/A)
end



""" 
Takes density, mach number, change in area, and area to calculate change in density for an arbitrary fluid
"""

function drho_from_mach(rho,M,dA,A)
    drho = rho * M^2/(1-M^2) * dA/A
end



"""
Takes velocity, mach number, change in area, and area to calculate dv for an arbitrary fluid
"""
function dv_from_mach(V,M,dA,A)
    dV = -V * (1/(1-M^2)) * (dA/A)
end



"""
Given Mach number and gamma, returns the relation of T / Tt
"""
function T_stgn(M,gamma = 1.4)
    Tt_ratio = 1 / (1+(gamma-1)/2 *M^2)
end



"""
given Mach number and gamma, returns the relation of P / Pt
"""
function P_stgn(M,gamma = 1.4)
    denom = 1 + (gamma-1)/2 * M^2
    Pt_ratio = (1 / denom ) ^ (gamma/(gamma-1))
end



"""
given Mach number and gamma, returns the relation of A / A*
"""
function Area_stgn(M,gamma = 1.4)
    A_star_ratio = 1/M * ((1 + (gamma-1)/2 * M^2) / ((gamma+1)/2)) ^ ((gamma+1)/(2*(gamma-1)))
end



"""
given Mach number and gamma, returns the relation of rho / rho_t
"""
function rho_stgn(M,gamma = 1.4)
    rho_t_ratio = (1 / (1 + (gamma-1)/2 * M^2 )) ^ (1 / (gamma-1))
end


#end statement for module
end



