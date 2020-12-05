from gas_dynamics.prandtl_meyer.prandtl_meyer import prandtl_meyer_mach, prandtl_meyer_turn, mach_wave_angle
from gas_dynamics.fluids import fluid, air
import numpy as np

def method_of_characteristics(exit_mach: float, characteristic_lines: 7, gas=air):
    """Generate a nozzle geometry using the method of characteristics

    Notes
    -----

    Parameters
    ----------
    exit_mach : float
        The design exit mach number
    characteristic_lines : int
        The number of characteristic lines used to generate the contour
    gas : fluid
        A user defined fluid object. Default is air

    Returns
    -------

    Examples
    --------

    """

    characteristic_lines
    if isinstance(characteristic_lines, int) == False:
        c_lines = round(characteristic_lines) + 1
    if characteristic_lines < 2 :
        characteristic_lines += 1

class Characteristic_line:
    def __init__(self, number):
        #this is K+ line, or "Right running line (bottom to top)"
        self.line = number          #the line number
        self.x = []                 #x coordinate of point
        self.y = []                 #y coordinate of point
        self.k_minus = []           #value of k-
        self.k_plus = []            #value of k+
        self.theta = []             #amount which the flow has turned
        self.nu = []                #the prandtl meyer turn angle corresponding to the mach number
        self.mach = []              #the mach number
        self.mach_wave_angle = []   #the mach wave angle


class Nozzle:
    def __init__(self, num_lines):
        self.lines = []
        self.make_lines(num_lines)

    def make_lines(self, num_lines):
        for n in range(num_lines):
            self.lines.append(Characteristic_line(n))


 #create a nozzle object to store all the information
#TODO: Variable me
nozzle = Nozzle(7)
#TODO: Variable me
exit_mach = 2.4

#get the number of lines
num_lines = len(nozzle.lines)

theta_max = prandtl_meyer_turn(mach=exit_mach) / 2

#TODO: Variable me
initial_turn = .375
theta_to_go = theta_max - initial_turn
theta_steps = theta_to_go / (num_lines - 1)

#initialize some values in the first line
[nozzle.lines[0].theta.append(initial_turn + (n)*theta_steps) for n in range(num_lines)]
[nozzle.lines[0].k_plus.append(0) for n in range(num_lines)] 
[nozzle.lines[0].k_minus.append(2*(initial_turn + (n)*theta_steps)) for n in range(num_lines)]
[nozzle.lines[0].nu.append(nu) for nu in nozzle.lines[0].theta]
[nozzle.lines[0].mach.append(prandtl_meyer_mach(nu)) for nu in nozzle.lines[0].nu]
[nozzle.lines[0].mach_wave_angle.append(mach_wave_angle(mach)) for mach in nozzle.lines[0].mach]

for n, line in enumerate(nozzle.lines):
    if n == 0:
        continue
    for m, k_minus in enumerate(nozzle.lines[n-1].k_minus[1:]):
        k_plus = -nozzle.lines[0].k_minus[n]

        line.k_plus.append(k_plus)
        line.k_minus.append(k_minus)

        if line.theta == []:
            line.theta.append(0)
        else:
            theta = .5 * (k_minus + k_plus)
            line.theta.append(theta)

        nu = .5 * (k_minus - k_plus)
        line.nu.append(nu)

        mach = prandtl_meyer_mach(nu)
        line.mach.append(mach)

        wave_angle = mach_wave_angle(mach)
        line.mach_wave_angle.append(wave_angle)

print()