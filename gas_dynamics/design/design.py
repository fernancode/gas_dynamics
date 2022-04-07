from gas_dynamics.prandtl_meyer.prandtl_meyer import prandtl_meyer_mach_from_angle, prandtl_meyer_angle_from_mach, mach_wave_angle
from gas_dynamics.fluids import fluid, air
from gas_dynamics.extra import tand
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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

###Create c class to contain all of the data for a line
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

#Create a class to store all of the lines
class Nozzle:
    def __init__(self, num_lines):
        self.lines = []
        self.make_lines(num_lines)
        self.x = [0]
        self.y = [1]

    def make_lines(self, num_lines):
        for n in range(num_lines):
            self.lines.append(Characteristic_line(n))


#add a number of lines to the nozzle class
#TODO: Variable me
nozzle = Nozzle(100)

#determine a mach number
#TODO: Variable me
exit_mach = 3.19

#get the number of lines
num_lines = len(nozzle.lines)
theta_max = prandtl_meyer_angle_from_mach(mach=exit_mach) / 2

#determine all the turn angles for the lines
#TODO: Variable me
initial_turn = .1
theta_to_go = theta_max - initial_turn
theta_steps = theta_to_go / (num_lines - 1)
fig, ax = plt.subplots()

#initialize some values in the first line
[nozzle.lines[0].theta.append(initial_turn + (n)*theta_steps) for n in range(num_lines)]
[nozzle.lines[0].k_plus.append(0) for n in range(num_lines)] 
[nozzle.lines[0].k_minus.append(2*(initial_turn + (n)*theta_steps)) for n in range(num_lines)]
[nozzle.lines[0].nu.append(nu) for nu in nozzle.lines[0].theta]
[nozzle.lines[0].mach.append(prandtl_meyer_mach_from_angle(nu)) for nu in nozzle.lines[0].nu]
[nozzle.lines[0].mach_wave_angle.append(mach_wave_angle(mach)) for mach in nozzle.lines[0].mach]

#create values on the first line
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
        mach = prandtl_meyer_mach_from_angle(nu)
        line.mach.append(mach)
        wave_angle = mach_wave_angle(mach)
        line.mach_wave_angle.append(wave_angle)

#get points on first line
for n, theta in enumerate(nozzle.lines[0].theta):
    if n == 0:
        slope = tand(theta - nozzle.lines[n].mach_wave_angle[n])
        nozzle.lines[n].x.append( nozzle.y / -slope)
        nozzle.lines[n].y.append(0)
        continue

    slope_down = tand(nozzle.lines[0].theta[n] - nozzle.lines[0].mach_wave_angle[n])
    line_down = lambda x : slope_down*x + 1

    slope_up = tand(nozzle.lines[0].theta[n-1] + nozzle.lines[0].mach_wave_angle[n-1])
    b = nozzle.lines[0].y[n-1] -(slope_up * nozzle.lines[0].x[n-1])
    
    line_up = lambda x : slope_up*x + b
    zero = lambda x : line_down(x) - line_up(x)
    x = fsolve(zero,1)

    nozzle.lines[0].x.append(x[0])
    nozzle.lines[0].y.append(line_down(x[0]))

for n, line in enumerate(nozzle.lines):
    if n == 0:
        continue
    for m, theta in enumerate(line.theta):
        if m == 0:
            slope = tand(0 - line.mach_wave_angle[m])
            b = nozzle.lines[n-1].y[m+1] - slope * nozzle.lines[n-1].x[m+1]
            line.x.append(-b/slope)
            line.y.append(0)
            continue
        
        slope_down = tand(theta - line.mach_wave_angle[m])
        b_down = nozzle.lines[n-1].y[m] - (slope_down * nozzle.lines[n-1].x[m])
        line_down = lambda x : slope_down*x + b_down

        slope_up = tand(line.theta[m-1] + line.mach_wave_angle[m-1])
        b_up = line.y[m-1] - (slope_up * line.x[m-1])
        line_up = lambda x : slope_up*x + b_up

        zero = lambda x : line_down(x) - line_up(x)
        x = fsolve(zero,1)

        line.x.append(x[0])
        line.y.append(line_down(x[0]))
        print()

#finally go for the contour
#duplicate the angle at every final point
for line in nozzle.lines:
    line.theta.append(line.theta[-1])
    line.mach_wave_angle.append(line.mach_wave_angle[-1])


for n, line in enumerate(nozzle.lines[:-1]):
    if n == 0:
        slope = tand((line.theta[-1]))
        b_up = 1
        contour_line = lambda x : x*slope + b_up

        slope2 = tand(line.mach_wave_angle[-1])
        b_up2 = line.y[-1] - slope2*line.x[-1]
        line_up = lambda x : x*slope2 + b_up2

        zero = lambda x : contour_line(x) - line_up(x)
        x = fsolve(zero,1)
        line.x.append(x[0])
        line.y.append(contour_line(x[0]))
        nozzle.x.append(x[0])
        nozzle.y.append(contour_line(x[0]))
        continue

    slope = .5 * tand(line.theta[-1])
    b_up = nozzle.lines[n-1].y[-1] - slope*nozzle.lines[n-1].x[-1]
    contour_line = lambda x : x*slope + b_up

    slope2 = tand(line.mach_wave_angle[-1])
    b_up2 = line.y[-1] - slope2*line.x[-1]
    line_up = lambda x : x*slope2 + b_up2

    zero = lambda x : contour_line(x) - line_up(x)
    x = fsolve(zero,1)
    line.x.append(x[0])
    line.y.append(contour_line(x[0]))
    nozzle.x.append(x[0])
    nozzle.y.append(contour_line(x[0]))

for line in nozzle.lines:
    ax.plot(line.x, line.y, color='k', alpha=0.5)

ax.plot(nozzle.x, nozzle.y)
ax.grid(which='both', axis='both')
ax.axis('equal')
ax.set_xlabel(r'$Z$')
ax.set_ylabel(r'$X$', rotation=0)


plt.title(r'$Method \ of \ Characteristics \ Nozzle$')
plt.legend(title=r'$Characteristic \ Lines$')
print('To get an appropriate resolution, the initial turn angle needs to be around .01 degrees and 150 characteristic lines compared with the area ratio required')
plt.show()