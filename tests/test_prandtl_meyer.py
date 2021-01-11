##############################
# Test prandtl_meyer functions
##############################
import gas_dynamics as gd
from gas_dynamics.fluids import air, methane
import random

class Test_prandtl_meyer_angle_from_mach:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.prandtl_meyer_angle_from_mach(a))

    def test_two(self):
        a = random.uniform(1.01,10)
        assert float(gd.prandtl_meyer_angle_from_mach(a, methane))

    def test_three(self):
        assert gd.prandtl_meyer_angle_from_mach(1) == 0.0

    def test_four(self):
        assert gd.prandtl_meyer_angle_from_mach(1, methane) == 0.0
        
    def test_five(self):
        zero = gd.prandtl_meyer_angle_from_mach(2) - 26.37976
        assert abs(zero) < 1e-4


class Test_prandtl_meyer_mach_from_angle:
    def test_one(self):
        a = random.uniform(0,130)
        assert float(gd.prandtl_meyer_mach_from_angle(a))

    def test_two(self):
        a = random.uniform(0,130)
        assert float(gd.prandtl_meyer_mach_from_angle(a, methane))

    def test_three(self):
        zero = gd.prandtl_meyer_mach_from_angle(0) - 1.0
        assert abs(zero) < 1e-5

    def test_four(self):
        zero = gd.prandtl_meyer_mach_from_angle(0, methane) - 1.0
        assert abs(zero) < 1e-5

    def test_five(self):
        zero = gd.prandtl_meyer_mach_from_angle(45) - 2.764452 
        assert abs(zero) < 1e-5


class Test_mach_wave_angle:
    def test_one(self):
        a = random.uniform(1,10)
        assert float(gd.mach_wave_angle(a))

    def test_two(self):
        zero = gd.mach_wave_angle(2.0) - 30
        assert abs(zero) < 1e-5