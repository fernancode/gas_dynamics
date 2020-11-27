######################
# Test shock functions
######################
import gas_dynamics as gd
from gas_dynamics.fluids import air, methane
import random

class Test_shock_mach:
    def test_one(self):
        assert gd.shock_mach(mach=1) == 1

    def test_two(self):
        a = random.uniform(1,10)
        assert float(gd.shock_mach(mach=a))

    def test_three(self):
        a = random.uniform(1,10)
        assert gd.shock_mach(mach=a) <= 1


class Test_shock_mach_before:
    def test_one(self):
        assert gd.shock_mach_before(mach=1) == 1

    def test_two(self):
        a = random.uniform(.38,1)
        print(a)
        assert float(gd.shock_mach_before(mach=a))

    def test_three(self):
        a = random.uniform(.38,1)
        assert gd.shock_mach_before(mach=a) >= 1


class Test_shock_pressure_ratio:
    def test_one(self):
        assert gd.shock_pressure_ratio(mach=1) == 1

    def test_two(self):
        a = random.uniform(1,10)
        assert gd.shock_pressure_ratio(mach=a) >= 1

    
class Test_shock_mach_from_pressure_ratio:
    def test_one(self):
        assert gd.shock_mach_from_pressure_ratio(pressure_ratio=1) == 1
    
    def test_two(self):
        a = random.uniform(1,10)
        assert gd.shock_mach_from_pressure_ratio(pressure_ratio=a) >= 1


class Test_shock_temperature_ratio:
    def test_one(self):
        assert gd.shock_temperature_ratio(mach=1) == 1
    
    def test_two(self):
        a = random.uniform(1,10)
        assert gd.shock_temperature_ratio(mach=a) >= 1


class Test_shock_dv_a:
    def test_one(self):
        assert gd.shock_dv_a(mach=1) == 0
    
    def test_two(self):
        a = random.uniform(1,10)
        assert gd.shock_dv_a(mach=a) >= 0


class Test_shock_stagnation_pressure_ratio:
    def test_one(self):
        assert gd.shock_stagnation_pressure_ratio(mach=1) == 1

    def test_two(self):
        a = random.uniform(1,10)
        assert gd.shock_stagnation_pressure_ratio(mach=1) <= 1


class Test_shock_flow_deflection:
    def test_one(self):
        a = random.uniform(1,10)
        assert abs(gd.shock_flow_deflection(mach=a, shock_angle=90)) < 1e-5

    def test_two(self):
        a = random.uniform(1,10)
        b = random.uniform(0,180)
        assert float(gd.shock_flow_deflection(mach=a, shock_angle=b))


class Test_shock_angle:
    def test_one(self):
        a = 3
        b = 25
        c = gd.shock_angle(mach=a, flow_deflection=b)
        assert isinstance(c, list)

    def test_two(self):
        a = 3
        b = 25
        c = gd.shock_angle(mach=a, flow_deflection=b)
        assert abs(c[0]-44.135928931) < 1e-5
        assert abs(c[1]-79.326212403) < 1e-5


class Test_shock_mach_given_angles:
    def test_one(self):
        a = 45
        b = 25
        c = gd.shock_mach_given_angles(shock_angle=44.1359289, flow_deflection=25)
        assert abs(c - 3.000000003555829) < 1e-5


class Test_shock_flow_deflection_from_machs:
    def test_one(self):
        assert gd.shock_flow_deflection_from_machs(mach_initial=1, mach_final=1) == 0

    def test_two(self):
        a = gd.shock_flow_deflection_from_machs(mach_initial=2, mach_final=1)
        assert a < 25 and a > 20