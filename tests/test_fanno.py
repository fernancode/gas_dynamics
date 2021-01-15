######################
# Test fanno functions
######################
import gas_dynamics as gd
from gas_dynamics.fluids import air, methane
import random

#TODO: these tests only test for float, not for actual correct values.
#could use more robust-ness and check versus tabulated values

class Test_fanno_temperature_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.fanno_temperature_ratio(a,b))

    def test_two(self):
        a = random.uniform(1.01,10)
        b = a #random.uniform(1.01,10)
        assert float(gd.fanno_temperature_ratio(a,b)) == 1


class Test_fanno_pressure_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.fanno_pressure_ratio(a,b))

    def test_two(self):
        a = random.uniform(1.01,10)
        b = a #random.uniform(1.01,10)
        assert float(gd.fanno_pressure_ratio(a,b)) == 1


class Test_fanno_density_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.fanno_density_ratio(a,b))

    def test_two(self):
        a = random.uniform(1.01,10)
        b = a #random.uniform(1.01,10)
        assert float(gd.fanno_density_ratio(a,b)) == 1


class Test_fanno_stagnation_pressure_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.fanno_stagnation_pressure_ratio(a,b))

    def test_two(self):
        a = random.uniform(1.01,10)
        b = a #random.uniform(1.01,10)
        assert float(gd.fanno_stagnation_pressure_ratio(a,b)) == 1


class Test_fanno_temperature_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.fanno_temperature_star_ratio(a))

    def test_two(self):
        assert float(gd.fanno_temperature_star_ratio(1)) == 1


class Test_fanno_pressure_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.fanno_pressure_star_ratio(a))

    def test_two(self):
        assert float(gd.fanno_pressure_star_ratio(1)) == 1


class Test_fanno_density_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.fanno_density_star_ratio(a))

    def test_two(self):
        assert float(gd.fanno_density_star_ratio(1)) == 1


class Test_fanno_velocity_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.fanno_velocity_star_ratio(a))

    def test_two(self):
        assert float(gd.fanno_velocity_star_ratio(1)) == 1


class Test_fanno_parameter:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.fanno_parameter(a,b))


class Test_fanno_parameter_max:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.fanno_parameter_max(a))

    def test_two(self):
        assert gd.fanno_parameter_max(1) == 0


#TODO: MAKE TEST
#class Test_mach_from_fanno:       