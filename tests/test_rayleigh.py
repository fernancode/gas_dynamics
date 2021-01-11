#########################
# Test rayleigh functions
#########################
import gas_dynamics as gd
from gas_dynamics.fluids import air, methane
import random

#TODO: these tests only test for float, not for actual correct values.
#could use more robust-ness

class Test_rayleigh_pressure_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.rayleigh_pressure_ratio(a,b))


class Test_rayleigh_temperature_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.rayleigh_temperature_ratio(a,b))


class Test_rayleigh_density_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.rayleigh_density_ratio(a,b))


class Test_rayleigh_stagnation_temperature_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.rayleigh_stagnation_temperature_ratio(a,b))


class Test_rayleigh_stagnation_pressure_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        assert float(gd.rayleigh_stagnation_pressure_ratio(a,b))


class Test_rayleigh_mach_from_pressure_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        m = random.uniform(1.01,10)
        assert float(gd.rayleigh_mach_from_pressure_ratio(m,a,b))


class Test_rayleigh_mach_from_temperature_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        b = random.uniform(1.01,10)
        m = random.uniform(1.01,10)
        assert float(gd.rayleigh_mach_from_temperature_ratio(m,a,b))


#TODO: this isn't great
class Test_rayleigh_mach_from_stagnation_temperature_ratio:
    def test_one(self):
        a = random.uniform(1.01,5)
        b = a + 1
        m = random.uniform(1.01,5)
        assert float(gd.rayleigh_mach_from_stagnation_temperature_ratio(m,b,a))


class Test_rayleigh_mach_from_stagnation_pressure_ratio:
    def test_one(self):
        a = random.uniform(1.01,5)
        b = a + 1
        m = random.uniform(1.01,5)
        assert float(gd.rayleigh_mach_from_stagnation_pressure_ratio(m,b,a))


class Test_rayleigh_pressure_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.rayleigh_pressure_star_ratio(a))

    def test_two(self):
        assert float(gd.rayleigh_pressure_star_ratio(1)) == 1


class Test_rayleigh_temperature_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.rayleigh_temperature_star_ratio(a))

    def test_two(self):
        assert float(gd.rayleigh_temperature_star_ratio(1)) == 1


class Test_rayleigh_density_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.rayleigh_density_star_ratio(a))

    def test_two(self):
        assert float(gd.rayleigh_density_star_ratio(1)) == 1


class Test_rayleigh_stagnation_pressure_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.rayleigh_stagnation_pressure_star_ratio(a))

    def test_two(self):
        assert float(gd.rayleigh_stagnation_pressure_star_ratio(1)) == 1


class Test_rayleigh_stagnation_temperature_star_ratio:
    def test_one(self):
        a = random.uniform(1.01,10)
        assert float(gd.rayleigh_stagnation_temperature_star_ratio(a))

    def test_two(self):
        assert float(gd.rayleigh_stagnation_temperature_star_ratio(1)) == 1


class Test_rayleigh_heat_flux:
    def test_one(self):
        a = random.uniform(1.01,100)
        b = random.uniform(1.01,500)
        assert float(gd.rayleigh_heat_flux(a,b))