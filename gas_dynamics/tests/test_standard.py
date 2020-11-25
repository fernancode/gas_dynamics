########################
# Test standard function
########################
import gas_dynamics as gd
from gas_dynamics.fluids import air, methane

class Test_sonic_velocity:
    def test_one(self):
        assert float(gd.sonic_velocity(300, air))
        assert float(gd.sonic_velocity(300, methane))

    def test_two(self):
        assert gd.sonic_velocity(0) == 0


class Test_stagnation_pressure:
    def test_one(self):
        #around the world...
        pt = 100
        mach = 1
        p = gd.stagnation_pressure(stagnation_pressure=pt, mach=mach)
        assert abs(gd.stagnation_pressure(pressure=p, stagnation_pressure=pt) - mach) < 1e-5
        assert abs(gd.stagnation_pressure(mach=mach, pressure=p) - pt) < 1e-5


class Test_stagnation_temperature:
    def test_one(self):
        #around the world
        Tt = 100
        mach = 1
        T = gd.stagnation_temperature(stagnation_temperature=Tt, mach=mach, gas=methane)
        assert abs(gd.stagnation_temperature(stagnation_temperature=Tt, temperature=T, gas=methane) - mach) < 1e-5 
        assert abs(gd.stagnation_temperature(temperature=T, mach=mach, gas=methane) - Tt) < 1e-5


class Test_stagnation_pressure_ratio:
    def test_one(self):
        assert gd.stagnation_pressure_ratio(mach=0, gas=methane) == 1
        assert abs(gd.stagnation_pressure_ratio(mach=1) - 0.5282817) < 1e-5


class Test_stagnation_temperature_ratio:
    def test_one(self):
        assert gd.stagnation_temperature_ratio(mach=0, gas=methane) == 1
        assert abs(gd.stagnation_temperature_ratio(mach=1) - 0.8333333) < 1e-5


class Test_stagnation_density_ratio:
    def test_one(self):
        assert gd.stagnation_density_ratio(mach=0, gas=methane) == 1
        assert abs(gd.stagnation_density_ratio(mach=1) - 0.633938145) < 1e-5


#class Test_stagnation_ratio:
#    def test_one(self):
#        assert gd.stagnation_ratio() == str


#class Test_stagnation_ratio_table:


class Test_mach_from_pressure_ratio:
    def test_one(self):
        m1 = 1
        assert gd.mach_from_pressure_ratio(pressure_initial=1, pressure_final=1, mach_initial=1) == 1
        assert gd.mach_from_pressure_ratio(pressure_initial=2, pressure_final=1, mach_initial=1) > m1


class Test_mach_from_temperature_ratio:
    def test_one(self):
        m1 = 1
        assert gd.mach_from_temperature_ratio(temperature_initial=1, temperature_final=1, mach_initial=1) == 1
        assert gd.mach_from_temperature_ratio(temperature_initial=2, temperature_final=1, mach_initial=1) > m1


class Test_pressure_from_mach_ratio:
    def test_one(self):
        assert gd.pressure_from_mach_ratio(mach_initial=1, mach_final=1, pressure_initial=1) == 1
        assert gd.pressure_from_mach_ratio(mach_initial=1, mach_final=2, pressure_initial=1) < 1


class Test_temperature_from_mach_ratio:
    def test_one(self):
        assert gd.temperature_from_mach_ratio(mach_initial=1, mach_final=1, temperature_initial=1) == 1
        assert gd.temperature_from_mach_ratio(mach_initial=1, mach_final=2, temperature_initial=1) < 1


#Test_entropy_produced:
#    def test_one(self):


class Test_mach_area_star_ratio:
    def test_one(self):
        assert gd.mach_area_star_ratio(mach=1) == 1


class Test_mach_area_ratio:
    def test_one(self):
        assert gd.mach_area_ratio(mach_initial=1, mach_final=1) == 1


class Test_mach_from_area_ratio:
    def test_one(self):
        assert abs(gd.mach_from_area_ratio(area_ratio=1)[0] - 1) < 1e-5
        assert abs(gd.mach_from_area_ratio(area_ratio=1)[1] - 1) < 1e-5


class Test_mass_flux_funcs:
    def test_one(self):
        assert gd.mass_flux_max(stagnation_pressure=0, stagnation_temperature=1) == 0
        assert gd.mass_flux(mach=1, stagnation_pressure=0, stagnation_temperature=1) == 0
        from_max = gd.mass_flux_max(stagnation_pressure=100, stagnation_temperature=100)
        manual = gd.mass_flux(mach=1, stagnation_pressure=100, stagnation_temperature=100)
        assert abs(manual - from_max) < 1e-5