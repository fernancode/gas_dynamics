from unittest import TestCase
import gas_dynamics as gd
from gas_dynamics.extra import methane


#==================================================
#Build a whole test suite for all funcs, default arguments, and other gasses
#==================================================


#==================================================
#sonic_velocity test
#==================================================
class test_sonic_velocity(TestCase):
    #test we get a number
    def test_is_float(self):
        a = gd.sonic_velocity()
        self.assertTrue(isinstance(a, float))


#==================================================
#stagnation pressure test
#==================================================
class test_stagnation_pressure(TestCase):
    #test we get a number
    def test1(self):
        p = gd.stagnation_pressure(pt=2, M=1)
        self.assertTrue(isinstance(p, float))
        
        pt = gd.stagnation_pressure(p=p, M=1)
        self.assertTrue(isinstance(pt, float))

        M = gd.stagnation_pressure(pt=pt, p=p)
        self.assertTrue(isinstance(M, float))

        self.assertNotEqual(p,pt)
        self.assertNotEqual(p,M)
        self.assertNotEqual(pt,M)        
        

#==================================================
#stagnation temperature test
#==================================================
class test_stagnation_temperature(TestCase):
    #test we get a number
    def test1(self):
        T = gd.stagnation_temperature(Tt=2, M=1)
        self.assertTrue(isinstance(p, float))
        
        Tt = gd.stagnation_temperature(T=T, M=1)
        self.assertTrue(isinstance(Tt, float))

        M = gd.stagnation_pressure(Tt=Tt, T=T)
        self.assertTrue(isinstance(M, float))

        self.assertNotEqual(T,Tt)
        self.assertNotEqual(T,M)
        self.assertNotEqual(Tt,M)        
        

#==================================================
#stagnation pressure ratio
#==================================================
class test_stagnation_pressure_ratio(TestCase):
    #test we get a number


class test_shock_mach(TestCase):
    def test_is_float(self):
        M2 = gd.shock_mach(M1=1)
        self.assertTrue(isinstance(M2, float))
