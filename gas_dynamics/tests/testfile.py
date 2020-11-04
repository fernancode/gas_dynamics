from unittest import TestCase

import gas_dynamics as gd

class test_shock_mach(TestCase):
    def test_is_float(self):
        M2 = gd.shock_mach(M1=1)
        self.assertTrue(isinstance(M2, float))