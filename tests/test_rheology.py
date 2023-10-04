"""
Test for rheology.py
"""
import unittest

import gdmate as gd
import numpy as np

class TestHelloWorld(unittest.TestCase):

    # TODO: Test with one layer (and maybe 2/4+ layers) if applicable
    # TODO: Test on layers of varying thicknesses
    # Testing on a model with 3 layers of thickness 1 km with A = 0 and depth=3
    # Other values have been adjusted to make math cleaner
    def test_cond_geotherm_one_km_layers(self):
        bd_temps, bd_heat_flows, z, cond_temps = \
            gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=3,
             radiogenic_heat=[0., 0., 0.], surface_t=273,
             heat_flow=0.05, thermal_conductivity=2.5)
        self.assertTrue((bd_temps == [273., 293., 313., 333.]).all())
        self.assertTrue((bd_heat_flows == [0.05, 0.05, 0.05, 0.05]).all())
        self.assertTrue((z == [0, 1000, 2000, 3000]).all())
        self.assertTrue((cond_temps == [273., 293., 313., 333.]).all())
    
    # Testing on a model with 3 layers of thickness 1 km with A = 0 and depth=5
    # Other values have been adjusted to make math cleaner
    def test_cond_geotherm_one_km_layers_extra_depth(self):
        bd_temps, bd_heat_flows, z, cond_temps = \
            gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=5,
             radiogenic_heat=[0., 0., 0.], surface_t=273,
             heat_flow=0.05, thermal_conductivity=2.5)
        self.assertTrue((bd_temps == [273., 293., 313., 333.]).all())
        self.assertTrue((bd_heat_flows == [0.05, 0.05, 0.05, 0.05]).all())
        self.assertTrue((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
        self.assertTrue((cond_temps == 
                         [273., 293., 313., 333., 333., 333.]).all())
    
    # Testing on a model with 3 layers of thickness 1 km with A = 1e-5, depth=5
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    def test_cond_geotherm_one_km_layers_fixed_A(self):
        bd_temps, bd_heat_flows, z, cond_temps = \
            gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=5,
             radiogenic_heat=[1e-5, 1e-5, 1e-5], surface_t=273,
             heat_flow=0.05, thermal_conductivity=2.5)
        self.assertTrue((bd_temps == [273., 291., 305., 315.]).all())
        self.assertTrue((np.round(bd_heat_flows, 2) == 
                         [0.05, 0.04, 0.03, 0.02]).all())
        self.assertTrue((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
        self.assertTrue((cond_temps == 
                         [273., 291., 305., 315., 315., 315.]).all())
        
    
    
    # Testing on a model with 3 layers of thickness 1 km with variable
    # radiogenic heat production and depth=5
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    def test_cond_geotherm_one_km_layers_variable_A(self):
        bd_temps, bd_heat_flows, z, cond_temps = \
            gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=5,
             radiogenic_heat=[1e-5, 2e-5, 5e-6], surface_t=273,
             heat_flow=0.05, thermal_conductivity=2.5)
        self.assertTrue((bd_temps == [273., 291., 303., 310.]).all())
        self.assertTrue((np.round(bd_heat_flows, 3) == 
                         [0.05, 0.04, 0.02, 0.015]).all())
        self.assertTrue((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
        self.assertTrue((cond_temps == 
                         [273., 291., 303., 310., 310., 310.]).all())
    
    def test_adiab_geotherm(self):
        self.assertEqual(gd.helloworld(),'Hello World')
    
    def test_geotherm(self):
        self.assertEqual(gd.helloworld(),'Hello World')


if __name__ == '__main__':
    unittest.main()