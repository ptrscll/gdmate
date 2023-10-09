"""
Test for rheology.py
"""

import gdmate as gd
import numpy as np


def test_cond_geotherm():

    # Testing on a model with 3 layers of thickness 1 km with A = 0 and depth=3
    # Other values have been adjusted to make math cleaner
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=3,
            radiogenic_heat=[0., 0., 0.], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)
    
    assert((bd_temps == [273., 293., 313., 333.]).all())
    assert((bd_heat_flows == [0.05, 0.05, 0.05, 0.05]).all())
    assert((z == [0, 1000, 2000, 3000]).all())
    assert((cond_temps == [273., 293., 313., 333.]).all())


    # Testing on a model with 3 layers of thickness 1 km with A = 0 and depth=5
    # Other values have been adjusted to make math cleaner
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=5,
            radiogenic_heat=[0., 0., 0.], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)

    assert((bd_temps == [273., 293., 313., 333.]).all())
    assert((bd_heat_flows == [0.05, 0.05, 0.05, 0.05]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 293., 313., 333., 333., 333.]).all())


    # Testing on a model with 3 layers of thickness 1 km with A = 1e-5, depth=5
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=5,
            radiogenic_heat=[1e-5, 1e-5, 1e-5], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)

    assert((bd_temps == [273., 291., 305., 315.]).all())
    assert((np.round(bd_heat_flows, 2) == [0.05, 0.04, 0.03, 0.02]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 291., 305., 315., 315., 315.]).all())
    

    # Testing on a model with 3 layers of thickness 1 km with variable
    # radiogenic heat production and depth=5
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[1, 1, 1], depth=5,
            radiogenic_heat=[1e-5, 2e-5, 5e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)
    
    assert((bd_temps == [273., 291., 303., 310.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.04, 0.02, 0.015]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 291., 303., 310., 310., 310.]).all())

    # Testing on a model with 3 layers of thickness 1-3 km with variable
    # radiogenic heat production and depth=8
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[3, 1, 2], depth=8,
            radiogenic_heat=[1e-5, 5e-6, 5e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)
    
    assert((bd_temps == [273., 315., 322., 330.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.02, 0.015, 0.005]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]).all())
    assert((cond_temps == 
            [273., 291., 305., 315., 322., 327., 330., 330., 330.]).all())

    # Testing on a model with 1 layer of thickness 3 km and depth=5
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[3], depth=5,
            radiogenic_heat=[1e-5], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)
    
    assert((bd_temps == [273., 315.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.02]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 291., 305., 315., 315., 315.]).all())

    # Testing on a model with 2 layers of thickness 1-3 km with variable
    # radiogenic heat production and depth=6
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[3, 1], depth=6,
            radiogenic_heat=[1e-5, 5e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)
    
    assert((bd_temps == [273., 315., 322.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.02, 0.015]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000, 6000]).all())
    assert((cond_temps == 
            [273., 291., 305., 315., 322., 322., 322.]).all())


    # Testing on a model with 4 layers of thickness 1-3 km with variable
    # radiogenic heat production and depth=9
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, cond_temps = \
        gd.rheology.cond_geotherm(thicknesses=[3, 1, 2, 2], depth=9,
            radiogenic_heat=[1e-5, 5e-6, 5e-6, 1e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5)
    
    assert((bd_temps == [273., 315., 322., 330., 333.2]).all())
    assert((np.round(bd_heat_flows, 3) == 
            [0.05, 0.02, 0.015, 0.005, 0.003]).all())
    assert((z == 
            [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]).all())
    assert((cond_temps == 
            [273., 291., 305., 315., 322., 327., 330., 331.8, 333.2, 
             333.2]).all())

#TODO: This
def test_adiab_geotherm():
    assert(gd.helloworld() == 'Hello World')

#TODO: This
def test_geotherm():
    assert(gd.helloworld() == 'Hello World')