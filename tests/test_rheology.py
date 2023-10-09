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

def test_adiab_geotherm():
    # Testing adiab_geotherm on a small scale with easy-to-handle numbers
    adiab_temps = gd.rheology.adiab_geotherm([0, 1000, 2000, 3000], ast=1573,
                                              gravity=10.0, 
                                              thermal_expansivity=0.01,
                                              heat_capacity=100, depth=3)
    assert((adiab_temps == [1573., 3146., 6292., 12584.]).all())

#TODO: This
def test_geotherm():
    
    # Testing on a model with 3 layers of thickness 1 km with A = 0, depth=3
    # and an adiabatic surface temperature of 0 K (i.e, no adiabatic geotherm)
    # Other values have been adjusted to make math cleaner
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[1, 1, 1], depth=3,
            radiogenic_heat=[0., 0., 0.], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=0, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)
    
    assert((bd_temps == [273., 293., 313., 333.]).all())
    assert((bd_heat_flows == [0.05, 0.05, 0.05, 0.05]).all())
    assert((z == [0, 1000, 2000, 3000]).all())
    assert((cond_temps == [273., 293., 313., 333.]).all())
    assert((adiab_temps == [0.0, 0.0, 0.0, 0.0]).all())
    assert((combined_temps == cond_temps).all())


    # Testing on a model with 3 layers of thickness 1 km with A = 0 and depth=5
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[1, 1, 1], depth=5,
            radiogenic_heat=[0., 0., 0.], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)

    assert((bd_temps == [273., 293., 313., 333.]).all())
    assert((bd_heat_flows == [0.05, 0.05, 0.05, 0.05]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 293., 313., 333., 333., 333.]).all())
    assert((adiab_temps == [1., 2., 4., 8., 16., 32.]).all())
    assert((combined_temps == [273., 294., 316., 340., 348., 364.]).all())


    # Testing on a model with 3 layers of thickness 1 km with A = 1e-5, depth=5
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[1, 1, 1], depth=5,
            radiogenic_heat=[1e-5, 1e-5, 1e-5], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)

    assert((bd_temps == [273., 291., 305., 315.]).all())
    assert((np.round(bd_heat_flows, 2) == [0.05, 0.04, 0.03, 0.02]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 291., 305., 315., 315., 315.]).all())
    assert((adiab_temps == [1., 2., 4., 8., 16., 32.]).all())
    assert((combined_temps == [273., 292., 308., 322., 330., 346.]).all())
    

    # Testing on a model with 3 layers of thickness 1 km with variable
    # radiogenic heat production and depth=5
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[1, 1, 1], depth=5,
            radiogenic_heat=[1e-5, 2e-5, 5e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)
    
    assert((bd_temps == [273., 291., 303., 310.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.04, 0.02, 0.015]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 291., 303., 310., 310., 310.]).all())
    assert((adiab_temps == [1., 2., 4., 8., 16., 32.]).all())
    assert((combined_temps == [273., 292., 306., 317., 325., 341.]).all())


    # Testing on a model with 3 layers of thickness 1-3 km with variable
    # radiogenic heat production and depth=8
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[3, 1, 2], depth=8,
            radiogenic_heat=[1e-5, 5e-6, 5e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)
    
    assert((bd_temps == [273., 315., 322., 330.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.02, 0.015, 0.005]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]).all())
    assert((cond_temps == 
            [273., 291., 305., 315., 322., 327., 330., 330., 330.]).all())
    assert((adiab_temps == [1., 2., 4., 8., 16., 32., 64., 128., 256.]).all())
    assert((combined_temps == 
            [273., 292., 308., 322., 337., 358., 393., 457., 585.]).all())

            
    # Testing on a model with 1 layer of thickness 3 km and depth=5
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[3], depth=5,
            radiogenic_heat=[1e-5], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)
    
    assert((bd_temps == [273., 315.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.02]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000]).all())
    assert((cond_temps == [273., 291., 305., 315., 315., 315.]).all())
    assert((adiab_temps == [1., 2., 4., 8., 16., 32.]).all())
    assert((combined_temps == [273., 292., 308., 322., 330., 346.]).all())


    # Testing on a model with 2 layers of thickness 1-3 km with variable
    # radiogenic heat production and depth=6
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[3, 1], depth=6,
            radiogenic_heat=[1e-5, 5e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)
    
    assert((bd_temps == [273., 315., 322.]).all())
    assert((np.round(bd_heat_flows, 3) == [0.05, 0.02, 0.015]).all())
    assert((z == [0, 1000, 2000, 3000, 4000, 5000, 6000]).all())
    assert((cond_temps == [273., 291., 305., 315., 322., 322., 322.]).all())
    assert((adiab_temps == [1., 2., 4., 8., 16., 32., 64.]).all())
    assert((combined_temps == [273., 292., 308., 322., 337., 353., 385.]).all())


    # Testing on a model with 4 layers of thickness 1-3 km with variable
    # radiogenic heat production and depth=9
    # Adiabatic geotherm has been fixed to start at 1 K and double every 1 km
    # Other values have been adjusted to make math cleaner (albeit unrealistic)
    bd_temps, bd_heat_flows, z, combined_temps, cond_temps, adiab_temps = \
        gd.rheology.geotherm(thicknesses=[3, 1, 2, 2], depth=9,
            radiogenic_heat=[1e-5, 5e-6, 5e-6, 1e-6], surface_t=273,
            heat_flow=0.05, thermal_conductivity=2.5, ast=1, gravity=10.0,
            thermal_expansivity=0.01, heat_capacity=100, plot=False, save=False)
    
    assert((bd_temps == [273., 315., 322., 330., 333.2]).all())
    assert((np.round(bd_heat_flows, 3) == 
            [0.05, 0.02, 0.015, 0.005, 0.003]).all())
    assert((z == 
            [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]).all())
    assert((cond_temps == 
            [273., 291., 305., 315., 322., 327., 330., 331.8, 333.2, 
             333.2]).all()) 
    assert((adiab_temps == 
            [1., 2., 4., 8., 16., 32., 64., 128., 256., 512.]).all())
    assert((combined_temps == 
            [273., 292., 308., 322., 337., 358., 393., 458.8, 588.2, 
             844.2]).all())
   