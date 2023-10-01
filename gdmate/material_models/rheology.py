import numpy as np
import matplotlib.pyplot as plt

'''
CHANGES TO MAKE
- Make sure code is written in most efficient way
- Make sure code is very well documented (see pyvista_vis - but give longer descriptions for function in docstring (fxn contract))
- Make tests (test_rheology.py; every function X in rheology.py should have a test_X function that does asserts)
    - Test rigorously, including for edge cases
    - OPTIONAL: Write CREs (or equivalent-  prolly with helpful error messages)
- Make notebook
- Also add gdmate.rheology to the API file
- (See code contribution guidelines)
    - Could also edit the markdown file for code contribution
- Don't need to hide adiabatic/conductive geotherm fxns
    - In the future there may be fxns that need to be more private
- May want to talk about some design choices in fxn (ex: do we want those print statements, do we want to plot things)
'''

def cond_geotherm(thicknesses=[20, 20, 60], depth=600,
             radiogenic_heat=[1.e-6, 2.5e-7, 0.], surface_t=273,
             heat_flow=0.05296, thermal_conductivity=2.5):
    """
    Calculate conductive continental geotherm values
    after Chapman86 and Naliboff scripts. Designed to be combined with
    adiabatic geotherm (i.e., asthenosphere temperature set as LAB
    temperature).
    
    Parameters:
        thicknesses:          List of ints representing the thicknesses of 
                              lithospheric units (km)

        depth:                Maximum depth of model (km)

        radiogenic_heat:      List of floats containing radiogenic heat 
                              production (W/m^3) of each lithospheric unit.
                              List should have same length as thicknesses

        surface_t:            Surface temperature (K)

        heat_flow:            Surface heat flow (W/m^3)

        thermal_conductivity: Thermal conductivity (W/m*K). Passed as a single
                              float (this value is assumed to be the same for 
                              all lithospheric units)
    
    Returns:
        boundary_temps:      Numpy array containing conductive temperatures (K)
                             at each layer boundary. First value is the
                             surface temperature, last value is the 
                             temperature at the bottom of the deepest layer.

        boundary_heat_flows: Numpy array containing heat flows (W/m^3) at each 
                             layer boundary. First value is the surface
                             heat flow, last value is the heat flow at the 
                             bottom of the deepest layer.

        z:                   Numpy array of depths (m). Depths are spaced
                             1000 m apart.

        tc:                  Numpy array of conductive temperatures (K) at each
                             depth given in z
    """

    # Convert thicknesses to meters
    thick_m = np.array(thicknesses) * 1000
    

    #### Set up arrays of heat flows and temperatures ####

    # Each value represents heat flow/temperature at a boundary between layers
    # Units: W/m^3 for heat_flows, K for temps
    boundary_heat_flows = np.zeros(len(thick_m) + 1)
    boundary_temps = np.zeros(len(thick_m) + 1)

    # At index 0, store the heat flow and temperature at the surface
    boundary_heat_flows[0] = heat_flow
    boundary_temps[0] = surface_t

    # Iteratively calculate heat flow and temperature at each layer boundary
    for i in range(len(thicknesses)):
    
        # Determine heat flows at each layer boundary
        boundary_heat_flows[i + 1] = \
            (boundary_heat_flows[i] - (radiogenic_heat[i] * thick_m[i]))
        
        # Determine temperatures at each layer boundary
        boundary_temps[i + 1] = boundary_temps[i] + \
            (boundary_heat_flows[i] / thermal_conductivity) * thick_m[i] - \
            (radiogenic_heat[i] * thick_m[i] ** 2)/(2. * thermal_conductivity)
        

    #### Calculate geotherm ####

    # Making array of depths (in m)
    z = np.arange(0,depth+1,1)*1000

    # Making empty array of temperatures (in K)
    all_temps = np.zeros(depth+1)
    
    #TODO: Figure out if this part is supposed to be hardcoded to only work with thicknesses of length 3
    # Set boundary locations for slicing tc
    boundaries = [0,thicknesses[0],thicknesses[0]+thicknesses[1],
                  thicknesses[0]+thicknesses[1]+thicknesses[2],
                  depth+1]
    
    # Split depth and temperature arrays into layers based on layer thicknesses
    layers = []
    temp_layers = []
    for i in range(len(thicknesses) + 1):
        layers.append(z[boundaries[i] : boundaries[i + 1]])
        temp_layers.append(all_temps[boundaries[i] : boundaries[i + 1]])

    # Assign appropriate temperature values for each set of depths
    for i in range(len(thicknesses)):
        temp_layers[i] = boundary_temps[i] + \
            (boundary_heat_flows[i] / thermal_conductivity) * \
                (layers[i] - boundaries[i] * 1000) - \
            (radiogenic_heat[i] * ((layers[i] - boundaries[i] * 1000) ** 2)) / \
                (2*thermal_conductivity)

    # Assign constant temperature to the asthenosphere (bottom layer)
    # To find the temperatures in the astenosphere, use adiab_geotherm
    temp_layers[-1] = temp_layers[-1] + boundary_temps[-1]

    # Combine temperatures at all depths into single array
    all_temps = np.concatenate(temp_layers)

    return (boundary_temps, boundary_heat_flows, z, all_temps)

def adiab_geotherm(z, ast=1573, gravity=9.81, thermal_expansivity=2.e-5,
                   heat_capacity=750, depth=600):
    """
    Calculate adiabatic geotherm. Assumes numpy array of depths (z) has
    already been calculated using conc_geotherm()
    
    Parameters:
        z:                   Numpy array of depths (m). Should be the same as
                             the one calculated in conc_geotherm()

        ast:                 Adiabatic surface temperature (K)

        gravity:             Gravitational acceleration (m/s^-2)

        thermal_expansivity: Thermal expansivity of asthenosphere (K^-1)

        heat_capacity:       Heat capacity of asthenosphere [TODO: is this right?] (J/K*kg)

        depth:               Maximum depth of model (km)
    
    Returns:
        adiab_temps: Adiabatic temperature for each depth (K). First
                     temperature is the surface temperature, last temperature
                     is the temperature at the deepest depth.
    """

    # Make empty array of adiabatic temperatures
    adiab_temps = np.zeros(depth+1)

    # Set first index to adiabatic surface temperature
    # By design, the adiabatic surface temperature is the LAB temp <-- TODO: What does this mean?
    adiab_temps[0] = ast

    # Iteratively calculating remaining adiabatic temperatures
    for i in range(1, np.size(z)):

        # See line 124 in source/adiabatic_conditions/compute_profile
        adiab_temps[i] = adiab_temps[i-1] * (
            1 + (thermal_expansivity * gravity * 1000 * 1./heat_capacity))
        
    return adiab_temps

def geotherm(thicknesses=[20,20,60],depth=600,
             radiogenic_heat=[1.e-6,2.5e-7,0.],surface_t=273,
             heat_flow=0.05296,thermal_conductivity=2.5,ast=1573,gravity=9.81,
             thermal_expansivity=2.e-5,heat_capacity=750,plot=True,
             save=True):
    """
    Calculate combined conductive and adiabatic geotherm, after Naliboff
    scripts.
    
    Parameters:
        thicknesses: Thicknesses of lithospheric units (km)
        depth: Depth of model (km)
        radiogenic_heat: Radiogenic heat production in each unit (W/m^3)
        surface_t: Surface temperature (K)
        heat_flow: Surface heat flow (W/m^3)
        thermal_conductivity: Thermal conductivity (W/m*K)
        ast: Adiabatic surface temperature (K)
        gravity: Gravitational acceleration (m/s^-2)
        thermal_expansivity: Thermal expansivity (K^-1)
        heat_capacity: Heat capacity (J/K*kg)
        plot: Whether to plot the geotherm
        save: Whether to save boundary temperatures and heat flows to
            separate csv file.
        
    Returns:
        temps: Conductive temperatures (K) at each layer boundary
        heat_flows: Heat flows (W/m^3) at each layer boundary
        z: Array of depths (m)
        tt: Temperature (K) of combined conductive and adiabatic geotherm at
            each depth.
        
    """
    # Conductive geotherm
    temps,heat_flows,z,tc = cond_geotherm(thicknesses=thicknesses,depth=depth,
                       radiogenic_heat=radiogenic_heat,surface_t=surface_t,
                       heat_flow=heat_flow,thermal_conductivity=
                       thermal_conductivity)
    
    # Adiabatic geotherm
    ta = adiab_geotherm(z=z,ast=ast,gravity=gravity,thermal_expansivity=
                        thermal_expansivity,heat_capacity=heat_capacity,
                        depth=depth)
    
    # Combined geotherm
    tt = tc + ta - ast
    
    # TODO: Consider deleting
    print('Conductive Boundary Temperatures: ', temps)
    print('Conductive Boundary Heat Flows: ', heat_flows)
    
    print('LAB Depth       = ', sum(thicknesses), 'km')
    print('LAB Temperature = ', tt[sum(thicknesses)], 'K')
    print('Bottom Temperature = ',tt[-1], 'K')
    
    # TODO: Consider getting rid of this
    #   If this is removed, no need for matplotlib
    if plot==True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tt,z/1000)
        ax.invert_yaxis()
        ax.set_xlabel('T (K)')
        ax.set_ylabel('Depth (km)')
        # TODO: Delete this later
        #plt.savefig("old_geotherm.png")
    
    '''
    if save==True:
        output = pd.Series(data=np.concatenate((temps,heat_flows[0:-1],tt[-1]),axis=None),
                           index=['ts1','ts2','ts3','ts4','qs1','qs2','qs3','base'])
        
        lith = np.sum(thicknesses)
        
        filename = 'thermal_' + str(lith) + '_' + str(depth) + '.csv'
                                     
        output.to_csv(filename)

    I commented this out for now because pandas is not a part of gdmate atm
    '''
    
    return(temps,heat_flows,z,tt,tc,ta)

temps,heat_flows,z,tc = cond_geotherm()
print(adiab_geotherm(z))
#print(geotherm(plot=False))
#geotherm(plot=True)
