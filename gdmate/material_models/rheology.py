import numpy as np
import matplotlib.pyplot as plt

'''
CHANGES TO MAKE
- Improve overall function descriptions in docstrings
- Make tests (test_rheology.py; every function X in rheology.py should have a test_X function that does asserts)
    - Test rigorously, including for edge cases
    - OPTIONAL: Write CREs (or equivalent-  prolly with helpful error messages)
- Make notebook
- Also add gdmate.rheology to the API file
- (See code contribution guidelines)
    - Could also edit the markdown file for code contribution

OTHER NOTES
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
    TODO: What does "after Chapman86 and Naliboff scripts" mean. (this also comes up in the geotherm() function)
    
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

        cond_temps:          Numpy array of conductive temperatures (K) at each
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
    cond_temps = np.zeros(depth+1)
    
    #TODO: Figure out if this part is supposed to be hardcoded to only work with thicknesses of length 3
    # Set boundary locations for slicing cond_temps
    boundaries = [0,thicknesses[0],thicknesses[0]+thicknesses[1],
                  thicknesses[0]+thicknesses[1]+thicknesses[2],
                  depth+1]
    
    # Split depth and temperature arrays into layers based on layer thicknesses
    layers = []
    temp_layers = []
    for i in range(len(thicknesses) + 1):
        layers.append(z[boundaries[i] : boundaries[i + 1]])
        temp_layers.append(cond_temps[boundaries[i] : boundaries[i + 1]])

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
    cond_temps = np.concatenate(temp_layers)

    return (boundary_temps, boundary_heat_flows, z, cond_temps)

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

def geotherm(thicknesses=[20, 20, 60], depth=600,
             radiogenic_heat=[1.e-6, 2.5e-7, 0.], surface_t=273, 
             heat_flow=0.05296, thermal_conductivity=2.5, ast=1573,
             gravity=9.81, thermal_expansivity=2.e-5, heat_capacity=750,
             plot=True, save=True):
    """
    Calculate combined conductive and adiabatic geotherm, after Naliboff
    scripts.

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

        ast:                  Adiabatic surface temperature (K)

        gravity:              Gravitational acceleration (m/s^-2)

        thermal_expansivity:  Thermal expansivity of asthenosphere (K^-1)

        heat_capacity:        Heat capacity of asthenosphere [TODO: is this right?] (J/K*kg)

        plot:                 Boolean indicating whether to produce a plot of
                              the geotherm. A value of True indicates that a
                              plot should be produced.

        save:                 Boolean indicating whether to save boundary 
                              temperatures and heat flows to separate csv file.
                              A value of True indicates these values should be
                              saved.
    
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

        combined_temps:      Numpy array containing the temperatures (K) of the
                             combined conductive and adiabatic geotherm at each
                             depth given in z

    """

    # Calculate conductive geotherm
    boundary_temps, boundary_heat_flows, z, cond_temps = cond_geotherm( \
        thicknesses=thicknesses, depth=depth, radiogenic_heat=radiogenic_heat, \
        surface_t=surface_t, heat_flow=heat_flow, \
        thermal_conductivity=thermal_conductivity)
    
    # Calculate adiabatic geotherm
    adiab_temps = adiab_geotherm(z=z, ast=ast, gravity=gravity, 
                                 thermal_expansivity=thermal_expansivity, 
                                 heat_capacity=heat_capacity, depth=depth)
    
    # Calculate combined geotherm
    combined_temps = cond_temps + adiab_temps - ast
    
    # TODO: Consider deleting
    # Printing relevant information from geotherm calculations
    print('Conductive Boundary Temperatures: ', boundary_temps)
    print('Conductive Boundary Heat Flows: ', boundary_heat_flows)
    
    print('LAB Depth       = ', sum(thicknesses), 'km')
    print('LAB Temperature = ', combined_temps[sum(thicknesses)], 'K')
    print('Bottom Temperature = ',combined_temps[-1], 'K')
    
    # TODO: Consider getting rid of this
    #   If this is removed, no need for matplotlib
    if plot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(combined_temps, z / 1000)
        ax.invert_yaxis()
        ax.set_xlabel('T (K)')
        ax.set_ylabel('Depth (km)')
    
    '''
    TODO: Do we keep this? (If yes, we need to add pandas to gdmate. We also may want to let user decide filename)
    if save==True:
        output = pd.Series(data=np.concatenate((temps,heat_flows[0:-1],tt[-1]),axis=None),
                           index=['ts1','ts2','ts3','ts4','qs1','qs2','qs3','base'])
        
        lith = np.sum(thicknesses)
        
        filename = 'thermal_' + str(lith) + '_' + str(depth) + '.csv'
                                     
        output.to_csv(filename)

    I commented this out for now because pandas is not a part of gdmate atm
    '''
    
    return (boundary_temps, boundary_heat_flows, z, combined_temps, cond_temps,
            adiab_temps)
