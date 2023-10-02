import numpy as np
import matplotlib.pyplot as plt

'''
CHANGES TO MAKE
- Make tests (test_rheology.py; every function X in rheology.py should have a test_X function that does asserts)
    - Test rigorously, including for edge cases
    - Check if floats/ints are interchangeable for functions
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
    temperature). Calculations are based on a model of the lithosphere with
    discrete layers of set thicknesses with different radiogenic heat 
    production. Heat flow and thermal conductivity are assumed to be constant
    at all depths. 
    This function returns the temperatures and heat flows at the boundaries 
    between layers, as well as an array of depths and the temperatures
    calculated for each of those depths.
    TODO: What does "after Chapman86 and Naliboff scripts" mean. (this also comes up in the geotherm() function)
    
    Parameters:
        thicknesses: list of ints
            A list of ints representing the thicknesses of lithospheric units
            in units of kilometers (default: [20, 20, 60])

        depth: int
            Maximum depth of model (km) (default: 600)

        radiogenic_heat: list of floats    
            A list of floats containing the radiogenic heat production (W/m^3)
            of each lithospheric unit. The list should have same length as 
            thicknesses. (default: [1.e-6, 2.5e-7, 0.])

        surface_t: int           
            Surface temperature (K) (default: 273)

        heat_flow: float
            Surface heat flow (W/m^3) (default: 0.05296)

        thermal_conductivity: float
            Thermal conductivity (W/m*K). This value is assumed to be the same
            for all lithospheric units (default: 2.5)
    
    Returns:
        boundary_temps: Numpy array of floats
            Numpy array containing conductive temperatures (K) at each layer 
            boundary. Values are ordered from least to greatest depth, with
            the first value being the surface temperature and the last value 
            being the temperature at the bottom of the deepest layer. This 
            array has a length one greater than the length of thicknesses.

        boundary_heat_flows: Numpy array of floats
            Numpy array containing heat flows (W/m^3) at each layer boundary.
            Values are ordered from least to greatest depth, with the first
            alue being the surface heat flow and the last value being the heat
            flow at the bottom of the deepest layer. This array has a length
            one greater than the length of thicknesses.

        z: Numpy array of ints
            Numpy array of depths in meters (not kilometers). Depths are spaced
            1000 m apart and start from 0 and end at the maximum depth given by
            the parameter depth.

        cond_temps: Numpy array of floats
            Numpy array of conductive temperatures (K) at each depth given in z.
            First temperature is the surface temperature, last temperature is 
            the temperature at the deepest depth.
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
    Function to calculate adiabatic geotherm. This function is generally 
    expected to be called after cond_geotherm() and uses the numpy array of
    depths (z) produced by cond_geotherm(). The function ultimately uses the
    inputted parameters (which are assumed to be constant throughout the 
    asthenosphere) to output an array of temperatures for each of
    the inputted depths based on the adiabatic geotherm. This function will
    not be likely to produce accurate values for the lithosphere because the
    adiabatic surface temperature is likely to be much higher than the actual
    surface temperature. [TODO: Explain why this is]
    
    Parameters:
        z: Numpy array of ints
            Numpy array of depths (in meters). Shallowest depth should be at
            the start of the array; deepest depth should be at the end of the
            array.  This array should be the same as the array returned by 
            conc_geotherm(). There is no default value for this array.

        ast: int                 
            Surface temperature (K) used for calculating adiabatic geotherm.
            (default: 1573)

        gravity: float            
            Gravitational acceleration (m/s^-2) (default: 9.81)

        thermal_expansivity: float
            Thermal expansivity of asthenosphere (K^-1) (default: 2.e-5)

        heat_capacity: int
            Heat capacity of asthenosphere [TODO: is this right?] (J/K*kg) (default: 750)

        depth: int              
            Maximum depth of model (km). This value should be the same as the 
            depth used in cond_geotherm() and should equal the last value in z
            divided by 1000 (since depth and z have different units). 
            (default: 600)
    
    Returns:
        adiab_temps: Numpy array of floats
            Numpy array of conductive temperatures (K) at each depth given in z.
            First temperature is the surface temperature, last temperature is 
            the temperature at the deepest depth.
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
    Function to calculate combined conductive and adiabatic geotherm, after 
    Naliboff scripts. For calculating conductive geotherm, assumes lithosphere
    is divided into discrete layers with differing radiogenic heat production.
    This function otherwise assumes most other variables (ex: thermal 
    conductivity, heat capacity) are constant at all depths.

    This function returns the temperatures and heat flows at the boundaries 
    between layers, as well as an array of depths and the temperatures
    calculated for each of those depths. Temperatures are calculated by adding
    the calculated temperatures from the conductive and adiabatic geotherms.

    TODO: Add discussion of save/plot options if needed

    Parameters:
        thicknesses: list of ints
            A list of ints representing the thicknesses of lithospheric units
            in units of kilometers (default: [20, 20, 60])

        depth: int
            Maximum depth of model (km) (default: 600)

        radiogenic_heat: list of floats    
            A list of floats containing the radiogenic heat production (W/m^3)
            of each lithospheric unit. The list should have same length as 
            thicknesses. (default: [1.e-6, 2.5e-7, 0.])

        surface_t: int           
            Surface temperature (K) (default: 273)

        heat_flow: float
            Surface heat flow (W/m^3) (default: 0.05296)

        thermal_conductivity: float
            Thermal conductivity (W/m*K). This value is assumed to be the same
            for all lithospheric units (default: 2.5)

        ast: int                 
            Adiabatic surface temperature (K) (default: 1573)

        gravity: float            
            Gravitational acceleration (m/s^-2) (default: 9.81)

        thermal_expansivity: float
            Thermal expansivity of asthenosphere (K^-1) (default: 2.e-5)

        heat_capacity: int
            Heat capacity of asthenosphere [TODO: is this right?] (J/K*kg) (default: 750)

        plot: bool            
            Boolean (True or False) indicating whether to produce a plot of the
            geotherm. (default: True)

        save: bool                
            Boolean indicating whether to save boundary temperatures and heat 
            flows to separate csv file. (default: True)
    
    Returns:
        boundary_temps: Numpy array of floats
            Numpy array containing conductive temperatures (K) at each layer 
            boundary. Values are ordered from least to greatest depth, with
            the first value being the surface temperature and the last value 
            being the temperature at the bottom of the deepest layer. This 
            array has a length one greater than the length of thicknesses.

        boundary_heat_flows: Numpy array of floats
            Numpy array containing heat flows (W/m^3) at each layer boundary.
            Values are ordered from least to greatest depth, with the first
            alue being the surface heat flow and the last value being the heat
            flow at the bottom of the deepest layer. This array has a length
            one greater than the length of thicknesses.

        z: Numpy array of ints
            Numpy array of depths in meters (not kilometers). Depths are spaced
            1000 m apart and start from 0 and end at the maximum depth given by
            the parameter depth.

        combined_temps: Numpy array of floats
            Numpy array containing the temperatures (K) of the combined 
            conductive and adiabatic geotherm at each depth given in z.

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
