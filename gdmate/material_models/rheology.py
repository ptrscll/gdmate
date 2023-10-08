import numpy as np
import matplotlib.pyplot as plt
import csv

'''
CHANGES TO MAKE
- Make tests (test_rheology.py; every function X in rheology.py should have a test_X function that does asserts)
    - Test rigorously, including for edge cases
    - Check if floats/ints are interchangeable for functions
    - OPTIONAL: Write CREs (or equivalent-  prolly with helpful error messages)
- Make notebook
- Also add gdmate.rheology to the API file (it won't affect the main API file)
- (See code contribution guidelines)
    - Could also edit the markdown file for code contribution

OTHER NOTES
- Don't need to hide adiabatic/conductive geotherm fxns
    - In the future there may be fxns that need to be more private
- May want to talk about some design choices in fxn (ex: do we want those print statements, do we want to plot things)
- Would it make sense to let users choose a depth interval?
'''

def cond_geotherm(thicknesses=[20, 20, 60], depth=600,
             radiogenic_heat=[1.e-6, 2.5e-7, 0.], surface_t=273,
             heat_flow=0.05296, thermal_conductivity=2.5):
    """
    Calculate conductive continental geotherm values.

    This function is based on Chapman86 and Naliboff scripts. It is designed to
    be combined with an adiabatic geotherm (which can be made using 
    adiab_geotherm()). As such, temperatures for the
    asthenosphere are all set to the LAB (Lithosphere-Asthenosphere Boundary) 
    temperature). Calculations are based on a model of the lithosphere with
    discrete layers of set thicknesses with different radiogenic heat 
    production. Heat flow and thermal conductivity are assumed to be constant
    at all depths.
    This function returns the temperatures and heat flows at the boundaries 
    between layers, as well as an array of depths and the temperatures
    calculated for each of those depths. If the user would like, this function
    can also output a graph of the geotherm and save t

    Parameters:
        thicknesses: list of ints
            A list of ints representing the thicknesses of lithospheric units
            in units of kilometers (default: [20, 20, 60]). These ints should
            sum to the total thickness of the lithosphere. The first value is
            the thickness of the uppermost lithospheric unit and each
            subsequent value represents the thickness of the next highest
            layer.

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
            the temperature at the deepest depth. Temperatures in the bottom
            layer (asthenosphere) remain constant.
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
        # Uses equation 6 in Chapman86
        boundary_heat_flows[i + 1] = \
            (boundary_heat_flows[i] - (radiogenic_heat[i] * thick_m[i]))
        
        # Determine temperatures at each layer boundary
        # Uses equation 5 in Chapman86
        boundary_temps[i + 1] = boundary_temps[i] + \
            (boundary_heat_flows[i] / thermal_conductivity) * thick_m[i] - \
            (radiogenic_heat[i] * thick_m[i] ** 2)/(2. * thermal_conductivity)
        

    #### Calculate geotherm ####

    # Making array of depths (in m)
    z = np.arange(0,depth+1,1)*1000

    # Making empty array of temperatures (in K)
    cond_temps = np.zeros(depth+1)
    
    #TODO: Figure out if this part is supposed to be hardcoded to only work with thicknesses of length 3
    # ANSWER: Yes, it can work for any number of thicknesses - but may have to change other code
    #   ex: pandas save code
    #           (this isn't necessary but is definitely helpful)
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

    # Calculate temperature values for each set of depths in each layer
    # For calculations, boundary depths must be multiplied by 1000 to convert
    # from km to m.
    # Calculations are based on equation 4 in Chapman86 
    for i in range(len(thicknesses)):
        temp_layers[i] = boundary_temps[i] + \
            (boundary_heat_flows[i] / thermal_conductivity) * \
                (layers[i] - boundaries[i] * 1000) - \
            (radiogenic_heat[i] * ((layers[i] - boundaries[i] * 1000) ** 2)) / \
                (2*thermal_conductivity)

    # Assign constant temperature to the asthenosphere (bottom layer)
    # To find the temperatures in the asthenosphere, use adiab_geotherm
    temp_layers[-1] = temp_layers[-1] + boundary_temps[-1]

    # Combine temperatures at all depths into single array
    cond_temps = np.concatenate(temp_layers)

    return (boundary_temps, boundary_heat_flows, z, cond_temps)

def adiab_geotherm(z, ast=1573, gravity=9.81, thermal_expansivity=2.e-5,
                   heat_capacity=750, depth=600):
    """
    Function to calculate adiabatic geotherm, based on equation 4.28 in 
    Turcotte and Schubert's Geodynamics textbook. This function is generally
    expected to be called after cond_geotherm() and uses the numpy array of
    depths (z) produced by cond_geotherm(). The function ultimately uses the
    inputted parameters (which are assumed to be constant throughout the 
    asthenosphere) to output an array of temperatures for each of
    the inputted depths based on the adiabatic geotherm. This function is not
    meant to predict temperatures in the lithosphere because the formulas it
    uses only account for adiabatic processes and ignore conductive heat
    transfer.
    
    Parameters:
        z: Numpy array of ints
            Numpy array of depths (in meters). Shallowest depth should be at
            the start of the array; deepest depth should be at the end of the
            array.  This array should be the same as the array returned by 
            conc_geotherm(). There is no default value for this array.

        ast: int                 
            Temperature at adiabatic surface (LAB) (K) used for calculating
            adiabatic geotherm (default: 1573)

        gravity: float            
            Gravitational acceleration (m/s^-2) (default: 9.81)

        thermal_expansivity: float
            Thermal expansivity of asthenosphere (K^-1) (default: 2.e-5)

        heat_capacity: int
            Heat capacity of asthenosphere (J/K*kg) (default: 750)

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
    # By design, the adiabatic surface temperature is the LAB temp
    adiab_temps[0] = ast

    # Iteratively calculating remaining adiabatic temperatures
    for i in range(1, np.size(z)):

        # See line 131 in source/adiabatic_conditions/compute_profile in ASPECT
        #   source code.
        # Numerator of fraction is multiplied by 1000 to account for delta z of
        # 1000 m
        adiab_temps[i] = adiab_temps[i-1] * \
            (1 + (thermal_expansivity * gravity * 1000 * 1./heat_capacity))
        
    return adiab_temps

def geotherm(thicknesses=[20, 20, 60], depth=600,
             radiogenic_heat=[1.e-6, 2.5e-7, 0.], surface_t=273, 
             heat_flow=0.05296, thermal_conductivity=2.5, ast=1573,
             gravity=9.81, thermal_expansivity=2.e-5, heat_capacity=750,
             plot=True, save=True):
    """
    Function to calculate combined conductive and adiabatic geotherm, based on 
    Naliboff scripts. For calculating conductive geotherm, assumes lithosphere
    is divided into discrete layers with differing radiogenic heat production.
    This function otherwise assumes most other variables (ex: thermal 
    conductivity, heat capacity) are constant at all depths.

    This function returns the temperatures and heat flows at the boundaries 
    between layers, as well as an array of depths and the temperatures
    calculated for each of those depths. Temperatures are calculated by adding
    the calculated temperatures from the conductive and adiabatic geotherms.

    In addition to returning the values discussed above, this function also
    prints out the boundary temperatures and heat flows, the bottom temperature,
    and the LAB temperature and depth. This information can be used to quickly
    verify that the inputted parameters are yielding reasonable results.

    Optionally, this function can also output a graph of the geotherm and save
    the printed data to a .csv file called 
    thermal_[lithosphere thickness (km)]_[maximum depth (km)].csv

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
            Temperature at adiabatic surface (LAB) (K) used for calculating
            adiabatic geotherm (default: 1573)

        gravity: float            
            Gravitational acceleration (m/s^-2) (default: 9.81)

        thermal_expansivity: float
            Thermal expansivity of asthenosphere (K^-1) (default: 2.e-5)

        heat_capacity: int
            Heat capacity of asthenosphere (J/K*kg) (default: 750)

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
    
    # Printing relevant information from geotherm calculations
    print('Conductive Boundary Temperatures: ', boundary_temps)
    print('Conductive Boundary Heat Flows: ', boundary_heat_flows)
    
    print('LAB Depth       = ', sum(thicknesses), 'km')
    print('LAB Temperature = ', combined_temps[sum(thicknesses)], 'K')
    print('Bottom Temperature = ',combined_temps[-1], 'K')
    
    # Plotting the temperature and depth arrays if the plot parameter is True
    if plot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(combined_temps, z / 1000)
        ax.invert_yaxis()
        ax.set_xlabel('T (K)')
        ax.set_ylabel('Depth (km)')

    # TODO: Give users control over name of csv file? Change contents of csv file to include full geotherm?
    # Saving data on boundary conditions to csv file if the user desires    
    if save==True:

        # Getting the row and fields
        data = np.concatenate((boundary_temps, boundary_heat_flows[0:-1], 
                               combined_temps[-1]), axis=None)
        # TODO: Rename these fields to be more descriptive?
        fields = ['ts1','ts2','ts3','ts4','qs1','qs2','qs3','base']
        
        # Getting the name of the file
        lith_thickness = np.sum(thicknesses)
        filename = 'thermal_' + str(lith_thickness) + '_' + str(depth) + '.csv'

        # Writing to the csv file

        # writing to csv file  
        with open(filename, 'w') as csvfile:  
            csvwriter = csv.writer(csvfile)  
            csvwriter.writerow(fields)
            csvwriter.writerow(data)

    return (boundary_temps, boundary_heat_flows, z, combined_temps, cond_temps,
            adiab_temps)

geotherm()
