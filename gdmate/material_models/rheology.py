import numpy as np
import matplotlib.pyplot as plt
import csv

def cond_geotherm(thicknesses=[20, 20, 60], depth=600,
             radiogenic_heat=[1.e-6, 2.5e-7, 0.], surface_t=273,
             heat_flow=0.05296, thermal_conductivity=2.5):
    """
    Calculate conductive continental geotherm values.

    This function's calculations are based on scripts by John Naliboff and
    equations from a 1986 article by D.S. Chapman. To use these equations,
    this function employs a model of the lithosphere with
    discrete layers of set thicknesses with different radiogenic heat 
    production. Thermal conductivity is assumed to be constant
    at all depths.
    
    This function is designed to be combined with an adiabatic geotherm (which
    can be found using adiab_geotherm()). As such, all temperatures calculated
    for the asthenosphere are set to the LAB (Lithosphere-Asthenosphere
    Boundary) temperature. 

    This function returns the temperatures and heat flows at the boundaries 
    between layers, as well as an array of depths and the temperatures
    calculated for each of those depths.

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
            1000 m apart. Depths start from 0 and end at the maximum depth given
            by the parameter depth.

        cond_temps: Numpy array of floats
            Numpy array of conductive temperatures (K) at each depth given in z.
            First temperature is the surface temperature, last temperature is 
            the temperature at the deepest depth. Temperatures in the
            asthenosphere remain constant.
    """

    # Convert thicknesses to meters
    thick_m = np.array(thicknesses) * 1000


    #### Set up arrays of heat flows and temperatures ####

    # Each value in one of these arrays represents heat flow/temperature at a 
    #   boundary between layers
    # Units: W/m^3 for heat_flows, K for temps
    boundary_heat_flows = np.zeros(len(thick_m) + 1)
    boundary_temps = np.zeros(len(thick_m) + 1)

    # At index 0, store the heat flow and temperature at the surface
    boundary_heat_flows[0] = heat_flow
    boundary_temps[0] = surface_t

    # Iteratively calculate heat flow and temperature at each layer boundary
    for i in range(len(thicknesses)):
    
        # Determine heat flows at each layer boundary
        # Uses equation 6 in Chapman 1986
        boundary_heat_flows[i + 1] = \
            (boundary_heat_flows[i] - (radiogenic_heat[i] * thick_m[i]))
        
        # Determine temperatures at each layer boundary
        # Uses equation 5 in Chapman 1986
        boundary_temps[i + 1] = boundary_temps[i] + \
            (boundary_heat_flows[i] / thermal_conductivity) * thick_m[i] - \
            (radiogenic_heat[i] * thick_m[i] ** 2)/(2. * thermal_conductivity)
        

    #### Calculate geotherm ####

    # Making array of depths (in m)
    z = np.arange(0,depth+1,1)*1000

    # Make empty array of temperatures (in K)
    cond_temps = np.zeros(depth+1)
    
    # Make empty array of depths of boundaries (in km) between layers
    boundaries = np.zeros(len(thicknesses) + 2, dtype=int)

    # The first boundary is the surface
    boundaries[0] = 0

    # The next boundaries are found between each lithospheric layer
    for i in range(len(thicknesses)):
        boundaries[i + 1] = boundaries[i] + thicknesses[i]

    # The bottom-most boundary is the maximum depth
    boundaries[-1] = depth + 1
    
    # Split depth and temperature arrays into layers based on layer boundaries
    layers = []
    temp_layers = []
    for i in range(len(thicknesses) + 1):
        layers.append(z[boundaries[i] : boundaries[i + 1]])
        temp_layers.append(cond_temps[boundaries[i] : boundaries[i + 1]])

    # Calculate temperature values for each set of depths in each layer
    # For calculations, boundary depths must be multiplied by 1000 to convert
    # from km to m.
    # Calculations are based on equation 4 in Chapman 1986 
    for i in range(len(thicknesses)):
        temp_layers[i] = boundary_temps[i] + \
            (boundary_heat_flows[i] / thermal_conductivity) * \
                (layers[i] - boundaries[i] * 1000) - \
            (radiogenic_heat[i] * ((layers[i] - boundaries[i] * 1000) ** 2)) / \
                (2*thermal_conductivity)

    # Assign the same temperature to the whole bottom layer (the asthenosphere)
    # To find the temperatures in the asthenosphere, use adiab_geotherm
    temp_layers[-1] = temp_layers[-1] + boundary_temps[-1]

    # Combine temperatures at all depths back into single array
    cond_temps = np.concatenate(temp_layers)

    return (boundary_temps, boundary_heat_flows, z, cond_temps)


def adiab_geotherm(z, ast=1573, gravity=9.81, thermal_expansivity=2.e-5,
                   heat_capacity=750, depth=600):
    """
    Function to calculate adiabatic geotherm using equation 4.254 in 
    Turcotte and Schubert's Geodynamics textbook. This function is generally
    expected to be called after cond_geotherm() and uses the numpy array of
    depths (z) produced by cond_geotherm(). The function uses the
    inputted parameters (which are assumed to be constant throughout the 
    asthenosphere) to calculate an array of adiabatic temperatures for each of
    the inputted depths. 
    
    This function is not meant to predict temperatures in the lithosphere 
    because the formulas it uses only account for adiabatic processes and
    ignore conductive heat transfer.
    
    Parameters:
        z: Numpy array of ints
            Numpy array of depths (in meters). Shallowest depth should be at
            the start of the array; deepest depth should be at the end of the
            array.  This array should be the same as the array returned by 
            cond_geotherm(). There is no default value for this array.

        ast: int                 
            Temperature at adiabatic surface (K) used for calculating
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
            Numpy array of adiabatic temperatures (K) at each depth given in z.
            First temperature is the surface temperature, last temperature is 
            the temperature at the deepest depth.
    """

    # Make empty array of adiabatic temperatures
    adiab_temps = np.zeros(depth+1)

    # Set first index to adiabatic surface temperature
    # By design, the adiabatic surface temperature is the LAB temp
    adiab_temps[0] = ast

    # Iteratively calculate remaining adiabatic temperatures
    for i in range(1, np.size(z)):

        # This line is based off line 131 in
        #   source/adiabatic_conditions/compute_profile in ASPECT source code.
        # Numerator of fraction is multiplied by 1000 to account for delta z =
        #   1000 m
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
    scripts by John Naliboff. For calculating conductive geotherm, assumes
    lithosphere is divided into discrete layers with differing radiogenic heat
    production. This function otherwise assumes most other variables (ex:
    thermal conductivity, heat capacity) are constant at all depths.

    This function returns the temperatures and heat flows at the boundaries 
    between layers, as well as an array of depths and the adiabatic, conductive,
    and combined adiabatic and conductive temperatures calculated for each of
    those depths.

    In addition to returning the values discussed above, this function also
    prints out the boundary temperatures and heat flows, the LAB temperature 
    and depth, and the temperature at the maximum depth analyzed by the model. 
    This information can be used to verify that inputted parameters are 
    yielding reasonable results.

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
            Temperature at adiabatic surface (K) used for calculating
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
            value being the surface heat flow and the last value being the heat
            flow at the bottom of the deepest layer. This array has a length
            one greater than the length of thicknesses.

        z: Numpy array of ints
            Numpy array of depths in meters (not kilometers). Depths are spaced
            1000 m apart and start from 0 and end at the maximum depth given by
            the parameter depth.

        combined_temps: Numpy array of floats
            Numpy array containing the combined conductive and adiabatic 
            temperatures (K) at each depth given in z.
        
        cond_temps: Numpy array of floats
            Numpy array containing the conductive temperature at each depth in z
        
        adiab_temps: Numpy array of floats
            Numpy array containing the adiabatic temperature at each depth in z

    """

    # Calculate conductive geotherm using cond_geotherm
    boundary_temps, boundary_heat_flows, z, cond_temps = cond_geotherm( \
        thicknesses=thicknesses, depth=depth, radiogenic_heat=radiogenic_heat, \
        surface_t=surface_t, heat_flow=heat_flow, \
        thermal_conductivity=thermal_conductivity)
    
    # Calculate adiabatic geotherm using adiab_geotherm
    adiab_temps = adiab_geotherm(z=z, ast=ast, gravity=gravity, 
                                 thermal_expansivity=thermal_expansivity, 
                                 heat_capacity=heat_capacity, depth=depth)
    
    # Calculate combined geotherm by adding conductive and adiabatic
    #   temperatures and subtracting the adiabatic surface temperature
    combined_temps = cond_temps + adiab_temps - ast
    
    # Printing relevant information from geotherm calculations
    print('Conductive Boundary Temperatures: ', boundary_temps)
    print('Conductive Boundary Heat Flows: ', boundary_heat_flows)
    
    print('LAB Depth       = ', sum(thicknesses), 'km')
    print('LAB Temperature = ', combined_temps[sum(thicknesses)], 'K')
    print('Bottom Temperature = ',combined_temps[-1], 'K')
    
    # Plotting combined geotherm if the plot parameter is True
    if plot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(combined_temps, z / 1000)
        ax.invert_yaxis()
        ax.set_xlabel('T (K)')
        ax.set_ylabel('Depth (km)')

    # Saving data on boundary conditions to csv file if save parameter is True  
    if save==True:

        # Getting the row of data
        data = np.concatenate((boundary_temps, boundary_heat_flows[0:-1], 
                               combined_temps[-1]), axis=None)
        
        # Getting the fields
        fields = []
        for i in range(len(boundary_temps)):
            fields.append("ts" + str(i + 1))
        for i in range(len(boundary_heat_flows) - 1):
            fields.append("qs" + str(i + 1))
        fields.append('base')
        
        # Getting the name of the file
        lith_thickness = np.sum(thicknesses)
        filename = 'thermal_' + str(lith_thickness) + '_' + str(depth) + '.csv'

        # writing the fields and row of data to the csv file  
        with open(filename, 'w') as csvfile:  
            csvwriter = csv.writer(csvfile)  
            csvwriter.writerow(fields)
            csvwriter.writerow(data)

    return (boundary_temps, boundary_heat_flows, z, combined_temps, cond_temps,
            adiab_temps)


def drucker_prager(pressure, internal_friction=30, cohesion=2e7):
    """
    Calculate failure strength of a material from a given pressure, angle of 
    internal friction, and cohesion using Drucker-Prager criterion.

    Parameters:
        pressure: int or float
            Pressure (Pa)

        internal_friction: int or float
            Angle of internal friction (degrees) (default: 30)

        cohesion: int or float
            Cohesion (Pa) (default: 2e7)
    
    Returns:
        strength: float
            Strength (Pa)

    """
    friction_rad = np.radians(internal_friction)
    failure_strength = pressure * np.sin(friction_rad) + \
        cohesion * np.cos(friction_rad)
    
    return failure_strength


def viscosity(A, n, d, m, E, P, V, T, strain_rate=1e-15, R=8.31451):
    """
    Calculate viscosity of a material according to equation from the ASPECT 
    manual.
    TODO: Does it make sense to still cite ASPECT manual here, or is this formula common knowledge
    (it matches the formula on the course slides but idk if there's alternate versions).
    If we still want to cite the manual, then Where in the manual is this equation?
    TODO: Confirm units for A
    
    For dislocation creep, m = 0. For diffusion creep, n = 1.
    
    Parameters:
        A: float
            Power-law constant (kg * m^-2 * s^-1)

        n: int
            Stress exponent. n = 1 for diffusion creep.

        d: float
            Grain size of the material (m)

        m: int
            Grain size exponent. m = 0 for dislocation creep.

        E: float
            Activation energy (J/mol)

        P: float
            Pressure (Pa)

        V: float
            Activation volume (m^3/mol)

        T: float
            Temperature (K)

        strain_rate: float
            square root of the second invariant of the strain rate tensor (s^-1)
            (default value: 1e-15)
            [TODO: What does this mean?]

        R: float
            Gas constant (J/K*mol) (default value: 8.31451)
        
    Returns:
        visc: float
            Viscosity of the material (Pa*s)
    """
    visc = (0.5 * A**(-1/n) * d**(m/n) * 
            (strain_rate)**((1-n)/n)*np.exp((E+P*V)/(n*R*T)))

    return(visc)


def visc_diffusion(A, d, m, E, P, V, T, strain_rate=1e-15, R=8.31451):
    """
    Calculate viscosity of a material undergoing diffusion creep.
    Calculations are based on equation from the ASPECT manual.
    TODO: Does it make sense to still cite ASPECT manual here, or is this formula common knowledge
    (it matches the formula on the course slides but idk if there's alternate versions).
    If we still want to cite the manual, then Where in the manual is this equation?
    TODO: Confirm units for A
    
    Parameters:
        A: float
            Power-law constant (kg * m^-2 * s^-1)

        d: float
            Grain size of the material (m)

        m: int
            Grain size exponent.

        E: float
            Activation energy (J/mol)

        P: float
            Pressure (Pa)

        V: float
            Activation volume (m^3/mol)

        T: float
            Temperature (K)

        strain_rate: float
            square root of the second invariant of the strain rate tensor (s^-1)
            (default value: 1e-15)
            [TODO: What does this mean?]

        R: float
            Gas constant (J/K*mol) (default value: 8.31451)
        
    Returns:
        visc: float
            Viscosity of the material (Pa*s)
    """
    
    visc = viscosity(A=A, n=1, d=d, m=m, E=E, P=P, V=V, T=T, 
                      strain_rate=strain_rate, R=R)
    
    return(visc)


def visc_dislocation(A, n, E, P, V, T, strain_rate=1e-15, R=8.31451):
    
    """
    Calculate viscosity for a material undergoing dislocation creep.
    Calculations based on equation from the ASPECT manual.
    TODO: Does it make sense to still cite ASPECT manual here, or is this formula common knowledge
    (it matches the formula on the course slides but idk if there's alternate versions).
    If we still want to cite the manual, then Where in the manual is this equation?
    TODO: Confirm units for A
    
    Parameters:
        A: float
            Power-law constant (kg * m^-2 * s^-1)

        n: int
            Stress exponent.

        E: float
            Activation energy (J/mol)

        P: float
            Pressure (Pa)

        V: float
            Activation volume (m^3/mol)

        T: float
            Temperature (K)

        strain_rate: float
            square root of the second invariant of the strain rate tensor (s^-1)
            (default value: 1e-15)
            [TODO: What does this mean?]

        R: float
            Gas constant (J/K*mol) (default value: 8.31451)
        
    Returns:
        visc: float
            Viscosity of the material (Pa*s)
    """
    
    # Value for d (grain size) is irrelevant when dealing with dislocation creep
    # Thus, we arbitrarily assign d to be 1 for this function
    visc = viscosity(A=A, n=n, d=1, m=0, E=E, P=P, V=V, T=T, 
                      strain_rate=strain_rate, R=R)
    
    return(visc)


def visc_composite(visc_dislocation, visc_diffusion):
    """
    Calculate composite viscosity of a material by combining the viscosities
    calculating from the diffusion creep and dislocation creep flow laws. Both
    viscosities must be pre-calculated and passed into this function as
    parameters.
    
    Parameters:
        visc_dislocation: float
            Dislocation creep viscosity (Pa*s). Can be calculated from function
            of the same name.
        
        visc_diffusion: float
            Diffusion creep viscosity (Pa*s). Can be calculated from function of
            the same name.
        
    Returns:
        visc: float
            Composite viscosity of the material (Pa*s)
    """
    
    # TODO: What is the basis of these calculations?
    visc = (
        (visc_dislocation * visc_diffusion) / 
        (visc_dislocation + visc_diffusion)
    )
    
    return(visc)


def adiab_density(input_density, thermal_expansivity, temperature, 
                   reference_temp):
    """
    Calculate array of adiabatic densities at different depths based on
    input densities and temperatures at different depths. The calculated
    density accounts for temperature variations in different depths/layers 
    within the Earth.

    TODO: Find reference source (and finish docstring)
    TODO: Replace with combined/adiabatic temps?

    Parameters:
        input_density: Numpy array of floats
            Input density (kg/m^2)

        thermal_expansivity: Float
            Thermal expansivity (K^-1)

        temperature: Temperature (K) <-- comb_temps in fxn call
        reference_temp: Reference temperature (K) <-- adiab_temps in fxn call
    
    Returns:
        output: Numpy array of floats
            Adiabatic densities (kg/m^2) corresponding to each input density

    """
    
    output = input_density * (
        1 - thermal_expansivity * (temperature - reference_temp))
    
    return(output)


def density_profile(z, thicknesses=[20, 20, 60], 
                    densities=[2800, 2900, 3300, 3300]):
    """
    Create numpy array of densities for every depth in a given depth array (z).
    
    Densities are made uniform for each layer in the model. Thicknesses and 
    densities of each individual layer may be passed into the function as
    parameters. Default values are for (in order) the upper crust, lower crust,
    mantle lithosphere, and mantle asthenosphere.
        
    Parameters:
        z: Numpy array of ints
            Sorted numpy array of depths in meters (not kilometers)

        thicknesses: Numpy array of ints
            Thicknesses (km) of each layer, ordered from top to bottom and
            excluding the bottom layer (which is assumed to extend to the
            greatest depth in z)
        
        densities: Numpy array of ints
            Densities (kg/m^2) for each layer, order from top to bottom (and
            including the bottom layer).
    
    Returns:
        rho: Numpy array of ints
            Numpy array of densities (kg/m^2) corresponding to each depth

    """
    rho = np.zeros(len(z))   # Density array to be returned
    depth_index = 0          # Variable storing the current index in z and rho
    layer_bottom = 0         # Variable storing the bottom of the current layer

    # For each layer, set the densities at every depth in the layer to the
    # same density (as given in the densities array)
    for i in range(len(thicknesses)):
        layer_bottom += thicknesses[i] * 1000
        while z[depth_index] < layer_bottom:
            rho[depth_index] = densities[i]
            depth_index += 1

    # Setting the densities in the depths found in the bottom layer to the last
    # provided density value
    rho[depth_index:] = densities[-1]
    
    return rho


def pressure(z, rho, g=9.81):
    """
    Produce an array of pressures at given depths in the Earth based on array
    of densities at different depths.
    
    Parameters:
        z: Numpy array of ints
            Sorted numpy array of depths in meters (not kilometers). First value
            should be 0, representing Earth's surface.

        rho: Numpy array of floats
            Numpy array of densities (kg/m^2) at each depth given in z.
        
        g: float
            Acceleration due to gravity (m/s^2) (default value: 9.81)
    
    Returns:
        P: Array of pressures (Pa) for each of the depths given in z. Note that
        first value is fixed to 0.
    """
    P = np.zeros(len(z))

    # Iteratively calculating the pressure at each depth based on the pressure
    # at the previous depth. 
    # Note that P[0] is already set to 0 and doesn't need to be recalculated
    for i in range(1,len(z)):
        P[i] = P[i-1] + rho[i] * g * (z[i]-z[i-1])
    
    return P

