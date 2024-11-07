function get_halfspace_time_viscous(options::Dict)
    """
    Arguments:
    lambda::Float64 --- Topography Wavelength in (m)
    
    Returns:
    time --- Time for viscous flow in (yr)
    
    Info:
    ***This is the halfspace solution***
    g --- Acceleration due to gravity in (m/s^2)
    ice_visosity --- Viscosity of ice in (Pa*s = kg/m*s)
    delta_rho --- Difference of density in (kg/m^3)
    """
    g = options["gravity of icy moon"]
    ice_viscosity = options["reference viscosity"]
    delta_rho = options["density of ocean"] - options["density of ice"]
    time = (4*pi)*(ice_viscosity)*(1/g)*(1/delta_rho)*(1/options["wavelength"])
    return time/3.15e7
end

function get_thickening_rate(options::Dict)
    """
    Arguments:
    hice::Float64 --- Ice shell thickness in (m)
    
    Returns:
    rate --- Rate of thickening in (m/s)

    Info:
    delta_T --- Difference of temperaure in (K)
    kthermal --- Thermal conductivity in (W/m*K)
    L --- Latent heat of fusion in (J/kg)
    rho --- Density in (kg/m^3)
    q --- Heat flux in (W/m^2)
    
    q = -k*(T-To)/h
    rate = q/rho*L
    *1W = 1J/s*
    """
    delta_T = 173
    kthermal = options["thermal conductivity of ice"]
    L = options["latent heat of fusion"]
    rho = options["density of ice"]
    q = kthermal*(delta_T/options["hice"])
    rate = q/(L*rho)   
    return rate
end

function get_thickening_time(options,rate::Float64)
    """
    Arguments:
    hice::Float64 --- Ice shell thickness in (m)
    
    Returns:
    time --- Time for thickening in (yr)

    Info:
    rate --- Rate of thickening in (m/s) 
    """
    time = options["amplitude"]/rate
    return time/3.15e7
end

function compute_numerical_thickening_time(h::Vector{Any},t::Vector{Any},options::Dict)
    thickening_rate = diff(h)./diff(t)  # Change in thickness / Change in time
    average_thickening_rate = mean(thickening_rate)
    t_thickening = options["amplitude"]/average_thickening_rate
    return t_thickening/3.15e7
end

function get_numerical_time_viscous(initial_amplitude::Float64,final_amplitude::Float64,time::Float64)
    """
    Arguments:
    time::Float64 --- model time run in (seconds)
    
    Returns:
    tr --- Characteristic time for viscous flow in (yr)

    Info:
    using equation (6.104) from Turcotte and Schubert
    w = w0*e^(-t/tr)
    
    i_amp --- Initial amplitude
    amp --- Final amplitude
    """   
    lnratio = log2(final_amplitude)-log2(initial_amplitude)
    tr = -time/lnratio
    return tr/3.15e7
end
