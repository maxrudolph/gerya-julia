function get_halfspace_time_viscous(lambda::Float64)
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
    g = 0.113 
    ice_viscosity = 1e15
    delta_rho = 80
    time = (4*pi)*(ice_viscosity)*(1/g)*(1/delta_rho)*(1/lambda)
    return time/3.15e7
end

function get_thickening_rate(hice::Float64)
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
    kthermal = 2.14
    L = 334*10^3 
    rho = 916.73
    q = kthermal*(delta_T/hice)
    rate = q/(L*rho)   
    return rate
end

function get_thickening_time(hice::Float64,rate::Float64)
    """
    Arguments:
    hice::Float64 --- Ice shell thickness in (m)
    
    Returns:
    time --- Time for thickening in (yr)

    Info:
    rate --- Rate of thickening in (m/s) 
    """
    time = hice/rate
    return time/3.15e7
end

function compute_numerical_thickening_time(h::Vector{Any},t::Vector{Any},hi::Float64)
    thickening_rate = diff(h)./diff(t)  # Change in thickness / Change in time
    average_thickening_rate = mean(thickening_rate)
    t_thickening = hi/average_thickening_rate
    return t_thickening/3.15e7
end

function get_numerical_time_viscous(initial_amplitude::Float64,final_amplitude::Float64,time::Float64)
# function get_numerical_time_viscous(i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},time::Float64)
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