# function initial_ice_depth(x::Float64,ice_thickness::Float64,wavelength::Float64,amplitude::Float64,initial_surface_depth::Float64)
#     """
#     Arguments:
#     x::Float64 -- A horizontal value of float64

#     Returns:
#     The vertical value (depth) corresponding to the horizontal argument value
    
#     Info:
#     The initial topography setup follows the following.
#     y = A*sin[k(x-b)]+c
#     |A| is the amplitude
#     k is the wave number (k=2π/wavelength)
#     The period is 2π/k
#     The horizontal shift is b
#     The vertical shift is c
#     """
#     return ice_thickness + initial_surface_depth + amplitude*sin( 2*pi/wavelength*x )
# end

# function initial_ice_depth(x::Float64,ice_thickness::Float64,wavelength::Float64,amplitude::Float64,initial_surface_depth::Float64)
#     return ice_thickness + initial_surface_depth + amplitude*sin( 2*pi/wavelength*x )
# end

# function get_interface(grid::CartesianGrid,mat::Matrix{Float64},contour_value::Float64)
#     # Finding interfaces 
#     interface_position = zeros(Float64,grid.nx+1);
#     for j in 1:grid.nx+1
#         i = 1
#         while i <= grid.ny
#             if mat[i,j] == contour_value
#                 interface_position[j] = grid.yc[i]
#                 break
#             elseif mat[i+1,j] < contour_value
#                 # interface is located within this cell.
#                 interface_position[j] = grid.yc[i] + (grid.yc[i+1]-grid.yc[i])/(mat[i+1,j]-mat[i,j])*(contour_value-mat[i,j])
#                 break
#             end
#             i = i+1
#         end
#     end
#     return interface_position
# end

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
    kthermal = 2.2
    L = 334*10^3 
    rho = 1000
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
    # initial_max_ice_shell_thickness = maximum(i_interface_2-i_interface_1)
    # initial_avg_ice_shell_thickness = mean(i_interface_2-i_interface_1)
    # initial_amplitude = initial_max_ice_shell_thickness-initial_avg_ice_shell_thickness
    # max_ice_shell_thickness = maximum(interface_2-interface_1)
    # avg_ice_shell_thickness = mean(interface_2-interface_1)
    # amplitude = max_ice_shell_thickness-avg_ice_shell_thickness
    lnratio = log2(final_amplitude)-log2(initial_amplitude)
    tr = -time/lnratio
    return tr/3.15e7
end

# function get_topography_info(interface_1::Vector{Float64},interface_2::Vector{Float64},time::Float64,time_type::String)
#     x_time = @sprintf("%.3g",time/3.15e7/1e3)
#     max_ice_shell_thickness = maximum(interface_2.-interface_1)
#     avg_ice_shell_thickness = mean(interface_2.-interface_1)
#     x_ice = interface_2.-interface_1
#     amplitude = max_ice_shell_thickness-avg_ice_shell_thickness
#     if time_type == "initial"
#         infot = printstyled("Intial:\n"; color = :blue, blink = true)
#         info1 = println("\tThe maximum total initial thickness of the icy shell is ",@sprintf("%.3g",max_ice_shell_thickness/1000),"(km)")
#         info2 = println("\tThe average initial thickness of the icy shell is ",@sprintf("%.3g",avg_ice_shell_thickness/1000),"(km)")
#         info3 = println("\tThe initial amplitude is ",@sprintf("%.3g",amplitude/1000),"(km)")
#         return infot,info1,info2,info3
#     elseif time_type == "final"
#         infot = printstyled("After:\n"; color = :blue, blink = true)
#         info1 = println("\tThe maximum total thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",max_ice_shell_thickness/1000),"(km)")
#         info2 = println("\tThe average thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",avg_ice_shell_thickness/1000),"(km)")
#         info3 = println("\tThe amplitude after $x_time kyr is ",@sprintf("%.3g",amplitude/1000),"(km)")
#         return infot,info1,info2,info3
#     end
# end 