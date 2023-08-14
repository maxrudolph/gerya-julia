using PyPlot
using Printf

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

function get_time_viscous(lambda::Float64)
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
    ice_viscosity = 10^15
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
    rate = q/L/rho   
    return rate
end

function get_thickening_time(hice::Float64)
    """
    Arguments:
    hice::Float64 --- Ice shell thickness in (m)
    
    Returns:
    time --- Time for thickening in (yr)

    Info:
    rate --- Rate of thickening in (m/s) 
    """
    rate = get_thickening_rate(hice)
    time = hice/rate
    return time/3.15e7
end

function get_numerical_time_viscous(i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},time::Float64)
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
    initial_max_ice_shell_thickness = maximum(i_interface_2.-i_interface_1)
    initial_avg_ice_shell_thickness = mean(i_interface_2.-i_interface_1)
    initial_amplitude = initial_max_ice_shell_thickness-initial_avg_ice_shell_thickness
    max_ice_shell_thickness = maximum(interface_2.-interface_1)
    avg_ice_shell_thickness = mean(interface_2.-interface_1)
    amplitude = max_ice_shell_thickness-avg_ice_shell_thickness
    t = time/3.15e7
    tr = -(t)/log(amplitude/initial_amplitude)
    return tr
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

function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},topography::Vector{Any},times::Vector{Any},time::Float64)
    x_time = @sprintf("%.3g",time/3.15e7/1e3)
    # Inital Model Schematic Profile
    figure() 
    pcolor(grid.xc/1000,grid.yc/1000,i_mat)
    colorbar(cmap="viridis")
    plot(grid.xc/1000,i_interface_1/1000,"m",label="air-ice interface")
    plot(grid.xc/1000,i_interface_2/1000,"r",label="ocean-ice interface")
    title(L"Initial\,\,Model\,\,Schematic")
    gca().invert_yaxis()
    gca().set_aspect("equal")
    gca().set_ylabel(L"Height\,(km)")
    gca().set_xlabel(L"Width\,(km)")
    # Legend is at the bottom
    legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5)
    show()

    # Initial ice shell thickness profile along the horizontal direction
    figure() 
    plot(grid.xc/1000,(i_interface_2.-i_interface_1)/1000)
    title(L"Initial\,\,Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction")
    gca().invert_yaxis()
    gca().set_ylabel(L"thickness\,(km)")
    gca().set_xlabel(L"x\,(km)")
    show()

    # Model Schematic Profile
    figure()
    pcolor(grid.xc/1000,grid.yc/1000,mat)
    colorbar(cmap="viridis")
    plot(grid.xc/1000,interface_1/1000,"m",label="air-ice interface")
    plot(grid.xc/1000,interface_2/1000,"r",label="ocean-ice interface")
    gca().invert_yaxis()
    gca().set_aspect("equal")
    gca().set_ylabel(L"Height\,(km)")
    gca().set_xlabel(L"Width\,(km)")
    title(L"Final\,\,Model\,\,Schematic")
    # Legend is at the bottom
    legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5)        
    show()

    # Ice shell thickness profile along the horizontal direction
    figure()
    plot(grid.xc/1000,(interface_2.-interface_1)./1000)
    title(L"Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction\,\,at\,\,%$x_time\,kyr")
    gca().invert_yaxis()
    gca().set_ylabel(L"thickness\,(km)")
    gca().set_xlabel(L"x\,(km)")
    show()

    # Profile of ice-water inferface topograpgy over time 
    figure()
    for i in 1:10:100
        plot(grid.xc,topography[i],label=(L"At",@sprintf("%.3g",times[i]/3.15e7/1e3),L"kyr"))
    end
    title(L"Profile\,\,of\,\,Ice-Water\,\,Interface\,\,Topography\,\,Over\,\,Time")
#     gca().set_xlim([0.0,1e4])
#     gca().set_ylim([1.8e4,2.2e4])
    gca().invert_yaxis()
    # Legend is at the bottom
    legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5)
    show()

#     # Profile of maximum topograpgy over time 
#     figure()
#     plot(times/3.15e7/1e3,max_topo/1000)
#     title(L"Maximum\,\,Topography\,\,After\,\,%$x_time\,ka")
#     gca().set_ylabel(L"Max.\,topography\,(km)")
#     gca().set_xlabel(L"Time\,(ka)")
#     show()
end 