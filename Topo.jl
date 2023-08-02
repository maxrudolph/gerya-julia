using Printf

function initial_ice_depth(x::Float64)
    """
    Arguments:
    x::Float64 -- A horizontal value of float64

    Returns:
    The vertical value (depth) corresponding to the horizontal argument value

    Info:
    The initial topography setup follows the following.
    y = A*sin[k(x-b)]+c
    |A| is the amplitude
    k is the wave number (k=2π/wavelength)
    The period is 2π/k
    The horizontal shift is b
    The vertical shift is c
    """
    w = W
    A = 0.1*H
    k = (2*pi)/w
    b = 0.0
    c = initial_hice + 0.5e4
    return (A)*(sin(k*(x-b)))+(c)
end

function initial_surface_depth(x::Float64)
    """
    Arguments:
    x::Float64 -- A horizontal value of float64

    Returns:
    The vertical value (depth) corresponding to the horizontal argument value
    """
     return 0.5e4
end

function get_interface(grid::CartesianGrid,mat,contour_value,time_type::String)
    if time_type == "initial"
        i_interface_position = zeros(Float64,grid.nx+1);
        for j in 1:grid.nx+1
            i = 1
            while i <= grid.ny
                if i_mat[i,j] == contour_value
                    i_interface_position[j] = grid.yc[j]
                    break
                elseif i_mat[i+1,j] < contour_value
                    # interface is located within this cell.
                    i_interface_position[j] = grid.yc[i] + (grid.yc[i+1]-grid.yc[i])/(i_mat[i+1,j]-i_mat[i,j])*(contour_value-i_mat[i,j])
                    break
                end
                i = i+1
            end
        end
        return i_interface_position
    elseif time_type == "after"
        interface_position = zeros(Float64,grid.nx+1);
        for j in 1:grid.nx+1
        i = 1
            while i <= grid.ny
                if mat[i,j] == contour_value
                    interface_position[j] = grid.yc[i]
                    break
                elseif mat[i+1,j] < contour_value
                    # interface is located within this cell.
                    interface_position[j] = grid.yc[i] + (grid.yc[i+1]-grid.yc[i])/(mat[i+1,j]-mat[i,j])*(contour_value-mat[i,j])
                    break
                end
                i = i+1
            end
        end
        return interface_position
    end
end

function get_ice_thickness(interface_1::Vector{Float64},interface_2::Vector{Float64},time_type::String)
    if time_type == "initial"
        im_ice = maximum(i_ocean_ice_interface.-i_air_ice_interface)
        i_ice_ave = mean(i_ocean_ice_interface.-i_air_ice_interface)
        ix_ice = i_ocean_ice_interface.-i_air_ice_interface
        return im_ice,i_ice_ave,ix_ice
    elseif time_type == "after"
        m_ice = maximum(ocean_ice_interface.-air_ice_interface)
        ice_ave = mean(ocean_ice_interface.-air_ice_interface)
        x_ice = ocean_ice_interface.-air_ice_interface
        return m_ice,ice_ave,x_ice
    end
end 

# function get_amplitude(it_data::Tuple,interface_1::Vector{Float64},interface_2::Vector{Float64},time_type::String)
#     if time_type == "initial"
#         i_amp = maximum(i_ocean_ice_interface-i_air_ice_interface)-i_it_data[2]
#         return i_amp
#     elseif time_type == "after"
#         amp = maximum(ocean_ice_interface-air_ice_interface)-it_data[2]
#         return amp
#     end
# end

function get_amplitude(interface_1::Vector{Float64},interface_2::Vector{Float64},time_type::String)
    i_it_data = get_ice_thickness(i_air_ice_interface,i_ocean_ice_interface,"initial")
    it_data = get_ice_thickness(i_air_ice_interface,i_ocean_ice_interface,"after");
    if time_type == "initial"
        i_amp = maximum(i_ocean_ice_interface-i_air_ice_interface)-i_it_data[2]
        return i_amp
    elseif time_type == "after"
        amp = maximum(ocean_ice_interface-air_ice_interface)-it_data[2]
        return amp
    end
end

# function get_topo_info(x_time::String,it_data::Tuple,amp::Float64,unit_type::String,time_type::String)
#     if time_type == "initial" && unit_type == "km"
#         infot = printstyled("Intial:\n"; color = :blue, blink = true)
#         info1 = println("\tThe maximum total initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[1]/1000),"(km)")
#         info2 = println("\tThe average initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[2]/1000),"(km)")
#         info3 = println("\tThe initial amplitude is ",@sprintf("%.3g",i_amp/1000),"(km)")
#         return infot,info1,info2,info3
#     elseif time_type == "initial" && unit_type == "m"
#         infot = printstyled("Intial:\n"; color = :blue, blink = true)
#         info1 = println("\tThe maximum total initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[1]),"(km)")
#         info2 = println("\tThe average initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[2]),"(m)")
#         info3 = println("\tThe initial amplitude is ",@sprintf("%.3g",i_amp),"(m)")
#         return infot,info1,info2,info3
#     elseif time_type == "after" && unit_type == "km"
#         infot = printstyled("After:\n"; color = :blue, blink = true)
#         info1 = println("\tThe maximum total thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[1]/1000),"(km)")
#         info2 = println("\tThe average thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[2]/1000),"(km)")
#         info3 = println("\tThe amplitude after $x_time kyr is ",@sprintf("%.3g",amp/1000),"(km)")
#         return infot,info1,info2,info3
#     elseif time_type == "after" && unit_type == "m"
#         infot = printstyled("After:\n"; color = :blue, blink = true)
#         info1 = println("\tThe maximum total thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[1]),"(km)")
#         info2 = println("\tThe average thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[2]),"(m)")
#         info3 = println("\tThe amplitude after $x_time kyr is ",@sprintf("%.3g",amp),"(m)")
#         return infot,info1,info2,info3
#     end
# end

function get_topo_info(time::Float64,unit_type::String,time_type::String)
    x_time = @sprintf("%.3g",time/3.15e7/1e3)
    i_it_data = get_ice_thickness(i_air_ice_interface,i_ocean_ice_interface,"initial")
    i_amp = get_amplitude(i_air_ice_interface,i_ocean_ice_interface,"initial")
    it_data = get_ice_thickness(air_ice_interface,ocean_ice_interface,"after")
    amp = get_amplitude(air_ice_interface,ocean_ice_interface,"after")
    if time_type == "initial" && unit_type == "km"
        infot = printstyled("Intial:\n"; color = :blue, blink = true)
        info1 = println("\tThe maximum total initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[1]/1000),"(km)")
        info2 = println("\tThe average initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[2]/1000),"(km)")
        info3 = println("\tThe initial amplitude is ",@sprintf("%.3g",i_amp/1000),"(km)")
        return infot,info1,info2,info3
    elseif time_type == "initial" && unit_type == "m"
        infot = printstyled("Intial:\n"; color = :blue, blink = true)
        info1 = println("\tThe maximum total initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[1]),"(km)")
        info2 = println("\tThe average initial thickness of the icy shell is ",@sprintf("%.3g",i_it_data[2]),"(m)")
        info3 = println("\tThe initial amplitude is ",@sprintf("%.3g",i_amp),"(m)")
        return infot,info1,info2,info3
    elseif time_type == "after" && unit_type == "km"
        infot = printstyled("After:\n"; color = :blue, blink = true)
        info1 = println("\tThe maximum total thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[1]/1000),"(km)")
        info2 = println("\tThe average thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[2]/1000),"(km)")
        info3 = println("\tThe amplitude after $x_time kyr is ",@sprintf("%.3g",amp/1000),"(km)")
        return infot,info1,info2,info3
    elseif time_type == "after" && unit_type == "m"
        infot = printstyled("After:\n"; color = :blue, blink = true)
        info1 = println("\tThe maximum total thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[1]),"(km)")
        info2 = println("\tThe average thickness of the icy shell after $x_time kyr is ",@sprintf("%.3g",it_data[2]),"(m)")
        info3 = println("\tThe amplitude after $x_time kyr is ",@sprintf("%.3g",amp),"(m)")
        return infot,info1,info2,info3
    end
end


# function get_time_viscous(materials::Materials,it_data::Tuple,gy::Float64)
#     """
#     Arguments:
#     materials::Material --- Dyictionary in which holds the material properties used in model.
#     ice_ave::Float64 --- Ice shell thickness in (m)
#     gy::Float64 -- Acceleration due to gravity in (m/s^2)
    
#     Returns:
#     time --- time for viscous flow in yrs
#     """
#     time = (4*pi)*(materials.eta[2])*(1/gy)*(1/(materials.rho0[1]-materials.rho0[2]))*(1/it_data[2])
#     return time/3.15e7
# end

# function get_time_thickening(materials::Materials,detla_T::Float64,L::Float64,it_data::Tuple)
#     """
#     Arguments:
#     materials::Material --- Dyictionary in which holds the material properties used in model.
#     detla_T::Float64 --- A temperature value in (Kelvin) of float64, which is the temperature difference in
#                          which you material conducts through.
#     L::Float64 --- Latent heat of fusion in (J/kg)
#     ice_ave::Float64 --- Ice shell thickness in (m)
    
#     Returns:
#     time --- time for thickening in yrs

#     Info:
#     dh_dt = rate of thickening, units are m/s
#     1W = 1J/s
#     q = -k*(T-To)/h
#     """
#     dh_dt = materials.kThermal[2]*detla_T*(1/it_data[2])*(1/materials.rho0[1])*(1/L)
#     time = it_data[2]/dh_dt
#     return time/3.15e7
# end

function time_viscous(lambda::Float64)
    """
    Arguments:
    lambda::Float64 --- Topography Wavelength (m)
    
    Returns:
    time --- Time for viscous flow (yrs)
    """
    g = 0.113
    ice_viscosity = 10^14
    delta_rho = 80
    time = (4*pi)*(ice_viscosity)*(1/g)*(1/delta_rho)*(1/lambda)
    return time/3.15e7
end

function thickening_rate(hice::Float64)
    """
    Arguments:
    hice::Float64 --- Ice shell thickness (m)
    
    Returns:
    rate --- Rate of thickening, units are (m/s)

    Info:
    1W = 1J/s
    q = -k*(T-To)/h
    """
    delta_T = 173 # Temperature (K)
    kthermal = 2.2 # Thermal conductivity (W/m*K)
    L = 334*10^3 # Latent heat of fusion (J/kg)
    rho = 1000 # Density (kg/m^3)
    q = kthermal*(delta_T/hice)
    rate = q/L/rho   
    return rate
end

function thickening_time(hice::Float64)
    """
    Arguments:
    hice::Float64 --- Ice shell thickness (m)
    
    Returns:
    time --- Time for thickening (yrs)

    Info:
    1W = 1J/s
    q = -k*(T-To)/h
    """
    rate = thickening_rate(hice)
    time = hice/rate
    return time/3.15e7
end

function numerical_time_viscous(time::Float64)
    """
    Arguments:
    time::Float64 --- model time run (seconds)
    
    Returns:
    tr --- Characteristic time (kyr)

    Info:
    using equation (6.104) from Turcotte and Schubert
    w = w0*e^(-t/tr)
    """   
    i_amp = get_amplitude(i_air_ice_interface,i_ocean_ice_interface,"initial")
    amp = get_amplitude(i_air_ice_interface,i_ocean_ice_interface,"after")
    t = time/3.15e7/1e3
    tr = -(t)/log(amp/i_amp)
    return tr
end