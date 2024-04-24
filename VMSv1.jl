"""
Viscous flow vs. Melting and Soldification model
Note: 
Model does not allow for the phase change of ice to water by solving the energy equation in terms of entropy.
This model includes sticky air.
"""

if length(ARGS) > 9
    error("specify proper input arguments for range function for ice shell thickness and wavelength topogaphy for ocean-ice interface")
else
    ice_start = parse(Float64,ARGS[1])
    ice_stop = parse(Float64,ARGS[2])
    ice_length = parse(Int64,ARGS[3])
    wavelength_start = parse(Float64,ARGS[4])
    wavelength_stop = parse(Float64,ARGS[5])
    wavelength_length = parse(Int64,ARGS[6])
    percent_amplitude_start = parse(Float64,ARGS[7])
    percent_amplitude_stop = parse(Float64,ARGS[8])
    percent_amplitude_length = parse(Int64,ARGS[9])
    println("Ice Shell Thickness Range(start:difference between two subsequent objects:stop): ",range(ice_start,ice_stop,ice_length))
    println("Wavelength Topography for Ocean-Ice Interface Range(start:difference between two subsequent objects:stop): ",range(wavelength_start,wavelength_stop,wavelength_length))
    println("Amplitude Percentage Range(start:difference between two subsequent objects:stop): ",range(percent_amplitude_start,percent_amplitude_stop,percent_amplitude_length))
end

# Importing (using/include) packages and files needed for the code to run
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using WriteVTK
using Printf
using Statistics
using Dates
# using EasyFit
using Printf
using HDF5
include("Grid.jl")
include("Markers.jl")
include("Stokes.jl")
include("Temperature.jl")
include("TemperatureEntropy.jl")
include("GridOperations.jl")

# Note: that we import pyplot last to avoid a name conflict with grid
using PyPlot
include("Visualization.jl")

function initial_ice_depth(x::Float64,options::Dict)
    return options["ice thickness"] + options["surface depth"] + options["amplitude"]*cos( pi/options["wavelength"]*x )
end

function ice_viscosity(T::Float64)
    Q = 40000.0 # Activation Enegry (J/mol)
    R_cont = 8.314 # Gas Constant (J/molK)
    meltingpoint_viscosity = 1e15
    ice_vis = meltingpoint_viscosity*exp((Q*(273.0-T))/(R_cont*(273.0*T)))
    upperlimit = 1e25
    lowerlimit = meltingpoint_viscosity
    if ice_vis < lowerlimit
        ice_vis = lowerlimit
    elseif ice_vis > upperlimit
        ice_vis = upperlimit
    else
        ice_vis = ice_vis
    end
    return ice_vis
end

struct Materials
    # 1 - subsurface global ocean
    # 2 - icy shell
    # 3 - sticky air
    alpha::Vector{Float64} # Thermal expansion (1/K)
    rho0::Vector{Float64} # Density (kg/m^3)
    Hr::Vector{Float64} # Radiogenic heat production (W/m^3)
    Cp::Vector{Float64} # Heat capacity (J/kg*K)
    kThermal::Vector{Float64} # Thermal conductivity (W/m*K)
    eta::Vector{Float64} # Viscosity (Pa*s)
    function Materials()
        new([0.0,0.0,0.0],[1000.0,920.0,1.0],[0.0,0.0,0.0],[4180.0,2100.0,1.0e6],[0.5610,2.2,0.024],[1e12,1e15,1e17])
    end
end

function update_marker_prop!(markers::Markers,materials::Materials)
    eta = markers.scalarFields["eta"]
    rho = markers.scalarFields["rho"]
    T = markers.scalarFields["T"]
    mmat = markers.integers[markers.integerFields["material"],:]
    for i in 1:markers.nmark
        markers.scalars[rho,i] = materials.rho0[mmat[i]]
        if mmat[i] == 2
            markers.scalars[eta,i] = ice_viscosity(markers.scalars[T,i])
        end
    end
end

function update_marker_temp!(markers::Markers,materials::Materials)
    T = markers.scalarFields["T"]
    mmat = markers.integers[markers.integerFields["material"],:]
    for i in 1:markers.nmark
        if mmat[i] == 1
            markers.scalars[T,i] = 273.0
        elseif mmat[i] == 3
            markers.scalars[T,i] = 100.0
        end
    end
end

function initial_conditions!(markers::Markers,materials::Materials,options::Dict)
    material = markers.integerFields["material"]
    T = markers.scalarFields["T"]
    eta = markers.scalarFields["eta"]
    alpha = markers.scalarFields["alpha"]
    Cp = markers.scalarFields["Cp"]
    Hr = markers.scalarFields["Hr"]
    kThermal = markers.scalarFields["kThermal"]
    for i in 1:markers.nmark
        mx = markers.x[1,i]
        my = markers.x[2,i]
        hice = initial_ice_depth(mx,options)
        hsurf = options["surface depth"]
        if my > hice
            # subsurface global ocean
            markers.integers[material,i] = 1
            markers.scalars[T,i] = 273.0
            markers.scalars[eta,i] = materials.eta[1]
            markers.scalars[alpha,i] = materials.alpha[1]
            markers.scalars[Cp,i] = materials.Cp[1]
            markers.scalars[Hr,i] = materials.Hr[1]
            markers.scalars[kThermal,i] = materials.kThermal[1]
        elseif my > hsurf
            # icy shell
            markers.integers[material,i] = 2
            markers.scalars[T,i] = 100.0+((273.0-100.0)/(hice-hsurf))*(my-hsurf)
            markers.scalars[alpha,i] = materials.alpha[2]
            markers.scalars[Cp,i] = materials.Cp[2]
            markers.scalars[Hr,i] = materials.Hr[2]
            markers.scalars[kThermal,i] = materials.kThermal[2]
        else
            # sticky air
            markers.integers[material,i] = 3
            markers.scalars[T,i] = 100.0
            markers.scalars[eta,i] = materials.eta[3]
            markers.scalars[alpha,i] = materials.alpha[3]
            markers.scalars[Cp,i] = materials.Cp[3]
            markers.scalars[Hr,i] = materials.Hr[3]
            markers.scalars[kThermal,i] = materials.kThermal[3]
        end
    end
    # end loop over markers
    update_marker_prop!(markers,materials)
end

function get_interface(grid::CartesianGrid,mat::Matrix{Float64},contour_value::Float64)
    # Finding interfaces
    interface_position = zeros(Float64,grid.nx+1);
    for j in 1:grid.nx+1
        i = 1
        while i <= grid.ny
            if mat[i,j] == contour_value
                interface_position[j] = grid.yc[i]
                break
            elseif mat[i+1,j] < contour_value
                # interface is located within this cell.
                interface_position[j] = grid.yc[i] + ((grid.yc[i+1]-grid.yc[i])/(mat[i+1,j]-mat[i,j]))*(contour_value-mat[i,j])
                break
            end
            i = i+1
        end
    end
    return interface_position
end

function run(options::Dict,io)
    W = options["wavelength"]
    H = options["ice thickness"] + options["surface depth"] + options["amplitude"] + options["ice thickness"]/2
    ny = 100 #50
    nx = 111 #Int64(ceil(W/H*ny))
    gx = 0.0
    gy = 0.113

    # -1 = insulating, 1 = constant temp 
    Tbctype = [-1,-1,1,1] #left, right, top, bottom
    Tbcval = [0.0,0.0,100.0,273.0] #left, right, top, bottom
    bc = BoundaryConditions(0,0,0,0) # currently does nothing but is required argument to stokes solver.
    materials = Materials()
    markx = 6
    marky = 6
    seconds_in_year = 3.15e7
    plot_interval = 1e6*seconds_in_year # plot interval in seconds
    end_time = 3e7*seconds_in_year
    dtmax = plot_interval
    grid = CartesianGrid(W,H,nx,ny)
    println(io,"Grid resolution(ny x nx) : $ny x $nx")
    println(io,"Cell size in the x-direction is $(grid.W/(grid.nx-1))")
    println(io,"Cell size in the y-direction is $(grid.H/(grid.ny-1))")

    materials = Materials()
    println(io,"Creating Markers...")
    markers = Markers(grid,["alpha","T","rho","eta","Cp","Hr","kThermal"],["material"] ; nmx=markx,nmy=marky,random=true)
    println(io,"Initial condition...")
    initial_conditions!(markers, materials,options)

    ### Setting up agruments for interface function ###
    Ai = options["amplitude"] # initial amplitude

    ### Setting up agruments for termination criteria ###
    max_step::Int64=-1
    max_time::Float64=-1.0
    max_time = max_time == -1.0 ? typemax(Float64) : max_time
    max_step = max_step == -1 ? typemax(Int64) : max_step

    time = 0.0
    iout= 0
    last_plot = 0.0
    dt = dtmax

    rho_c = nothing
    rho_vx = nothing
    rho_vy = nothing
    alpha = nothing
    Hr = nothing
    Cp_c = nothing
    eta_s = nothing
    eta_n = nothing
    vxc = nothing
    vyc = nothing
    T = nothing
    dTmax = nothing
    dTemp = nothing
    Tnew = nothing
    Tlast = nothing
    kThermal = nothing
    mat = nothing
    Af = nothing

    itime = 1
    output_dir = options["visualization file path"]*"/Data"

    # touch(options["text file path"]*"/Time_Data.txt")
    # # Writing to a file using write() method
    # start_data = open(options["text file path"]*"/Time_Data.txt","w")
    # write(start_data,"Table:\n")
    # write(start_data,"Amplitude\t| Thickening Time(Myr)\n")
    # # Closing the file in order to write the content from the disk to file
    # close(start_data)

    terminate = false
    while !terminate
        ### update the markers properties ###
        update_marker_prop!(markers,materials)
        update_marker_temp!(markers,materials)
        ### Transfer properties markers -> nodes ###
        # Basic Nodes
        eta_s_new, = marker_to_stag(markers,grid,["eta",],"basic",method="logarithmic")
        # Cell Centers
        rho_c_new,Cp_c_new,alpha_new,Hr_new,kThermal_new = marker_to_stag(markers,grid,["rho","Cp","alpha","Hr","kThermal"],"center")
        eta_n_new, = marker_to_stag(markers,grid,["eta",],"center",method="logarithmic")
        # interpolate temperature using rhocp as weight
        rhocp = markers.scalars[markers.scalarFields["rho"],:] .* markers.scalars[markers.scalarFields["Cp"],:]
        Tlast_new, = marker_to_stag(markers,grid,["T"],"center",extra_weight = rhocp)
        # Vx and Vy nodes:
        rho_vx_new, = marker_to_stag(markers,grid,["rho",],"vx")
        rho_vy_new, = marker_to_stag(markers,grid,["rho",],"vy")

        # deal with any NaN values from interpolation:
        if itime > 1
            if any(isnan.(eta_s_new))
                println("found nan values")
            end
            replace_nan!(eta_s,eta_s_new)
            replace_nan!(rho_c,rho_c_new)
            replace_nan!(Hr,Hr_new)
            replace_nan!(Cp_c,Cp_c_new)
            replace_nan!(alpha,alpha_new)
            replace_nan!(eta_n,eta_n_new)
            replace_nan!(Tlast,Tlast_new)
            replace_nan!(rho_vx,rho_vx_new)
            replace_nan!(rho_vy,rho_vy_new)
            replace_nan!(kThermal,kThermal_new)
        end
        # Copy field data
        kThermal = copy(kThermal_new)
        rho_vx = copy(rho_vx_new)
        rho_vy = copy(rho_vy_new)
        rho_c = copy(rho_c_new)
        Hr = copy(Hr_new)
        Cp_c = copy(Cp_c_new)
        alpha = copy(alpha_new)
        eta_s = copy(eta_s_new)
        eta_n = copy(eta_n_new)
        Tlast = copy(Tlast_new)

        if itime == 1
            ghost_temperature_center(grid,Tlast,Tbctype,Tbcval)
            #cell_center_to_markers!(markers,grid,Tlast,markers.scalars[[markers.scalarFields["T"],],:])
        else
            ghost_temperature_center(grid,Tlast,Tbctype,Tbcval)
        end

        # Assembling and solving the stokes equations
        L,R = form_stokes(grid,eta_s,eta_n,rho_vx,rho_vy,bc,gx,gy;dt=0.0)
        stokes_solution = L\R
        vx,vy,P = unpack(stokes_solution,grid;ghost=true)

        # Obtaining velocity at the cell centers
        vxc,vyc = velocity_to_centers(grid,vx,vy)
        adiabatic_heating = compute_adiabatic_heating(grid,rho_c,Tlast,alpha,gx,gy,vxc,vyc)*0.0
        shear_heating = compute_shear_heating(grid,vx,vy,eta_n,eta_s)*0.0
        H = (adiabatic_heating .+ shear_heating .+ Hr).*0.0

        diffusion_timestep = calculate_diffusion_timestep(grid,options)
        if itime > 1
            this_dtmax = min(1.2*dt,dtmax)
        else
            this_dtmax = dtmax
        end
        if this_dtmax > diffusion_timestep
            dt = diffusion_timestep
        else 
            dt = this_dtmax
        end
        dt = compute_timestep(grid,vxc,vyc;dtmax=this_dtmax,cfl=0.25)
        if dt > diffusion_timestep
            dt = diffusion_timestep
        end

        if itime == 1   
            println(io,"Diffusion timestep is ",diffusion_timestep/seconds_in_year," yr, ",diffusion_timestep/seconds_in_year/1e3," Kyr, ",diffusion_timestep/seconds_in_year/1e6," Myr")
        end
        println(io,"Starting step $itime, with dt = ",dt/seconds_in_year," yr, ",dt/seconds_in_year/1e3," Kyr, ",dt/seconds_in_year/1e6," Myr")

        dTmax = Inf
        dTemp = nothing
        Tnew = nothing
        titer = 1
        for titer=1:2
            # println("Trying with timestep ",dt/seconds_in_year/1e3," kyr")
            L,R = assemble_energy_equation_center(grid,rho_c,Cp_c,kThermal,H,Tlast,dt,Tbctype,Tbcval)
            Tnew = L\R;
            Tnew = reshape(Tnew,grid.ny,grid.nx);
            Tnew = ghost_temperature_center(grid,Tnew,Tbctype,Tbcval)
            T = copy(Tnew)
            dTemp = Tnew-Tlast
            # compute the maximum temperature change
            dTmax = maximum(abs.(dTemp[2:end-1,2:end-1]))
            println(io,"dTmax=",dTmax," dt=",dt/seconds_in_year/1e3," kyr")
            dt = min(dt,dTmax < 10.0 ? dt : dt*10.0/dTmax)
            if dTmax < 10.0
                break
            end
        end

        dT_subgrid_node = subgrid_temperature_relaxation_center!(markers,grid,Tlast,dt)
        dT_remaining = dTemp - dT_subgrid_node
        cell_center_change_to_markers!(markers,grid,dT_remaining,"T")

        ### Setting up agruments for Termination Criteria ###
        mat, = marker_to_stag(markers,grid,markers.integers[[markers.integerFields["material"]],:],"center")
        ocean_ice_interface = get_interface(grid,mat,1.5)
        air_ice_interface = get_interface(grid,mat,2.5)
        max_ice_shell_thickness = maximum(ocean_ice_interface)-maximum(air_ice_interface)
        avg_ice_shell_thickness = mean(ocean_ice_interface)-mean(air_ice_interface)
        Af = max_ice_shell_thickness-avg_ice_shell_thickness
        i_A = @sprintf("%.6g",Ai/1e3)
        f_A = @sprintf("%.6g",Af/1e3)
        println(io,"Initial Amplitude: $i_A (km), Final amplitude: $f_A (km)")

        # if mod(itime,20) == 0
        #     rate = get_thickening_rate(avg_ice_shell_thickness)
        #     t_tic = get_thickening_time(avg_ice_shell_thickness,rate)
        #     i_time = @sprintf("%4.6g",t_tic/1e6)
        #     # Open file in append mode
        #     start_data = open(options["text file path"]*"/Time_Data.txt","a")
        #     write(start_data,"$f_A\t\t\t\t")
        #     write(start_data,"$i_time\n")
        #     # Closing the file in order to write the content from the disk to file
        #     close(start_data)
        # end

        # Checking Termination Criteria, time is in Myr, amplitude is in meters
        if time >= max_time || itime >= max_step || Af/Ai <= 1/exp(1)
            terminate = true
        end

        # if time == 0.0 || mod(itime,100) == 0 || terminate
        #     last_plot = time 
        #     # Gird output
        #     name1 = @sprintf("%s/viz.%04d.vtr",output_dir,iout)
        #     println("Writing visualization file = ",name1)
        #     vn = velocity_to_basic_nodes(grid,vxc,vyc)
        #     visualization(grid,rho_c,eta_s,vn,P,Tnew[2:end-1,2:end-1],time/seconds_in_year;filename=name1)
        #     # Markers output
        #     name2 = @sprintf("%s/markers.%04d.vtp",output_dir,iout)
        #     println("Writing visualization file = ",name2)
        #     visualization(markers,time/seconds_in_year;filename=name2)
        #     # # Hydrostatic Pressure output
        #     # name3 = @sprintf("%s/hp.%04d.vtr",output_dir,iout)
        #     # println("Writing visualization file = ",name3)
        #     # output_fields = Dict("Hydrostatic Pressure"=>hp[2:end-1,2:end-1])
        #     # visualization(grid,output_fields,time/seconds_in_year;filename=name3)
        #     iout += 1
        # end

        # println(io,"Min/Max velocity: ",minimum(vyc)," ",maximum(vyc))
        # Moving the markers and advancing to the next timestep
        move_markers_rk4!(markers,grid,vx,vy,dt,continuity_weight=1/3)
        time += dt
        println(io,"Finished step $itime")
        println(io,"time = ",time/seconds_in_year," yr, ",time/seconds_in_year/1e3," Kyr, ",time/seconds_in_year/1e6," Myr") 
        itime += 1
        # mat, = marker_to_stag(markers,grid,markers.integers[[markers.integerFields["material"]],:],"center")
        # ocean_ice_interface = get_interface(grid,mat,1.5)
        # air_ice_interface = get_interface(grid,mat,2.5)
        # append!(topography,[ocean_ice_interface])
        # append!(amplitude,Af)
    end
    return grid,time,itime,Af
    #return grid,i_mat,mat,time,itime
    # return grid,i_mat,mat,time_plot,time,itime,amplitude
end

include("Topo.jl")
include("Outputs.jl")

function model_run()
    options = Dict()
    top_dir = mk_modelrun_dir()
    nlambda = wavelength_length
    nhice = ice_length
    namp = percent_amplitude_length
    local lambda = range(wavelength_start,wavelength_stop,nlambda)
    local hice =  range(ice_start,ice_stop,nhice)
    local per_amp = range(percent_amplitude_start,percent_amplitude_stop,namp)
    options["specific heat of ice"] = 2.1e3 # J/kg*K (ice)
    options["density of ice"] = 920.0 # kg/m^3
    options["thermal conductivity of ice"] = 2.2 # W/m*K
    options["thermal diffusivity"] = options["thermal conductivity of ice"] / (options["density of ice"]*options["specific heat of ice"]) # m^2/s
    irun = 1
    for k in 1:namp
        t_halfspace = zeros(nlambda,nhice)
        t_rel = zeros(nlambda,nhice)
        t_tic = zeros(nlambda,nhice)
        amplitude_percentage = per_amp[k]
        amp_decimal = amplitude_percentage/100
        sub_dir = mk_sub_dir(top_dir,amplitude_percentage)
        for i in 1:nlambda
            for j in 1:nhice
                options["wavelength"] = lambda[i]*1e3
                options["ice thickness"] = hice[j]*1e3
                options["amplitude"] = amp_decimal*options["ice thickness"]
                options["surface depth"] = options["amplitude"]
                println("Starting model execution for model run $irun...")
                sub_dir_by_run,sub_dir_plots,sub_dir_data = mk_output_dir(sub_dir,irun)
                options["visualization file path"] = sub_dir_by_run
                options["text file path"] = sub_dir_by_run
                data_info(ice_start,ice_stop,nhice,wavelength_start,wavelength_stop,nlambda,sub_dir,amplitude_percentage)
                println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)"," , ","Using Amplitde Percentage: $amplitude_percentage%")
                io = open(sub_dir_by_run*"/output.txt","w")
                # model_runtime(sub_dir_by_run,"Start")
                grid,time,itime,Af = run(options,io);
                # model_runtime(sub_dir_by_run,"End")
                t_rel[i,j] = get_numerical_time_viscous(options["amplitude"],Af,time)
                t_halfspace[i,j] = get_halfspace_time_viscous(options["wavelength"])
                rate = get_thickening_rate(options["ice thickness"])
                t_tic[i,j] = get_thickening_time(options["amplitude"],rate)
                println(io,"Analytic relaxation time: ",t_halfspace[i,j],"(yr)",t_halfspace[i,j]/1e3,"(kyr) or ",t_halfspace[i,j]/1e6,"(Myr)")
                println(io,"Numerical relaxation time: ",t_rel[i,j],"(yr)",t_rel[i,j]/1e3,"(kyr) or ",t_rel[i,j]/1e6,"(Myr)")
                println(io,"Thickening time: ",t_tic[i,j]/1e3,"(kyr) or ",t_tic[i,j]/1e6,"(Myr)")
                close(io)
                println("Model ran successfully for model run $irun. Outputs saved to output.txt")
                irun += 1
            end
        end
        # return lambda,hice,per_amp,t_halfspace,t_rel,t_tic,sub_dir
        data_to_hdf5_file(lambda,hice,t_halfspace,t_rel,t_tic,sub_dir)
        if namp != 1
            irun = 1
        end
    end
end

# function thickening_run()
#     options = Dict()
#     top_dir = mk_modelrun_dir()
#     nlambda = wavelength_length
#     nhice = ice_length
#     namp = percent_amplitude_length
#     local lambda = range(wavelength_start,wavelength_stop,nlambda)
#     local hice =  range(ice_start,ice_stop,nhice)
#     local per_amp = range(percent_amplitude_start,percent_amplitude_stop,namp)
#     irun = 1
#     for k in 1:namp
#         t_tic = zeros(nlambda,nhice)
#         amplitude_percentage = per_amp[k]
#         amp_decimal = amplitude_percentage/100
#         sub_dir = mk_sub_dir(top_dir,amplitude_percentage)
#         for i in 1:nlambda
#             for j in 1:nhice
#                 options["wavelength"] = lambda[i]*1e3
#                 options["ice thickness"] = hice[j]*1e3
#                 options["amplitude"] = amp_decimal*options["ice thickness"]
#                 options["surface depth"] = options["amplitude"] 
#                 println("Starting model execution for model run $irun...")
#                 sub_dir_by_run,sub_dir_plots,sub_dir_data = mk_output_dir(sub_dir,irun)
#                 options["visualization file path"] = sub_dir_by_run
#                 options["text file path"] = sub_dir_by_run
#                 data_info(ice_start,ice_stop,nhice,wavelength_start,wavelength_stop,nlambda,sub_dir,amplitude_percentage)
#                 println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)"," , ","Using Amplitde Percentage: $amplitude_percentage%")
#                 open(sub_dir_by_run*"/output.txt", "w") do out
#                     redirect_stdout(out) do 
#                         model_runtime(sub_dir_by_run,"Start")
#                         rate = get_thickening_rate(options["ice thickness"])
#                         t_tic[i,j] = get_thickening_time(options["amplitude"],rate)
#                         model_runtime(sub_dir_by_run,"End")
#                         println("Thickening time :",t_tic[i,j]/1e3,"(kyr) or ",t_tic[i,j]/1e6,"(Myr)")
#                     end
#                 end
#                 println("Model ran successfully for model run $irun. Outputs saved to output.txt")
#                 irun += 1
#             end
#         end
#         # thickening_data_to_hdf5_file(lambda,hice,t_tic,sub_dir)
#         if namp != 1
#             irun = 1
#         end
#     end
# end

try
    model_run()
catch e
    println("Model encountered an error. Error details saved to error_log.txt")
    open("error_log.txt","w")do out
        redirect_stdout(out) do
            showerror(stdout,e,catch_backtrace())
        end
    end
    @error "Something went wrong" exception=(e, catch_backtrace())
end
