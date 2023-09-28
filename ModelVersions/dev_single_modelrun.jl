# println("Input the arugments for a julia range function for ice shell thickness")
# print("Type in arugments separated by a space: ")
# userinput_ice = readline()
# userinput_ice = rsplit(userinput_ice," ")
# userinput_ice = map(x->parse(Float64,x),userinput_ice)
# ice_length = convert(Int64,userinput_ice[3])
# # range(start=userinput_ice[1],stop=userinput_ice[2],length=ice_length)
# if length(userinput_ice) != 3
#     error("specify proper input arguments for range function")
# else
#     println(range(start=userinput_ice[1],stop=userinput_ice[2],length=ice_length))
# end
# println("Input the arugments for a julia range function for wavelength topography")
# print("Type in arugments separated by a space: ")
# userinput_wavelength = readline()
# userinput_wavelength = rsplit(userinput_wavelength," ")
# userinput_wavelength = map(x->parse(Float64,x),userinput_wavelength)
# wavelength_length = convert(Int64,userinput_wavelength[3])
# if length(userinput_wavelength) != 3
#     error("specify proper input arguments for range function")
# else
#     println(range(start=userinput_wavelength[1],stop=userinput_wavelength[2],length=wavelength_length))
# end

# Importing (using/include) packages and files needed for the code to run
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using WriteVTK
using Printf
using Statistics 
using Dates
using EasyFit
using Printf
using HDF5
include("Grid.jl")
include("Markers.jl")
include("Stokes.jl")
include("Temperature.jl")
include("GridOperations.jl")

# Note: that we import pyplot last to avoid a name conflict with grid
using PyPlot
include("Visualization.jl")

function initial_ice_depth(x::Float64,ice_thickness::Float64,wavelength::Float64,amplitude::Float64,initial_surface_depth::Float64)
    return ice_thickness + initial_surface_depth + amplitude*sin( 2pi/wavelength*x )
end

function ice_viscosity(T::Float64)
    Q = 40000.0 # Activation Enegry (J/mol)
    R_cont = 8.314 # Gas Constant (J/molK)
    ice_vis = (1e15)*exp((Q*(273.0-T))/(R_cont*(273.0*T)))
    upperlimit = 1e25
    lowerlimit = 1e12
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
        new([0.0,0.0,0.0],[1000.0,920.0,1.0],[0.0,0.0,0.0],[4180.0,2100.0,1.0e6],[0.5610,2.1,0.024],[1e12,1e15,1e17])
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
        hice = initial_ice_depth(mx,options["ice thickness"],options["wavelength"],options["amplitude"],options["surface depth"])
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
            # markers.scalars[eta,i] = eta_i[i]
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
    interface_position = zeros(Float64,grid.nx-1);
    for j in 1:grid.nx-1
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

function run(options::Dict)
    W = options["wavelength"]
    H = options["ice thickness"] + options["surface depth"] + options["amplitude"] + 1e4
    ny = 251
    nx = Int64(ceil(W/H*ny))
    gx = 0.0
    gy = 0.113
    
    println("nx: $nx, ny: $ny")

    # Tbctype = [-1,-1,1,1] #left, right, top, bottom
    Tbctype = [-1,-1,1,1]
    # Tbcval = [0.0,0.0,100.0,273.0] #left, right, top, bottom
    Tbcval = [0.0,0.0,100.0,273.0]
    bc = BoundaryConditions(0,0,0,0) # currently does nothing but is required argument to stokes solver.
    materials = Materials()

    markx = 6
    marky = 6
    seconds_in_year = 3.15e7
    plot_interval = 1e6*seconds_in_year # plot interval in seconds
    end_time = 3e7*seconds_in_year
    dtmax = plot_interval
    grid = CartesianGrid(W,H,nx,ny)
    println("Creating Markers...")
    @time markers = Markers(grid,["alpha","T","rho","eta","Cp","Hr","kThermal"],["material"] ; nmx=markx,nmy=marky,random=true)
    println("Initial condition...")
    @time initial_conditions!(markers, materials,options)

    local time_plot = []
    # local topography = []
    local amplitude = []
    
    ### Setting up agruments for interface function ###
    # initial 
    i_mat, = marker_to_stag(markers,grid,markers.integers[[markers.integerFields["material"]],:],"center")
    i_air_ice_interface = get_interface(grid,i_mat,2.5)
    i_ocean_ice_interface = get_interface(grid,i_mat,1.5)
    Ai = options["amplitude"] 

    ### Setting up agruments for termination criteria ###
    max_step::Int64=-1
    max_time::Float64=-1.0
    max_time = max_time == -1.0 ? typemax(Float64) : max_time
    max_step = max_step == -1 ? typemax(Int64) : max_step
    
    time = 0.0
    iout= 0
    last_plot = 0.0
    dt = 1e10

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
    x_time = nothing
    kThermal = nothing
    # ocean_ice_interface = nothing
    mat = nothing

    itime = 1
    output_dir = "test"

    terminate = false
    while !terminate
        # 0. update the markers properties  
        update_marker_prop!(markers,materials)
        update_marker_temp!(markers,materials)
        # 1. Transfer properties markers -> nodes
        # 1a. Basic Nodes
        eta_s_new, = marker_to_stag(markers,grid,["eta",],"basic")
        # 1b. Cell Centers
        rho_c_new,Cp_c_new,alpha_new,eta_n_new,Tlast_new,Hr_new,kThermal_new = marker_to_stag(markers,grid,["rho","Cp","alpha","eta","T","Hr","kThermal"],"center")
        # 1c. Vx and Vy nodes:
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
            cell_center_to_markers!(markers,grid,Tlast,markers.scalars[[markers.scalarFields["T"],],:])
        else
            ghost_temperature_center(grid,Tlast,Tbctype,Tbcval)
        end

        # 2. Assemble and solve the stokes equations
        L,R = form_stokes(grid,eta_s,eta_n,rho_vx,rho_vy,bc,gx,gy;dt=dt)
        stokes_solution = L\R
        vx,vy,P = unpack(stokes_solution,grid;ghost=true)

        # Get the velocity at the cell centers:
        vxc,vyc = velocity_to_centers(grid,vx,vy)
        adiabatic_heating = compute_adiabatic_heating(grid,rho_c,Tlast,alpha,gx,gy,vxc,vyc)*0.0
        shear_heating = compute_shear_heating(grid,vx,vy,eta_n,eta_s)*0.0
        H = (adiabatic_heating .+ shear_heating .+ Hr).*0.0

        # 3. Compute the advection timestep
        if itime > 1
            this_dtmax = min(1.2*dt,dtmax)
        else
            this_dtmax = dtmax
        end
        dt = compute_timestep(grid,vxc,vyc;dtmax=this_dtmax,cfl=0.25)
        diffusion_timestep = (grid.x[2]-grid.x[1])^2 / 1e-6
        if dt > diffusion_timestep
            dt = diffusion_timestep
        end

        dTmax = Inf
        dTemp = nothing
        Tnew = nothing
        titer = 1
        for titer=1:2
            # assemble and solve the energy equation
            # println("Trying with timestep ",dt/seconds_in_year/1e3," kyr")
            L,R = assemble_energy_equation_center(grid,rho_c,Cp_c,kThermal,H,Tlast,dt,Tbctype,Tbcval)
            Tnew = L\R;
            Tnew = reshape(Tnew,grid.ny,grid.nx);
            Tnew = ghost_temperature_center(grid,Tnew,Tbctype,Tbcval)
            T = copy(Tnew)

            dTemp = Tnew-Tlast
            # compute the maximum temperature change
            dTmax = maximum(abs.(dTemp[2:end-1,2:end-1]))
            # println("dTmax=",dTmax," dt=",dt/seconds_in_year/1e3," kyr")
            dt = min(dt,dTmax < 10.0 ? dt : dt*10.0/dTmax)
            if dTmax < 10.0
                break
            end
        end

        dT_subgrid_node = subgrid_temperature_relaxation_center!(markers,grid,Tlast,Cp_c[1,1],kThermal[1,1],dt)
        dT_remaining = dTemp - dT_subgrid_node

        cell_center_change_to_markers!(markers,grid,dT_remaining,"T")
        
        if time == 0.0 || mod(itime,10) == 0 || true
            last_plot = time 
            # Gird output
            name1 = @sprintf("%s/viz.%04d.vtr",output_dir,iout)
            println("Writing visualization file = ",name1)
            vn = velocity_to_basic_nodes(grid,vxc,vyc)
            visualization(grid,rho_c,eta_s,vn,P,Tnew[2:end-1,2:end-1],time/seconds_in_year;filename=name1)
            # Markers output
            # name2 = @sprintf("%s/markers.%04d.vtp",output_dir,iout)
            # println("Writing visualization file = ",name2)
            # visualization(markers,time/seconds_in_year;filename=name2)
            iout += 1
        end

        mat, = marker_to_stag(markers,grid,markers.integers[[markers.integerFields["material"]],:],"center")
        ocean_ice_interface = get_interface(grid,mat,1.5)
        air_ice_interface = get_interface(grid,mat,2.5)
        max_ice_shell_thickness = maximum(ocean_ice_interface)-maximum(air_ice_interface)
        avg_ice_shell_thickness = mean(ocean_ice_interface)-mean(air_ice_interface)
        Af = max_ice_shell_thickness-avg_ice_shell_thickness
        i_A = @sprintf("%.6g",Ai/1e3)
        f_A = @sprintf("%.6g",Af/1e3)
        println("Initial Amplitude: $i_A (km), Finial Amplitude: $f_A (km)")
        # Checking Termination Criteria, time is in Myr
        if time >= max_time || itime >= max_step || Af/Ai <= 1/exp(1)
            terminate = true
            # println("Finished Step ",itime," time=",time/seconds_in_year/1e3," kyr")
        end     

        # println("Min/Max velocity: ",minimum(vyc)," ",maximum(vyc))            
        # Moving the markers and advancing to the next timestep
        move_markers_rk4!(markers,grid,vx,vy,dt,continuity_weight=1/3)
        time += dt
        itime += 1
        println("Finished Step ",itime," time=",time/seconds_in_year/1e3," kyr")
        mat, = marker_to_stag(markers,grid,markers.integers[[markers.integerFields["material"]],:],"center")
        ocean_ice_interface = get_interface(grid,mat,1.5)
        air_ice_interface = get_interface(grid,mat,2.5)
        append!(time_plot,time)
        # append!(topography,[ocean_ice_interface])
        append!(amplitude,Af)
    end
    return grid,i_mat,mat,time_plot,time,itime,amplitude
end

include("Topo.jl")
include("Outputs.jl")

function model_run()
    irun = 1
    top_dir = mk_modelrun_dir()
    sub_dir,sub_dir_plots,sub_dir_data = mk_output_dir(top_dir,irun)  
    # nlambda = wavelength_length
    # nhice = ice_length
    
    local lambda = 1e5
    local hice = 1e4
    
    # local lambda = range(userinput_wavelength[1],userinput_wavelength[2],nlambda)
    # local hice =  range(userinput_ice[1],userinput_ice[2],nhice)
    # t_halfspace = zeros(nlambda,nhice)
    # t_vis = zeros(nlambda,nhice)
    # t_tic = zeros(nlambda,nhice)
    # t_vis_fitted_amp_time = zeros(nlambda,nhice)
    # t_vis_fitted_amp = zeros(nlambda,nhice)
    
    options = Dict()
    options["wavelength"] = 1e5
    options["ice thickness"] = 1e4
    options["amplitude"] = 0.10*options["ice thickness"]
    options["surface depth"] = options["amplitude"]    
    
    t_halfspace = []
    t_vis = []
    t_tic = []
    # t_vis_fitted_time = []
    # t_vis_fitted_amp = []
    
    println("Starting model execution...")
    open(sub_dir*"/output.txt", "w") do out
        redirect_stdout(out) do
            println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)")   
            model_runtime(sub_dir,"Start")
            grid,i_mat,mat,times,time,itime,amplitude = run(options)
            model_runtime(sub_dir,"End")
            air_ice_interface = get_interface(grid,mat,2.5)
            ocean_ice_interface = get_interface(grid,mat,1.5)
            i_air_ice_interface = get_interface(grid,i_mat,2.5)
            i_ocean_ice_interface = get_interface(grid,i_mat,1.5)
            max_ice_shell_thickness = maximum(ocean_ice_interface)-maximum(air_ice_interface)
            avg_ice_shell_thickness = mean(ocean_ice_interface)-mean(air_ice_interface)
            Af = max_ice_shell_thickness-avg_ice_shell_thickness
            ths = get_halfspace_time_viscous(options["wavelength"])
            append!(t_halfspace,ths)
            tvis = get_numerical_time_viscous(options["amplitude"],Af,last(times))
            append!(t_vis,tvis)
            rate = get_thickening_rate(options["ice thickness"])
            ttic = get_thickening_time(options["ice thickness"],rate)
            append!(t_tic,ttic)
            # fitted_time,fitted_amp = amp_fitting_data(amplitude,times,itime,sub_dir_plots)
            # append!(t_vis_fitted_time,fitted_time)
            # fitted_amp_time = get_numerical_time_viscous(options["amplitude"],fitted_amp,last(times))
            # append!(t_vis_fitted_amp,fitted_amp_time)
            get_topography_plots(grid,i_mat,mat,i_air_ice_interface,air_ice_interface,i_ocean_ice_interface,ocean_ice_interface,times,time,itime,sub_dir_plots)
        end
    end
    println("Model ran successfully. Outputs saved to output.txt")
    # for i in 1:nlambda
    #     for j in 1:nhice
    #         options["wavelength"] = lambda[i]*1e3
    #         options["ice thickness"] = hice[j]*1e3
    #         if options["wavelength"] >= options["ice thickness"]
    #             sub_dir,sub_dir_plots,sub_dir_data = mk_output_dir(top_dir,irun)             
    #             options["amplitude"] = 0.10*options["ice thickness"]
    #             options["surface depth"] = options["amplitude"] 
    #             open(sub_dir*"/output.txt", "w") do out
    #                 redirect_stdout(out) do
    #                     println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)")   
    #                     model_runtime(sub_dir,"Start")
    #                     grid,i_mat,mat,times,time,itime,amplitude = run(options)
    #                     model_runtime(sub_dir,"End")
    #                     air_ice_interface = get_interface(grid,mat,2.5)
    #                     ocean_ice_interface = get_interface(grid,mat,1.5)
    #                     i_air_ice_interface = get_interface(grid,i_mat,2.5)
    #                     i_ocean_ice_interface = get_interface(grid,i_mat,1.5)
    #                     max_ice_shell_thickness = maximum(ocean_ice_interface)-maximum(air_ice_interface)
    #                     avg_ice_shell_thickness = mean(ocean_ice_interface)-mean(air_ice_interface)
    #                     Af = max_ice_shell_thickness-avg_ice_shell_thickness
    #                     ths = get_halfspace_time_viscous(options["wavelength"])
    #                     t_halfspace[i,j] = ths
    #                     tvis = get_numerical_time_viscous(options["amplitude"],Af,last(times))
    #                     t_vis[i,j] = tvis
    #                     ttic = get_thickening_time(options["ice thickness"])
    #                     t_tic[i,j] = ttic
    #                     fitted_time,fitted_amp = amp_fitting_data(amplitude,times,itime,sub_dir_plots)
    #                     t_vis_fitted_amp_time[i,j] = fitted_time
    #                     t_vis_fitted_amp[i,j] = get_numerical_time_viscous(options["amplitude"],fitted_amp,last(times))
    #                     get_topography_plots(grid,i_mat,mat,i_air_ice_interface,air_ice_interface,i_ocean_ice_interface,ocean_ice_interface,times,time,itime,sub_dir_plots)
    #                 end
    #             end
    #             irun += 1
    #         end
    #     end
    # end
    return lambda,hice,t_halfspace,t_vis,t_tic,t_vis_fitted_time,t_vis_fitted_amp,sub_dir_data
end

model_run()

# try
#     results = model_run()
#     data_to_hdf5_file(results[1],results[2],results[3],results[4],results[5],results[6],results[7],results[8])
# catch e
#     println("Model encountered an error. Error details saved to error_log.txt")
#     open("error_log.txt","w") do out
#         redirect_stdout(out) do
#             showerror(stdout,e,catch_backtrace())
#         end
#     end
#     @error "Something went wrong" exception=(e, catch_backtrace())
# end

# open(sub_dir*"/Output.txt","w") do io
#     redirect_stdout(io) do
#         options["amplitude"] = 0.10*options["ice thickness"]
#         options["surface depth"] = options["amplitude"] 
#         println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)")                     
#         model_runtime(sub_dir,"Start")
#         grid,i_mat,mat,times,time,itime,amplitude = run(options)
#         model_runtime(sub_dir,"End")
#         air_ice_interface = get_interface(grid,mat,2.5)
#         ocean_ice_interface = get_interface(grid,mat,1.5)
#         i_air_ice_interface = get_interface(grid,i_mat,2.5)
#         i_ocean_ice_interface = get_interface(grid,i_mat,1.5)
#         max_ice_shell_thickness = maximum(ocean_ice_interface)-maximum(air_ice_interface)
#         avg_ice_shell_thickness = mean(ocean_ice_interface)-mean(air_ice_interface)
#         Af = max_ice_shell_thickness-avg_ice_shell_thickness
#         ths = get_halfspace_time_viscous(options["wavelength"])
#         t_halfspace[i,j] = ths
#         tvis = get_numerical_time_viscous(options["amplitude"],Af,last(times))
#         t_vis[i,j] = tvis
#         ttic = get_thickening_time(options["ice thickness"])
#         t_tic[i,j] = ttic
#         fitted_time,fitted_amp = amp_fitting_data(amplitude,times,itime,sub_dir_plots)
#         t_vis_fitted_amp_time[i,j] = fitted_time
#         t_vis_fitted_amp[i,j] = get_numerical_time_viscous(options["amplitude"],fitted_amp,last(times))
#         get_topography_plots(grid,i_mat,mat,i_air_ice_interface,air_ice_interface,i_ocean_ice_interface,ocean_ice_interface,times,time,itime,sub_dir_plots)
#         # topography,time,itime,sub_dir_plots)     
#         showerror(io,catch_backtrace())
#     end
# end 