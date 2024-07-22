if length(ARGS) > 4
    error("specify proper input arguments for range function for ice shell thickness and wavelength topogaphy for ocean-ice interface")
else
    ice_shell_thickness = parse(Float64,ARGS[1])
    wavelength = parse(Float64,ARGS[2])
    percent_amplitude = parse(Float64,ARGS[3])
    top_dir = ARGS[4]
    println("Model run for ice shell thickness of $ice_shell_thickness, wavelength of $wavelength, amplitude percentage of $percent_amplitude")
end

### Model agruments ###
options = Dict()
options["latent heat of fusion"] = 3.34e5 #J/kg
options["specific heat of ice"] = 2.1e3 # J/kg*K (ice)
options["density of ice"] = 1e3 # kg/m^3
options["thermal conductivity of ice"] = 2.2 # W/m*K
options["thermal diffusivity"] = options["thermal conductivity of ice"] / (options["density of ice"]*options["specific heat of ice"]) # m^2/s
options["Tm"] = 273.0 # K
options["ny"] = 51
options["markx"] = 6
options["marky"] = 6
options["hice"] = ice_shell_thickness*1e3
options["wavelength"] = wavelength*1e3

# Importing (using/include) packages and files needed for the code to run
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using WriteVTK
using Printf
using Statistics
using Dates
using SpecialFunctions
using Roots
using NLsolve
using Printf
using HDF5
using EasyFit
include("Grid.jl")
include("Markers.jl")
include("Stokes.jl")
include("GridOperations.jl")
# Note: that we import pyplot last to avoid a name conflict with grid
using PyPlot
include("Visualization.jl")
include("Topo.jl")
include("Outputs.jl")
include("FittingData.jl")
include("TemperatureEntropy.jl")
include("ModelPlots.jl")
include("Outputs_sbatch.jl")

### Setting Up the Initial Conditions ###
include("InitialConditions.jl")

function ice_viscosity(T::Float64)
    Q = 40000.0 # Activation Enegry (J/mol)
    R_cont = 8.314 # Gas Constant (J/molK)
    meltingpoint_viscosity = 1e14
    ice_vis = meltingpoint_viscosity*exp((Q*(273.0-T))/(R_cont*(273.0*T)))
    upperlimit = 1e25
    lowerlimit = meltingpoint_viscosity
    if ice_vis < lowerlimit
        ice_vis = lowerlimit
    elseif ice_vis > upperlimit
        ice_vis = upperlimit
    end
    return ice_vis
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

### Updating Schemes ###
## starts here ##
function update_marker_prop!(markers::Markers)
    rho = markers.scalarFields["rho"]
    eta = markers.scalarFields["eta"]
    T = markers.scalarFields["T"]
    S = markers.scalarFields["S"]
    X = markers.scalarFields["X"]
    Threads.@threads for i in 1:markers.nmark
        if markers.scalars[X,i] <= 0.0
            markers.scalars[rho,i] = 920.0
        elseif markers.scalars[X,i] >= 1.0
            markers.scalars[rho,i] = 1000.0 # kg/m^3
        else
            markers.scalars[rho,i] = 920.0 + (1000.0-920.0)*markers.scalars[X,i] # kg/m^3
        end
        if markers.scalars[S,i] < 0.0
            markers.scalars[eta,i] = ice_viscosity(markers.scalars[T,i])
        else
            markers.scalars[eta,i] = 1e12
        end
    end
end

function update_marker_T_X!(markers::Markers,options::Dict)
    T = markers.scalarFields["T"]
    X = markers.scalarFields["X"]
    S = markers.scalarFields["S"]
    Threads.@threads for i in 1:markers.nmark
        markers.scalars[T,i],markers.scalars[X,i] = compute_T_X_from_S((markers.scalars[S,i]),options)
    end
end
## ends here ##

### Model Setup ###
## starts here ##
function model_setup(options::Dict,plot_dir::String,io)
    W = options["wavelength"]
    H = options["hice"] + options["amplitude"] + options["hice"]/2
    ny = options["ny"]
    nx::Int64 = ceil(ny/H*W)
    gx = 0.0
    gy = 0.113

    # -1 = insulating, 1 = constant temp
    Tbctype = [-1,-1,1,1] #left, right, top, bottom
    Tbcval = [0.0,0.0,100.0,273.0] #left, right, top, bottom
    bc = BoundaryConditions(0,0,0,0) # currently does nothing but is required argument to stokes solver.
    markx = options["markx"]
    marky = options["marky"]
    seconds_in_year = 3.15e7
    plot_interval = 1e3*seconds_in_year # 1 kyr
    end_time = 3e7*seconds_in_year
    dtmax = plot_interval
    grid = CartesianGrid(W,H,nx,ny)
    println(io,"Grid resolution(ny x nx) : $ny x $nx")
    println(io,"Cell size in the x-direction is $(grid.W/(grid.nx-1))")
    println(io,"Cell size in the y-direction is $(grid.H/(grid.ny-1))")

    materials = Materials()
    println(io,"Creating Markers...")
    @time markers = Markers(grid,["alpha","T","rho","eta","Cp","Hr","kThermal","S","X"],["material"] ; nmx=markx,nmy=marky,random=false)
    println(io,"Initial condition...")
    @time initial_conditions!(markers,materials,options)

    local ice_shell_thickness = []
    local interface_topography_array = []
    local time_plot = []

    ### Setting up agruments for interface function ###
    Ai = options["amplitude"]
    local stopping_ratio::Float64

    ### Setting up agruments for termination criteria ###
    max_step::Int64=-1
    max_time::Float64=-1.0
    max_time = max_time == -1.0 ? typemax(Float64) : max_time
    max_step = max_step == -1 ? typemax(Int64) : max_step

    time = 0.0
    iout= 0
    last_plot = 0.0
    dt = seconds_in_year

    rho_c = nothing
    rho_vx = nothing
    rho_vy = nothing
    kThermal_vx = nothing
    kThermal_vy = nothing
    eta_s = nothing
    eta_n = nothing
    vxc = nothing
    vyc = nothing
    alpha = nothing
    Hr = nothing
    dTmax = nothing
    dTemp = nothing
    Tlast = nothing
    Slast = nothing
    Xlast = nothing
    q_vx = nothing
    q_vy = nothing
    Af = nothing
    melt_fraction_contour = nothing
    Tnew = nothing
    Snew = nothing
    Xnew = nothing
    P = nothing
    vx = nothing
    vy = nothing

    # Initial
    ## Transfer properties markers -> nodes ##
    # Cell Centers
    Xlast,Tlast = marker_to_stag(markers,grid,["X","T"],"center")
    # Interpolating entropy using rhoT as weight
    rhoT = markers.scalars[markers.scalarFields["rho"],:] .* markers.scalars[markers.scalarFields["T"],:]
    Slast, = marker_to_stag(markers,grid,["S"],"center",extra_weight=rhoT)
    # Projecting Slast from the cell centers to the markers
    cell_center_to_markers!(markers,grid,Slast,"S")
    # Updating Temperature and Melt fraction on the markers
    update_marker_T_X!(markers,options)

    ### Initial Plots ###
    get_plots_new(grid,Slast,Tlast,Xlast,"initial",plot_dir)

    itime = 1
    output_dir = options["visualization file path"]
    terminate = false
    while !terminate

        ## update the markers properties ##
        update_marker_prop!(markers)
        ## Transfer properties markers -> nodes ##
        # Basic Nodes
        eta_s_new, = marker_to_stag(markers,grid,["eta",],"basic",method="logarithmic")
        # Cell Centers
        rho_c_new,Cp_c_new,alpha_new,Hr_new,kThermal_new,Slast_new,Xlast_new = marker_to_stag(markers,grid,["rho","Cp","alpha","Hr","kThermal","S","X"],"center")
        eta_n_new, = marker_to_stag(markers,grid,["eta",],"center",method="logarithmic")
        # interpolate temperature using rhocp as weight
        rhocp = markers.scalars[markers.scalarFields["rho"],:] .* markers.scalars[markers.scalarFields["Cp"],:]
        Tlast_new, = marker_to_stag(markers,grid,["T"],"center",extra_weight = rhocp)
        # Vx and Vy nodes:
        rho_vx_new, = marker_to_stag(markers,grid,["rho",],"vx")
        rho_vy_new, = marker_to_stag(markers,grid,["rho",],"vy")
        kThermal_vx_new, = marker_to_stag(markers,grid,["kThermal",],"vx")
        kThermal_vy_new, = marker_to_stag(markers,grid,["kThermal",],"vy")

        # deal with any NaN values from interpolation:
        if itime > 1
            if any(isnan.(eta_s_new))
                println(io,"found nan values")
            end
            replace_nan!(eta_s,eta_s_new)
            replace_nan!(eta_n,eta_n_new)
            replace_nan!(rho_c,rho_c_new)
            replace_nan!(rho_vx,rho_vx_new)
            replace_nan!(rho_vy,rho_vy_new)
            replace_nan!(Hr,Hr_new)
            replace_nan!(alpha,alpha_new)
            replace_nan!(kThermal_vx,kThermal_vx_new)
            replace_nan!(kThermal_vy,kThermal_vy_new)
            replace_nan!(Tlast,Tlast_new)
            replace_nan!(Slast,Slast_new)
            replace_nan!(Xlast,Xlast_new)
        end

        # Copy field data
        kThermal_vx = copy(kThermal_vx_new)
        kThermal_vy = copy(kThermal_vy_new)
        rho_vx = copy(rho_vx_new)
        rho_vy = copy(rho_vy_new)
        rho_c = copy(rho_c_new)
        Hr = copy(Hr_new)
        alpha = copy(alpha_new)
        eta_s = copy(eta_s_new)
        eta_n = copy(eta_n_new)
        Tlast = copy(Tlast_new)
        Slast = copy(Slast_new)
        Xlast = copy(Xlast_new)

        for itime in 1:(itime==1 ? 2 : 1)
            # Assembling and solving the stokes equations
            L,R = form_stokes(grid,eta_s,eta_n,rho_vx,rho_vy,bc,gx,gy;dt=dt)
            stokes_solution = L\R
            vx,vy,P = unpack(stokes_solution,grid;ghost=true)

            # Obtaining velocity at the cell centers
            vxc,vyc = velocity_to_centers(grid,vx,vy)
            adiabatic_heating = compute_adiabatic_heating(grid,rho_c,Tlast,alpha,gx,gy,vxc,vyc) * 0.0
            shear_heating = compute_shear_heating(grid,vx,vy,eta_n,eta_s) * 0.0
            H = (adiabatic_heating .+ shear_heating .+ Hr) .* 0.0

            # Computing the advection timestep
            this_dtmax = min(1.2*dt,dtmax)
            dt = compute_timestep(grid,vxc,vyc;dtmax=this_dtmax,cfl=0.1)
            diff_timestep = calculate_diffusion_timestep(grid,options)
            if dt > diff_timestep
                dt = diff_timestep
            end
       end

        last_T_norm = NaN
        T_norm = NaN
        dT = nothing
        dTmax = Inf
        dS = nothing
        dSmax = Inf
        tolerance = 1e-8
        dTnorm = []
        ititer = []
        titer = 1
        max_titer = 300
        for titer=1:max_titer
            # Computing conductive heat flux (using Tlast yields an explicit scheme)
            q_vx,q_vy = compute_q_cond(grid,Tlast,kThermal_vx,kThermal_vy)

            # Computing the new entropy (using Tlast yields an explicit scheme)
            Snew = compute_S_new(grid,Tlast,rho_c,Hr,q_vx,q_vy,Slast,dt);

            # Updating the new temperature and new melt fraction from the new entropy
            Tnew,Xnew = update_T_X_from_S(Snew,options)
            Tnew,Xnew,Snew = ghost_nodes_center_TXS(grid,Tnew,Xnew,Snew,Tbctype,Tbcval,options)

            # Computing residual for entropy
            residual_S = compute_entropy_residual(grid,Tlast,rho_c,Hr,q_vx,q_vy,Slast,Snew,dt)
            Snorm = norm(residual_S[2:ny,2:nx])

            if titer > 1
                last_T_norm = T_norm;
                T_norm = norm(Tnew)
                push!(dTnorm,abs(T_norm-last_T_norm))
                push!(ititer,titer)
            else
                T_norm = norm(Tnew)
                push!(dTnorm,abs(T_norm-last_T_norm))
                push!(ititer,titer)
            end

            # Computing the maximum temperature change
            dT = Tnew - Tlast
            dTmax = maximum(abs.(dT[2:grid.ny,2:grid.nx]))

            # Computing the maximum entropy change
            dS = Snew - Slast
            dSmax = maximum(abs.(dS[2:grid.ny,2:grid.nx]))

            # Checking for convergence:
            if Snorm < tolerance
                break
            elseif titer == max_titer
                terminate = true
                @error(io,"Did not converged")
            iout += 1
            elseif any(isnan.(dT))
                terminate = true
                @error(io,"NaN or Inf apperred")
            end
        end

        # Updating entropy on the markers by projecting dS from the cell centers to the markers
        cell_center_change_to_markers!(markers,grid,dS,"S")
        update_marker_T_X!(markers,options)

        Slast_new = copy(Snew)
        Xlast_new = copy(Xnew)
        Tlast_new = copy(Tnew)

        ### Setting up agruments for Termination Criteria ###
        melt_fraction_contour = get_interface(grid,-Xnew,-0.5)
        max_ice_shell_thickness = maximum(melt_fraction_contour)
        avg_ice_shell_thickness = mean(melt_fraction_contour)
        append!(ice_shell_thickness,avg_ice_shell_thickness)
        append!(interface_topography_array,[melt_fraction_contour])
        append!(time_plot,time)

        Af = max_ice_shell_thickness-avg_ice_shell_thickness
        i_A = @sprintf("%.6g",Ai/1e3)
        f_A = @sprintf("%.6g",Af/1e3)

        if itime == 1
            stopping_ratio = Af/Ai
        end
        term_crit = 1/exp(1)
        #if isapprox(Af/Ai,(stopping_ratio*0.75)+(term_crit*0.25);atol=1e-3)
         #   println("Model progress (~25%)")
        #elseif isapprox(Af/Ai,(stopping_ratio*0.50)+(term_crit*0.50);atol=1e-3)
         #   println("Model progress (~50%)")
        #elseif isapprox(Af/Ai,(stopping_ratio*0.25)+(term_crit*0.75);atol=1e-3)
         #   println("Model progress (~75%)")
        #end

        # Checking Termination Criteria, time is in Myr, amplitude is in meters
        if time >= max_time || itime >= max_step || Af/Ai <= term_crit
            terminate = true
            ### Final Plots ###
            get_plots_new(grid,Snew,Tnew,Xnew,"final",plot_dir)
        end

        #if time == 0.0 || mod(itime,10) == 0 || terminate
        if itime == 0.0 || Af/Ai <= term_crit ||terminate
            last_plot = time
            # Gird output
            name1 = @sprintf("%s/viz.%04d.vtr",output_dir,iout)
            vn = velocity_to_basic_nodes(grid,vxc,vyc)
            viz_fields = Dict("rho"=>rho_c[2:end-1,2:end-1],"eta"=>eta_s,"velocity"=>vn,"Pressure"=>P)
            visualization(grid,viz_fields,time/seconds_in_year;filename=name1);
            # Markers output
            name2 = @sprintf("%s/markers.%04d.vtp",output_dir,iout)
            visualization(markers,time/seconds_in_year;filename=name2);
            # Additional output
            name3 = @sprintf("%s/outfields.%04d.vtr",output_dir,iout)
            output_fields = Dict("Temperature"=>Tnew[2:end-1,2:end-1],"Entropy"=>Snew[2:end-1,2:end-1],"Melt Fraction"=>Xnew[2:end-1,2:end-1]);
            visualization(grid,output_fields,time/seconds_in_year;filename=name3);
            iout += 1
        end

        # Moving the markers and advancing to the next timestep
        move_markers_rk4!(markers,grid,vx,vy,dt,continuity_weight=1.0/3.0)
        time += dt
        if mod(itime,200) == 0
            ice_shell = (ice_shell_thickness[itime] - ice_shell_thickness[1])
            ice_shell = @sprintf("%.8g",ice_shell/1e3)
            println(io,"Ice shell as thicken by $ice_shell (km)")
            println(io,"time = ",time/seconds_in_year," yr, ",time/seconds_in_year/1e3," Kyr, ",time/seconds_in_year/1e6," Myr")
            println(io,"Finished step $itime")
        end
        itime += 1
    end
    return grid,time,itime,Af,interface_topography_array,time_plot
end

function modelrun()
    sub_plots,sub_data = mk_sub_dir(top_dir)
    amp_decimal = percent_amplitude/100
    options["amplitude"] = amp_decimal*options["hice"]
    options["surface depth"] = options["amplitude"]
    options["visualization file path"] = sub_data
    io = open(top_dir*"/output.txt","w")
    println(io,"Using Wavelength: ", options["wavelength"] / 1e3, "(km)", ", ", "Using Ice Shell Thickness: ", options["hice"] / 1e3, "(km)", ", ", "Using Amplitude Percentage: $percent_amplitude%")
    grid,time,itime,Af,interface_topograhy_array,time_plot = model_setup(options,sub_plots,io);
    thickness_over_time(grid,interface_topograhy_array,time_plot,itime,sub_plots)
    t_rel = get_numerical_time_viscous(options["amplitude"],Af,time)
    t_halfspace = get_halfspace_time_viscous(options["wavelength"])
    rate = get_thickening_rate(options["hice"])
    t_tic = get_thickening_time(options["amplitude"],rate)
    t_rel_fitted = fitting_data(interface_topograhy_array,time_plot,itime,sub_plots)
    println(io,"Analytic relaxation time: ",t_halfspace,"(yr)",t_halfspace/1e3,"(kyr) or ",t_halfspace/1e6,"(Myr)")
    println(io,"Numerical relaxation time: ",t_rel,"(yr)",t_rel/1e3,"(kyr) or ",t_rel/1e6,"(Myr)")
    close(io)
    println("Model ran successfully")
    hdf5_file(options,t_halfspace,t_rel,t_tic,t_rel_fitted,top_dir)
end

try
    modelrun();
catch e
    println("Model encountered an error. Error details saved to error_log.txt")
    open("error_log.txt","w")do out
        redirect_stdout(out) do
            showerror(stdout,e,catch_backtrace())
        end
    end
    @error "Something went wrong" exception=(e,catch_backtrace())
end
