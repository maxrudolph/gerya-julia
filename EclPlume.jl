Threads.nthreads()=24
# Define options and parse command-line arguments:
if length( ARGS ) < 3
    error("specify excess temperature (K)/lithosphere thickness (m)/plume radius (m)")
else
    Tex = parse( Float64, ARGS[1] )
    h = parse( Float64, ARGS[2] )
    r = parse( Float64, ARGS[3] )
end

seconds_in_year = 3.15e7

options = Dict()
options["nx"] = 285+1 #201
options["ny"] = 286 #571
options["markx"] = 7
options["marky"] = 7
options["W"] = 2.850e6
options["H"] = 2.850e6
options["g"] = 10.0

options["Tcmb"] = 2590.33 + 1200.
options["lithosphere thickness"] = h
options["plume radius"] = r

options["plot interval"] = 2e6*seconds_in_year
options["melting plot interval"] = 4e5*seconds_in_year
options["output directory"] = "plume_test_" * string(Tex) * "_" * string(h)
options["max time"] = 1e8*seconds_in_year
options["max step"] = 1

println("Options: ", options )

# Import necessary packages
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using WriteVTK
using Printf
using Statistics 
using Random
using NearestNeighbors
using DataFrames
using Printf
using CSV,DataFrames

include("Grid.jl")
include("GridOperations.jl")
include("Markers.jl")
include("Stokes.jl")
include("StokesCylindrical.jl")

include("Temperature.jl")
include("TemperatureCylindrical.jl")
include("melting/yasuda.jl")

# note that we import pyplot last to avoid a name conflict with grid.
using PyPlot
include("Visualization.jl")

#
# Functions related to model setup
#

using SpecialFunctions
# function to define plate cooling solution
function plate_cooling(Tsurf,Tbtm,yL0,kappa,y,t;return_gradient=false)
   # y is a vector of y coordinates where we need the solution
    T = Tsurf .+ (Tbtm-Tsurf) .* (y/yL0)
    dTdz = 0.0
    for n in 1:20 # note 20 is an index in a summation from 1 to Inf
        T += (Tbtm-Tsurf)*2/pi .* 1/n*exp(-kappa*n^2*pi^2*t/yL0^2)*sin(n*pi*y/yL0)
        dTdz += return_gradient ? (Tbtm-Tsurf)*2/pi .* 1/n*exp(-kappa*n^2*pi^2*t/yL0^2)*(n*pi/yL0)*cos(n*pi*y/yL0) : 0.0
    end
    if return_gradient
        return T,dTdz
    else
        return T
    end
end

function halfspace_cooling(Tsurf,Tmantle,kappa,y,t;return_gradient=false)
    if t == 0.0
        if y==0.0
            if return_gradient
                return Tsurf,Inf64
            else
                return Tsurf
            end
        else
            if return_gradient
                return Tsurf,0.0
            else
                return Tmantle
            end
        end
    else
        T = Tsurf + (Tmantle-Tsurf)*erf(y/(2*sqrt(kappa*t)))
        dTdz = return_gradient ? (Tmantle-Tsurf) * 1.0/2.0/sqrt(kappa*t) * 2.0/sqrt(pi) * exp(-(y/2/sqrt(kappa*t))^2) : 0.0
        if return_gradient
            return T, dTdz
        else
            return T
        end
    end
end

function halfspace_cooling_from_thickness(Tsurf,Tmantle,kappa,y,thickness)
    t = (thickness/4)^2/kappa
    return halfspace_cooling(Tsurf,Tmantle,kappa,y,t)
end

function setup_steinberger_viscosity()
   filename = "data/visc.sh08"
   tmp = CSV.read(filename,DataFrame,comment="#");
   radius_nondimensional = tmp[:,1]
   depth_m = 6.371e6 .* (1.0 .- radius_nondimensional)
   viscosity = tmp[:,2]
   return linear_interpolation(reverse(depth_m),reverse(viscosity),extrapolation_bc=Line())
end

viscosity_depth_function = setup_steinberger_viscosity()

# function to compute viscosity
function viscosity(eta0::Float64,depth::Float64,T::Float64,E::Float64 ; visc_max=1.0e25)
   # E should be given in J/mol/K
   # Expect all temperatures in kelvin.
   Tref = adiabatic_temperature(depth)
   R = 8.314 #J/mol/K
   #depth_factor = depth > 6.6e5 ? 30.0 : 1.0
    depth_factor = viscosity_depth_function(depth)
   viscosity = depth_factor*exp(E/R/Tref*(Tref/T-1))
   if viscosity > visc_max
      viscosity = visc_max
   end
   return viscosity
end

struct Materials
    # 1 - ambient mantle
    # 2 - enhanced in eclogite
    alpha::Vector{Float64}
    rho0::Vector{Float64}
    Hr::Vector{Float64}
    Cp::Vector{Float64}
    kThermal::Vector{Float64}
    eta::Vector{Float64}
    Ea::Vector{Float64}
    frac::Vector{Float64}
    function Materials()
         new([3e-5,3e-5],[3300.,3300.],[0.0,0.0],[1000.,1000.],[3.0,3.0],[1e21,1e21],[3e5,3e5],[0.0,0.0])
    end
end

# adiabat temperature
function setup_alpha()
   filename = "data/alpha.csv"
   tmp = CSV.read(filename,DataFrame,delim=",");
   depth = tmp[:,2] .* 1.e3 
   alpha = tmp[:,1]
   return linear_interpolation(depth,alpha,extrapolation_bc=Line())
end
alpha_profile = setup_alpha()
function adiabatic_temperature_profile()
    T_surf = 1400 # potential temperature
    H = options["H"]
    Cp = 1000.
    g = options["g"]
    d = LinRange(0.1,H,4000)
    T = zeros(Float64,4000)
    T[1] = T_surf
    for i = 2:4000
        depth::Float64 = d[i-1]
        alpha = alpha_profile(d[i-1])
        dTdz = alpha*g/Cp*T[i-1]
        T[i] = T[i-1] + dTdz * (d[i]-d[i-1])
    end
   return linear_interpolation(d,T,extrapolation_bc=Line())
end
adiabatic_temperature = adiabatic_temperature_profile()

# density profile
function setup_prem_density()
   filename = "data/prem_den.csv"
   tmp = CSV.read(filename,DataFrame,delim=",");
   depth = tmp[:,2] .* 1.e3 
   prem = tmp[:,1]
   return linear_interpolation(depth,prem,extrapolation_bc=Line())
end
# eclogite profile
function setup_eclogite_density()
   filename = "data/ecl_diff.csv"
   tmp = CSV.read(filename,DataFrame,delim=",");
   depth = tmp[:,2] .* 1.e3 
   ecl_diff = tmp[:,1]
   return linear_interpolation(depth,ecl_diff,extrapolation_bc=Line())
end

prem_profile = setup_prem_density()
ecl_diff = setup_eclogite_density()

# eclogite content as a function of depth
function ecl_fraction(depth::Float64,h::Float64,max_frac::Float64,min_frac::Float64,ind::Int16)
    slope = (max_frac - min_frac)/h
    b = max_frac - slope*2.85e6
    if ind == 2
        if depth >= 2.85e6-h
            ecl_fraction = slope*depth + b
        else
            ecl_fraction = min_frac
        end
    else
        ecl_fraction = 0.0
    end
end

# function to compute density
function density(mmat::Int16,alpha::Float64,depth::Float64,T::Float64,frac::Float64)
   # Reference temperature is the 1600 K adiabat.
   Tref = adiabatic_temperature(depth)
   # bottom mantle has higher eclogite content
   density =  prem_profile(depth)*(1.0-alpha*(T-Tref)) + frac*ecl_diff(depth)
   return density
end

function update_melt!(markers::Markers,dt::Float64,mask::BitVector)
    # update the melt fraction on the markers.
    T = markers.scalarFields["T"]
    melt = markers.scalarFields["Xmelt"]
    dxdt = markers.scalarFields["dXdt"]
    mat = markers.integerFields["material"]
#    carbon = markers.scalarFields["carbon"]
#    dcarbon = markers.scalarFields["dC"]
    
    Threads.@threads for i in 1:markers.nmark
        if mask[i]            
            # note - melt fraction assumes celsius temperature:  
            new_melt = markers.integers[mat,i] == 2 ? melt_fraction(melt_model,markers.x[2,i],adiabatic_temperature(markers.x[2,i])-273.0) : 0.0
            old_melt = markers.scalars[melt,i]
            markers.scalars[dxdt,i] = new_melt > old_melt ? (new_melt - old_melt)/dt : 0.0
            markers.scalars[melt,i] = new_melt
        end
    end
end

function update_melt!(markers::Markers,dt::Float64)
    mask = BitArray{1}(undef,markers.nmark)
    mask[:] .= true
    update_melt!(markers,dt,mask)
end

function update_marker_properties!(markers::Markers,materials::Materials)
    rho = markers.scalarFields["rho"]
    frac = markers.scalarFields["frac"]
    T = markers.scalarFields["T"]
    mmat = markers.integers[markers.integerFields["material"],:]
    eta = markers.scalarFields["eta"]
    alpha = markers.scalarFields["alpha"]
    Threads.@threads for i in 1:markers.nmark        
        # re-compute density using the current temperature value
        # assume reference temperature is 273.0
        # markers.scalars[rho,i] = materials.rho0[mmat[i]] # don't update density - for comparison with gerya
        # define density based on mixing of eclogite and pyrolite
        markers.scalars[alpha,i] = alpha_profile(markers.x[2,i])
        markers.scalars[rho,i] = density(mmat[i],materials.alpha[mmat[i]],markers.x[2,i],markers.scalars[T,i],markers.scalars[frac,i])
        markers.scalars[eta,i] = viscosity(materials.eta[mmat[i]],markers.x[2,i],markers.scalars[T,i],materials.Ea[mmat[i]])
    end
end

function reassimilate_lithosphere!(markers::Markers,options::Dict)
    # This function assuimilates (i.e. overprints the temperature structure of the lithosphere)
    # It does so by re-setting the temperature on the markers.
    mantle_temperature = 1660.0
    lithosphere_thickness = options["lithosphere thickness"]
    T = markers.scalarFields["T"]

    Threads.@threads for i in 1:markers.nmark
        my::Float64 = markers.x[2,i]
        if my < 2e5
            markers.scalars[T,i] = halfspace_cooling_from_thickness(273.0,mantle_temperature,1e-6,my,lithosphere_thickness)
            # plate_cooling(273.0,mantle_temperature,1.5e5,1e-6,my,50e6*3.15e7)
        end
    end
end

function initial_conditions!(markers::Markers,materials::Materials,options::Dict)
    # Define geometric properties
    lithosphere_thickness = options["lithosphere thickness"]
    
    material = markers.integerFields["material"]
    T = markers.scalarFields["T"]
    eta = markers.scalarFields["eta"]
    alpha = markers.scalarFields["alpha"]
    cp = markers.scalarFields["Cp"]
    Hr = markers.scalarFields["Hr"]
    dxdt = markers.scalarFields["dXdt"]
    frac = markers.scalarFields["frac"]
    rho = markers.scalarFields["rho"]
    #carbon = markers.scalarFields["carbon"]
    #dC = markers.scalarFields["dC"]
    h = 400.e3 # boundary layer
    r_plume = options["plume radius"]
    Threads.@threads for i in 1:markers.nmark
        mx = markers.x[1,i]
        my = markers.x[2,i]
        mr = ((mx-0)^2 + (2.85e6-my)^2)^0.5 
        
        #set material - eclogite at cmb
        if my > 2.85e6-h || mr <= r_plume
           markers.integers[material,i] = 2
        else 
           markers.integers[material,i] = 1
        end
        if my < lithosphere_thickness
            markers.scalars[T,i] = halfspace_cooling_from_thickness(273.0,adiabatic_temperature(lithosphere_thickness),1e-6,my,lithosphere_thickness)
        elseif my > 2.85e6-h
            markers.scalars[T,i] = halfspace_cooling_from_thickness(options["Tcmb"],adiabatic_temperature(options["H"]-h),1e-6,2.85e6-my,h)
        else
            markers.scalars[T,i] = adiabatic_temperature(markers.x[2,i])
        end
        T_plume = adiabatic_temperature(markers.x[2,i])+400.0        
        if mr <= r_plume && markers.scalars[T,i]<T_plume
            # plume material
            markers.scalars[T,i] = T_plume
        end        
                        
        ind = markers.integers[material,i]
        markers.scalars[alpha,i] = alpha_profile(markers.x[2,i])
        markers.scalars[eta,i] = viscosity(materials.eta[ind],markers.x[2,i],markers.scalars[T,i],materials.Ea[ind])
        markers.scalars[rho,i] = density(ind,materials.alpha[ind],markers.x[2,i],markers.scalars[T,i],ecl_fraction(markers.x[2,i],h,0.5,0.15,ind))
        markers.scalars[frac,i] = ecl_fraction(markers.x[2,i],h,0.5,0.15,ind)
        #markers.scalars[alpha,i] = materials.alpha[ind]            
        markers.scalars[cp,i] = materials.Cp[ind]  
        markers.scalars[Hr,i] = materials.Hr[ind]  
        markers.scalars[dxdt,i] = 0.0
        #markers.scalars[dC,i] = 0.0
    end
    update_marker_properties!(markers,materials)
end

# Data and functions related to melting model
melt_model = yasuda()
# melt_model.pressure(5e5)
# melt_fraction(melt_model,3e5,1400.0)
melt_model.pressure(0.99e5)

function total_melt_rate(grid::CartesianGrid,dXdt::Matrix{Float64})
    total_melt = 0.0
     for j in 2:grid.nx
        for i in 2:grid.ny
            # 2 pi r dr dz
            total_melt += dXdt[i,j] > 0.0 ? 2.0*pi*grid.xc[j]*(grid.x[j]-grid.x[j-1])*(grid.y[i]-grid.y[i-1])*dXdt[i,j] : 0.0
        end
    end
    return total_melt
end

function update_statistics(stats_file,step,time,total_melt)
    str = @sprintf("%d\t%e\t%e\n",step,time,total_melt);
    write(stats_file, str) 
    flush(stats_file)
end

function plume_model(options::Dict;max_step::Int64=-1,max_time::Float64=-1.0)
    nx = options["nx"]
    ny = options["ny"]
    W = options["W"]
    H = options["H"]
    grid = CartesianGrid(W,H,nx,ny)
    gx = 0.0
    gy = options["g"]

    Tbctype = [-1,-1,1,1] #left, right, top, bottom
    Tbcval = [0.0,0.0,273.0,options["Tcmb"]]
    bc = BoundaryConditions(0,0,0,0) # currently does nothing but is required argument to stokes solver.
    materials = Materials()
    melt_model = yasuda()
    markx = options["markx"]
    marky = options["marky"]
    target_markers = markx*marky
    min_markers = Int(floor(target_markers*0.75))
    max_markers = Int(ceil(target_markers*2.0))

    plot_interval = options["plot interval"] # plot interval in seconds
    max_time::Float64 = max_time == -1.0 ? typemax(Float64) : max_time
    max_step::Int64 = max_step == -1 ? typemax(Int64) : max_step
    
    dtmax = plot_interval
    
    println("Creating Markers...")
    @time markers = Markers(grid,["alpha","Cp","T","rho","eta","Hr","Xmelt","dXdt","frac"],["material"] ; nmx=markx,nmy=marky,random=false)
    println("Initial condition...")
    @time initial_conditions!(markers, materials, options)

    # define arrays for k, rho, cp, H at the basic nodes. Fill them with constant values for now.
    kThermal = 3.0 .*ones(grid.ny+1,grid.nx+1);

    time = 0.0
    iout=0
    last_plot = 0.0
    rho_c = nothing

    local rho_c
    local rho_vx
    local rho_vy
    local alpha
    #       local Hr
    local Cp_c
    local eta_s
    local eta_n
    local eta_vx
    local eta_vy
    local vx,vy
    local vxc=nothing
    local vyc=nothing
    local T
    local dt=1e10
    local dTmax
    local dTemp
    local Tnew=nothing
    local Tlast
    local Xlast
    local Xnew
    local dXdt
    local frac
    local reset_temperature=false

    #pardiso_solver = MKLPardisoSolver()
    #set_nprocs!(pardiso_solver, 20)
    
    output_dir = options["output directory"]
    try
        mkdir(output_dir)
    catch
        println("output directory ",output_dir," already exists.")
    end
    
    itime=1
    stats_file = open(output_dir * "/statistics.txt","w")
    for key in keys(options)
        println(stats_file,"# ",key," => ",options[key])
    end
    flush(stats_file)
    
    depth = LinRange(0,H,ny)
    T_adiabat=zeros(Float64,ny,nx)
    rho_ecl = zeros(Float64,ny,nx)
    rho_ref = zeros(Float64,ny,nx)
    #print(type(depth[1]))
    for i in 1:ny
        d::Float64 = depth[i]
        a = adiabatic_temperature(d)
        T_adiabat[i,:] .= a
        b = ecl_diff(d)
        rho_ecl[i,:] .= b
        c = prem_profile(d)
        rho_ref[i,:] .= c
    end
    
    local terminate = false
    while !terminate        
        update_marker_properties!(markers,materials)#itime==1 ? 0.0 : dt)
        # 1. Transfer properties markers -> nodes
        visc_method = "logarithmic"
        # 1a. Basic Nodes
        eta_s_new, = marker_to_stag(markers,grid,["eta",],"basic",method=visc_method);
        # 1b. Cell Centers
        rho_c_new,Cp_c_new,alpha_new,Tlast_new = marker_to_stag(markers,grid,["rho","Cp","alpha","T"],"center")
        eta_n_new, = marker_to_stag(markers,grid,["eta",],"center",method=visc_method);
        # 1c. Vx and Vy nodes:        
        eta_vx_new, = marker_to_stag(markers,grid,["eta"],"vx",method=visc_method)
        eta_vy_new, = marker_to_stag(markers,grid,["eta"],"vy",method=visc_method)
        rho_vx_new, = marker_to_stag(markers,grid,["rho",],"vx")
        rho_vy_new, = marker_to_stag(markers,grid,["rho",],"vy")
        
        # deal with any NaN values from interpolation:
        if itime > 1
            if any(isnan.(eta_s_new))
                println("found nan values")
            end
            replace_nan!(eta_s,eta_s_new)
            replace_nan!(rho_c,rho_c_new)
            #replace_nan!(Hr,Hr_new)
            replace_nan!(Cp_c,Cp_c_new)
            replace_nan!(alpha,alpha_new)
            replace_nan!(eta_n,eta_n_new)
            replace_nan!(Tlast,Tlast_new)
            # replace_nan!(Xlast,Xlast_new)            
            replace_nan!(rho_vx,rho_vx_new)
            replace_nan!(rho_vy,rho_vy_new)
            replace_nan!(eta_vx,eta_vx_new)
            replace_nan!(eta_vy,eta_vy_new)
        end
        # Copy field data 
        rho_vx = copy(rho_vx_new)
        rho_vy = copy(rho_vy_new)
        eta_vx = copy(eta_vx_new)
        eta_vy = copy(eta_vy_new)
        rho_c = copy(rho_c_new)
        #Hr = copy(Hr_new)
        Cp_c = copy(Cp_c_new)
        alpha = copy(alpha_new)
        eta_s = copy(eta_s_new)
        eta_n = copy(eta_n_new)
        Tlast = copy(Tlast_new)
        # Xlast = copy(Xlast_new)

        Tlast = ghost_temperature_center(grid,Tlast,Tbctype,Tbcval)
        
        # 2. Assemble and solve the stokes equations
        #L,R = form_stokes(grid,eta_s,eta_n,rho_vx,rho_vy,bc,gx,gy,dt=dt)
        L,R = form_stokes_cylindrical(grid,eta_s,eta_n,eta_vx,eta_vy,rho_vx,rho_vy,bc,gx,gy)
        stokes_solution = L\R
        #stokes_solution = solve(pardiso_solver,L,R) # note - problems with accuracy using pardiso
        vx,vy,P = unpack(stokes_solution,grid;ghost=true)
    
        # Get the velocity at the cell centers:
        vxc,vyc = velocity_to_centers(grid,vx,vy)
        adiabatic_heating = compute_adiabatic_heating(grid,rho_c,Tlast,alpha,gx,gy,vxc,vyc)
        shear_heating = compute_shear_heating(grid,vx,vy,eta_n,eta_s)
        H = (adiabatic_heating .+ shear_heating.*0.0)
    
        # 3. Compute the advection timestep:
        if itime > 1
            # try to increase the timestep by 20% over the previous timestep.
            this_dtmax = min(1.2*dt,dtmax) 
        else
            this_dtmax = 1e10
        end
        dt = compute_timestep(grid,vxc,vyc ; dtmax=this_dtmax,cfl=0.25)
        if dt < 0.1*seconds_in_year
            terminate=true
        end
        dTmax = Inf
        dTemp = nothing
        Tnew = nothing
        titer=1
        for titer=1:2# limit maximum temperature change
            # assemble and solve the energy equation
            println("Trying with timestep ",dt/3.15e7/1e6," Myr")
            L,R = assemble_energy_equation_cylindrical(grid,rho_c,Cp_c,kThermal,H,Tlast,dt,Tbctype,Tbcval);
            Tnew = L\R;
            #Tnew = solve(pardiso_solver,L,R);
            Tnew = reshape(Tnew,grid.ny,grid.nx);
            Tnew = ghost_temperature_center(grid,Tnew,Tbctype,Tbcval);

            T = copy(Tnew)

            dTemp = Tnew-Tlast
            # compute the maximum temperature change
            dTmax = maximum(abs.(dTemp[2:end-1,2:end-1]))
            println("dTmax=",dTmax," dt=",dt/3.15e7/1e6)
            dt = min(dt,dTmax < 10.0 ? dt : dt*10.0/dTmax)
            if dTmax < 10.0
                break
            end
        end
        if any(isnan.(markers.scalars[markers.scalarFields["rho"],:]))
            println("nan in density")
            break
        end
        dT_subgrid_node = subgrid_temperature_relaxation_center!(markers,grid,dTemp,Cp_c[1,1],kThermal[1,1],dt)
        dT_remaining = dTemp - dT_subgrid_node

        cell_center_change_to_markers!(markers,grid,dT_remaining,"T")
        
        # compute the melt fraction and carbon release on the markers using the NEW temperature.
        if itime > 1
            update_melt!(markers,dt) # compute melt on the markers (using new temperature)
            dXdt_new, = marker_to_stag(markers,grid,["dXdt"],"center")
            replace_nan!(dXdt,dXdt_new)
            #replace_nan!(dC,dC_new)
            dXdt = copy(dXdt_new)
            #dC = copy(dC_new)
            total_melt = total_melt_rate(grid,dXdt)
            println("total melt rate: ",total_melt)
        else
            total_melt = 0.0
            #total_carbon = 0.0
            dXdt = zeros(grid.ny+1,grid.nx+1)
            #dC = zeros(grid.ny+1,grid.nx+1)
        end
        
        if itime > 20 && total_melt > 0 && !reset_temperature
            # when melting begins, re-assimilate the temperature in the lithosphere.
            # only do this once.
            reassimilate_lithosphere!(markers,options)
            reset_temperature = true
        end

        compute_boundary_heat_flow(grid,Tnew,kThermal)
        # Add/remove markers. When markers are added, give them temperature using the nodal temp.
        new_markers = add_remove_markers!(markers,grid,Tnew,min_markers,target_markers,max_markers);
        update_melt!(markers,dt,new_markers); # set the correct melt fraction on new markers.

        # Check Termination Criteria
        if time >= max_time || itime >= max_step
            terminate = true
        end
                
        # Visualization Output
        this_plot_interval = total_melt > 0.0 ? options["melting plot interval"] : options["plot interval"]
        if time == 0.0 || time - last_plot >= this_plot_interval || terminate
            last_plot = time 
            # Eulerian grid output:
            name = @sprintf("%s/viz.%04d.vtr",output_dir,iout)
            println("Writing visualization fle ",name)
            vn = velocity_to_basic_nodes(grid,vxc,vyc)
	        Tn = temperature_to_basic_nodes(grid,Tnew)
            delta_T = Tn .- T_adiabat
            frac_n, = marker_to_stag(markers,grid,["frac",],"basic",method=visc_method);
            alpha_n, = marker_to_stag(markers,grid,["alpha",],"basic",method=visc_method);
            delta_rho = frac_n .* rho_ecl .- alpha_n .* delta_T .* rho_ref
            output_fields = Dict("rho"=>rho_c[2:end-1,2:end-1],"eta"=>eta_s,"velocity"=>vn,"pressure"=>P[2:end-1,2:end-1],"T"=>Tn,"dT"=>delta_T,"frac"=>frac_n,"drho"=>delta_rho,"dXdt"=>dXdt[2:end-1,2:end-1])
            @time visualization(grid,output_fields,time/seconds_in_year;filename=name)
            # Markers output:
            name1 = @sprintf("%s/markers.%04d.vtp",output_dir,iout)
            println("Writing visualization fle ",name1)
            @time visualization(markers,time/seconds_in_year;filename=name1)
            
            iout += 1
        end
        update_statistics(stats_file,itime,time,total_melt)
        
        # Move the markers and advance to the next timestep
        println("Min/Max velocity: ",minimum(vyc)," ",maximum(vyc))
        move_markers_rk4!(markers,grid,vx,vy,dt,continuity_weight=1.0/3.0)
        time += dt
        println("Finished Step ",itime," time=",time/seconds_in_year/1e6," Myr")
        itime += 1
     end
     close(stats_file)
     return grid,markers,vx,vy,vxc,vyc,rho_c,dTemp,Tnew,Tlast,time
 end

@time grid,markers,vx,vy,vxc,vyc,rho_c,dTemp,Tnew,Tlast,time = plume_model(options,max_time=options["max time"]);
#@time grid,markers,vx,vy,vxc,vyc,rho_c,dTemp,Tnew,Tlast,time = plume_model(options,max_step=1)