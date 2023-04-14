# Define options and parse command-line arguments:
if length( ARGS ) < 2
    error("specify excess temperature and lithosphere thickness (in m)")
else
    Tex = parse( Float64, ARGS[1] )
    h = parse( Float64, ARGS[2] )
end

seconds_in_year = 3.15e7

options = Dict()
options["nx"] = 134 #201
options["ny"] = 381 #571
options["markx"] = 10
options["marky"] = 10
options["W"] = 1e6
options["H"] = 2.850e6
options["g"] = 10.0

options["Tcmb"] = 1300.0 + Tex + 273.0
options["lithosphere thickness"] = h
options["mantle temperature"] = 1300.0 + 273.0

options["plot interval"] = 1e6*seconds_in_year
options["melting plot interval"] = 1e5*seconds_in_year
options["output directory"] = "plume_test_" * string(Tex) * "_" * string(h)
options["max time"] = 1e8*seconds_in_year

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
    t = (thickness/2)^2/kappa
    return halfspace_cooling(Tsurf,Tmantle,kappa,y,t)
end

# function setup_steinberger_viscosity()
#    filename = "data/visc.sh08"
#    tmp = CSV.read(filename,DataFrame,comment="#");
#    radius_nondimensional = tmp[:,1]
#    depth_m = 6.371e6 .* (1.0 .- radius_nondimensional)
#    viscosity = tmp[:,2]
#    return linear_interpolation(reverse(depth_m),reverse(viscosity),extrapolation_bc=Line())
# end

# viscosity_depth_function = setup_steinberger_viscosity()

# function to compute viscosity
function viscosity(eta0::Float64,depth::Float64,T::Float64,E::Float64 ; visc_max=1.0e25)
   # E should be given in J/mol/K
   # Expect all temperatures in kelvin.
   Tref = 1300.0+273.0
   R = 8.314 #J/mol/K
   depth_factor = depth > 6.6e5 ? 20.0 : 1.0
   #depth_factor = viscosity_depth_function(depth)
   # note that eta0 = 5e19 in Leitch and Davies - quite low viscosity.
   viscosity = 5.0e19*depth_factor*exp(E/R/Tref*(Tref/T-1))
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
    function Materials()
         new([3e-5,3e-5],[3300.,3300.],[0.0,0.0],[1000.,1000.],[3.0,3.0],[1e21,1e21],[3e5,3e5])
    end
end

function update_melt!(markers::Markers,dt::Float64,mask::BitVector)
    # update the melt fraction on the markers.
    T = markers.scalarFields["T"]
    melt = markers.scalarFields["Xmelt"]
    dxdt = markers.scalarFields["dXdt"]
    mat = markers.integerFields["material"]
    carbon = markers.scalarFields["carbon"]
    dcarbon = markers.scalarFields["dC"]
    
    Threads.@threads for i in 1:markers.nmark
        if mask[i]            
            # note - melt fraction assumes celsius temperature:  
            new_melt = markers.integers[mat,i] == 2 ? melt_fraction(melt_model,markers.x[2,i],adiabatic_temperature(markers.x[2,i],markers.scalars[T,i])-273.0) : 0.0
            old_melt = markers.scalars[melt,i]
            markers.scalars[dxdt,i] = new_melt > old_melt ? (new_melt - old_melt)/dt : 0.0
            markers.scalars[melt,i] = new_melt
            if new_melt > .01
                # if melt fraction exceeds 1%, liberate all carbon
                markers.scalars[dcarbon,i] = markers.scalars[carbon,i] # change in carbon is equal to amount of carbon
                markers.scalars[carbon,i] = 0.0 # reset carbon to zero
            end
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
    T = markers.scalarFields["T"]
    mmat = markers.integers[markers.integerFields["material"],:]
    eta = markers.scalarFields["eta"]
    
    Threads.@threads for i in 1:markers.nmark        
        # re-compute density using the current temperature value
        # assume reference temperature is 273.0
        # markers.scalars[rho,i] = materials.rho0[mmat[i]] # don't update density - for comparison with gerya
        markers.scalars[rho,i] = materials.rho0[mmat[i]]*(1.0-materials.alpha[mmat[i]]*(markers.scalars[T,i]-273.0)) 
        markers.scalars[eta,i] = viscosity(materials.eta[mmat[i]],markers.x[2,i],markers.scalars[T,i],materials.Ea[mmat[i]])
    end
end

function reassimilate_lithosphere!(markers::Markers,options::Dict)
    # This function assuimilates (i.e. overprints the temperature structure of the lithosphere)
    # It does so by re-setting the temperature on the markers.
    mantle_temperature = options["mantle temperature"]
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
    mantle_temperature = options["mantle temperature"]
    
    material = markers.integerFields["material"]
    T = markers.scalarFields["T"]
    eta = markers.scalarFields["eta"]
    alpha = markers.scalarFields["alpha"]
    cp = markers.scalarFields["Cp"]
    Hr = markers.scalarFields["Hr"]
    dxdt = markers.scalarFields["dXdt"]
    carbon = markers.scalarFields["carbon"]
    dC = markers.scalarFields["dC"]

    Threads.@threads for i in 1:markers.nmark
        mx = markers.x[1,i]
        my = markers.x[2,i]
        mr = ((mx-0)^2 + (my-1.2e6)^2)^0.5 
        
        #define initial cmb hot layer geometry
        h = 1.5e5 - (1.5e5-1.12e5)*(mx)/options["W"]
                
        #set material - eclogite at cmb
        if my > 2.85e6-h # eclogite-enriched material
           markers.integers[material,i] = 2
           markers.scalars[carbon,i] = 900.0*0.15 + 137.0*0.85 # assume 900 ppm for eclogite, 137 ppm for peridotite
        else # background mantle
           markers.integers[material,i] = 1
           markers.scalars[carbon,i] = 1.0*137.0
        end
        
        if my < 6.6e5
            markers.scalars[T,i] = halfspace_cooling_from_thickness(273.0,mantle_temperature,1e-6,my,options["lithosphere thickness"])
            #plate_cooling(273.0,mantle_temperature,1.5e5,1e-6,my,50e6*3.15e7)
        else
            #markers.scalars[T,i] = 1300.0+273.0
            markers.scalars[T,i] = halfspace_cooling_from_thickness(options["Tcmb"],mantle_temperature,1e-6,options["H"]-my,h)
        end
                        
        ind = markers.integers[material,i]
        markers.scalars[eta,i] = viscosity(materials.eta[ind],markers.x[2,i],markers.scalars[T,i],materials.Ea[ind])
        markers.scalars[alpha,i] = materials.alpha[ind]            
        markers.scalars[cp,i] = materials.Cp[ind]  
        markers.scalars[Hr,i] = materials.Hr[ind]  
        markers.scalars[dxdt,i] = 0.0
        markers.scalars[dC,i] = 0.0
    end
    update_marker_properties!(markers,materials)
end

# Data and functions related to melting model
melt_model = yasuda()
# melt_model.pressure(5e5)
# melt_fraction(melt_model,3e5,1400.0)
melt_model.pressure(0.99e5)
function adiabatic_temperature(depth::Float64,T::Float64)
   # assume representative upper mantle adiabat of 0.4 K/km
   return T + 0.4/1000.0 * depth
end

function total_melt_rate(grid::CartesianGrid,dXdt::Matrix{Float64},dC::Matrix{Float64})
    total_melt = 0.0
    total_carbon = 0.0
     for j in 2:grid.nx
        for i in 2:grid.ny
            # 2 pi r dr dz
            total_melt += dXdt[i,j] > 0.0 ? 2.0*pi*grid.xc[j]*(grid.x[j]-grid.x[j-1])*(grid.y[i]-grid.y[i-1])*dXdt[i,j] : 0.0
            total_carbon += dC[i,j] > 0.0 ? 2.0*pi*grid.xc[j]*(grid.x[j]-grid.x[j-1])*(grid.y[i]-grid.y[i-1])*dC[i,j] : 0.0
        end
    end
    return total_melt,total_carbon
end

function update_statistics(stats_file,step,time,total_melt,total_carbon)
    str = @sprintf("%d\t%e\t%e\t%e\n",step,time,total_melt,total_carbon);
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
    @time markers = Markers(grid,["alpha","Cp","T","rho","eta","Hr","Xmelt","dXdt","carbon","dC"],["material"] ; nmx=markx,nmy=marky,random=false)
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
    local dC
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
        H = (adiabatic_heating .+ shear_heating).*0.0
    
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
            dXdt_new,dC_new = marker_to_stag(markers,grid,["dXdt","dC"],"center")
            replace_nan!(dXdt,dXdt_new)
            replace_nan!(dC,dC_new)
            dXdt = copy(dXdt_new)
            dC = copy(dC_new)
            total_melt,total_carbon = total_melt_rate(grid,dXdt,dC)
            println("total melt rate: ",total_melt)
        else
            total_melt = 0.0
            total_carbon = 0.0
            dXdt = zeros(grid.ny+1,grid.nx+1)
            dC = zeros(grid.ny+1,grid.nx+1)
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
            output_fields = Dict("rho"=>rho_c,"eta"=>eta_s,"velocity"=>vn,"pressure"=>P[2:end-1,2:end-1],"T"=>Tn,"dXdt"=>dXdt[2:end-1,2:end-1])
            @time visualization(grid,output_fields,time/seconds_in_year;filename=name)
            # Markers output:
            name1 = @sprintf("%s/markers.%04d.vtp",output_dir,iout)
            println("Writing visualization fle ",name1)
            @time visualization(markers,time/seconds_in_year;filename=name1)
            
            iout += 1
        end
        update_statistics(stats_file,itime,time,total_melt,total_carbon)
        
        # Move the markers and advance to the next timestep
        println("Min/Max velocity: ",minimum(vyc)," ",maximum(vyc))
        move_markers_rk4!(markers,grid,vx,vy,dt,continuity_weight=1.0/3.0)
        time += dt
        itime += 1
        println("Finished Step ",itime," time=",time/seconds_in_year/1e6," Myr")
     end
     close(stats_file)
     return grid,markers,vx,vy,vxc,vyc,rho_c,dTemp,Tnew,Tlast,time
 end

@time grid,markers,vx,vy,vxc,vyc,rho_c,dTemp,Tnew,Tlast,time = plume_model(options,max_time=options["max time"]);
#@time grid,markers,vx,vy,vxc,vyc,rho_c,dTemp,Tnew,Tlast,time = plume_model(options,max_step=1)
