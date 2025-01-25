#### start ####
##### Equations for setting up Stefan conidtion #####
function get_lambda1(options::Dict)
    """
    Arguments:
        options - allow the simulation of optional keyword arguments from a dictonary
    """
    L = options["latent heat of fusion"] # J/kg
    c = options["specific heat of ice"] # J/kg*K
    dT = options["Tm"]-options["To"] # K
    f!(lambda1) = L*sqrt(pi)/(c*dT)-exp(-lambda1^2)/(lambda1*erf(lambda1))
    initial_guess = 0.1
    lambda1_solution = fzero(f!,initial_guess)
    return lambda1_solution
end

function get_t(lambda1::Float64,options::Dict)
    """
    Arguments:
        lambda1 - constant
        options - allow the simulation of optional keyword arguments from a dictonary

    Returns:
        t - time in units of (seconds)
    """
    kappa = options["thermal diffusivity"] # m^2/s
    ym = options["ym"]^2 # m
    lambda = lambda1^2
    t = ym/(4*lambda*kappa) # seconds
    return t
end

function get_y(lambda1::Float64,t::Float64,options::Dict)
    kappa = options["thermal diffusivity"] # m^2/s
    y = 2*lambda1*sqrt(kappa*t)
    return y
end

function get_theta(y::Float64,t::Float64,lambda1::Float64)
    """
    Arguments:
        y - depth position in units of (meters)
        t - time in units of (seconds)
        lambda1 - constant
    """
    kappa = options["thermal diffusivity"] # m^2/s
    eta = y/(2*sqrt(kappa*t))
    theta = erf(eta)/erf(lambda1)
    return theta
end

function stefan_initial_condition(theta::Float64,options::Dict)
    """
    Arguments:
        theta - constant
        options - allow the simulation of optional keyword arguments from a dictonary
    Returns:
        T - temperature in units of (Kelvin)
    """
    Tm = options["Tm"] # K
    To = options["To"] # K
    dT = Tm-To
    T = (theta*dT)+To
    return T
end

function compute_stefan_temp_solution(grid::CartesianGrid,Numerical_T::Matrix{Float64},X_contour::Vector{Any},time::Vector{Any},itime::Int64)
    # Setting up Stefan conidtion
    lambda1 = get_lambda1(options)
    t = get_t(lambda1,options)
    y = get_y.(lambda1,time_plot.+t,Ref(options))

    stefan_T = zeros(grid.ny,grid.nx)
    for k in 1:itime-1
        for j in 1:grid.nx
            for i in 1:grid.ny
                theta = get_theta(grid.yc[i],time[k]+t,lambda1)
                stefan_T[i,j] = stefan_initial_condition(theta,options)
            end
        end
    end

    figure()
    title("Comparison of Temperature")
    plot(grid.yc[1:end-1]/1e3,stefan_T[:,1],"b-",label="Stefan Temperature")
    plot(grid.yc[1:end-1]/1e3,Numerical_T[1:end-1,1],"r--",label="Numerical Temperature")
    gca().set_xlabel("Depth(km)")
    gca().set_ylabel("Temperature(K)")
    legend()
    show()

    figure()
    title("Comparison of Thickness Over Time")
    plot(time_plot/3.15e7/1e6,y/1e3,"b-",label="Stefan Solution")
    plot(time_plot/3.15e7/1e6,X_contour/1e3,"r--",label="Numerical Solution")
    gca().set_ylabel(L"Ice\,Thickness\,(km)")
    gca().set_xlabel(L"Time\,(Myr)")
    legend()
    # legend(loc="upper center",bbox_to_anchor=(0.5,-0.15),ncol=5)
    show()
end
#### end ####

#### start ####
##### Equations for T,X = fcn(S) and S = fcn(T,X) #####
function compute_T_X_from_S(S::Float64,options::Dict)
    """
    Arguments:
        S - entropy in untis of (J/kg*K)
        options - allow the simulation of optional keyword arguments from a dictonary
    Returns:
        T - temperature in units of (Kelvin)
    """
    Hfus = options["latent heat of fusion"] # J/kg
    Tm = options["Tm"] # K
    Cv = options["specific heat of ice"] # J/kg*K
    if S < 0
        T = exp(S/Cv) * Tm
        X = 0.0
    elseif S > (Hfus/Tm)
        T = exp(((S*Tm)-Hfus)/(Cv*Tm)) * Tm
        X = 1.0
    else
        X = (S*Tm)/Hfus
        T = Tm
    end
    return T,X
end

function compute_S_from_T_X(X::Float64,T::Float64,options::Dict)
    """
    Arguments:
        X - melt fraction (unitless)
        T - temperature in units of (Kelvin)
        options - allow the simulation of optional keyword arguments from a dictonary
    Returns:
        S - entropy in untis of (J/kg*K)
    """
    Hfus = options["latent heat of fusion"] # J/kg
    Tm = options["Tm"] # K
    Cv = options["specific heat of ice"] # J/kg*K
    if T < Tm
        S = Cv*(log(T)-log(Tm))
    elseif T > Tm
        S = Cv*(log(T)-log(Tm))+(Hfus/Tm)
    else
        S = (Hfus/Tm)*X
    end
    return S
end
#### end ####

#### start ####
##### Numerical equations for conductive heat flux, entropy #####
function compute_q_cond(grid::CartesianGrid,T::Matrix{Float64},k_vx::Matrix{Float64},k_vy::Matrix{Float64})
    # Note - this function expects T to include ghost values on all sides of the domain.
    q_vx = zeros(grid.ny+1,grid.nx)
    q_vy = zeros(grid.ny,grid.nx+1)
    for j in 1:grid.nx
        for i in 2:grid.ny
            q_vx[i,j] = -k_vx[i,j] * ((T[i,j+1]-T[i,j])/(grid.xc[j+1]-grid.xc[j]))
        end
    end
    for j in 2:grid.nx
        for i in 1:grid.ny
            q_vy[i,j] = -k_vy[i,j] * ((T[i+1,j]-T[i,j])/(grid.yc[i+1]-grid.yc[i]))
        end
    end
    return q_vx,q_vy
end

function compute_S_new(grid::CartesianGrid,T::Matrix{Float64},rho::Matrix{Float64},H::Matrix{Float64},qx::Matrix{Float64},qy::Matrix{Float64},S_old::Matrix{Float64},dt::Float64)
    S = zeros(grid.ny+1,grid.nx+1)
    # for j in 2:grid.nx
    #     for i in 2:grid.ny
    #         S[i,j] = begin (dt/(rho[i,j]*Tlast[i,j])) * (-((qx[i,j]-qx[i,j-1])/(grid.x[j]-grid.x[j-1])
    #                     + (qy[i,j]-qy[i-1,j])/(grid.y[i]-grid.y[i-1])) + H[i,j]) + S_old[i,j] end
    #     end
    # end
    # return S
    for j in 2:grid.nx
        for i in 2:grid.ny
            S[i,j] = (dt/( rho[i,j] * T[i,j]) ) *
                ( -( (qx[i,j]-qx[i,j-1] )/(grid.x[j]-grid.x[j-1]) +
                (qy[i,j]-qy[i-1,j])/(grid.y[i]-grid.y[i-1]) ) +
                H[i,j] ) +
            S_old[i,j]
        end
    end
    return S
end

function compute_entropy_residual(grid::CartesianGrid,T::Matrix{Float64},rho::Matrix{Float64},H::Matrix{Float64},qx::Matrix{Float64},qy::Matrix{Float64},S_old::Matrix{Float64},S_new::Matrix{Float64},dt::Float64)
    return S_new .- compute_S_new(grid,T,rho,H,qx,qy,S_old,dt)
end
#### end ####

#### start ###
#### Update functions ####
function update_T_X_from_S(Snew::Matrix{Float64},options::Dict)
    """
    Arguments:
        S - a matrix of entropy in untis of (J/kg*K)
        options - allow the simulation of optional keyword arguments from a dictonary

    Returns:
        X - a martix of melt fraction (unitless)
        T - a matrix of temperature in units of (Kelvin)
    """
    # Broadcasting the function to each element in S
    results = compute_T_X_from_S.(Snew,Ref(options))
    # Extracting T and X matrices from the results
    T = [result[1] for result in results]
    X = [result[2] for result in results]
    return T,X
end
#### end ####

function calculate_diffusion_timestep(grid::CartesianGrid,options::Dict)
    # dx = grid.W / (grid.nx-1)  # Cell size in the x-direction
    # dy = grid.H / (grid.ny-1)  # Cell size in the y-direction
    dx = grid.x[2]-grid.x[1]
    dy = grid.y[2]-grid.y[1]
    diffusion_timestep_x = dx^2 / options["thermal diffusivity"]
    diffusion_timestep_y = dy^2 / options["thermal diffusivity"]
    return min(diffusion_timestep_x, diffusion_timestep_y) / 6
end

#### start ####
##### Function to compute the ghost nodes #####
function ghost_nodes_center_TXS(grid::CartesianGrid,T::Matrix{Float64},X::Matrix{Float64},S::Matrix{Float64},bctype,bcval,options::Dict)
    # along the left, right, top, and bottom (in that order)
    # -1 = insulating, 1 = constant temp
    # Extracting the boundary condition types for the left, right, top, and bottom
    bcleft  = bctype[1]
    bcright = bctype[2]
    bctop   = bctype[3]
    bcbottom  = bctype[4]

    # Creating a matrix then copying the values from the original matrix into the interior of new matrix
    Tpad = Array{Float64,2}(undef,grid.ny+1,grid.nx+1)
    Tpad[1:grid.ny,1:grid.nx] = T[1:grid.ny,1:grid.nx]
    Xpad = Array{Float64,2}(undef,grid.ny+1,grid.nx+1)
    Xpad[1:grid.ny,1:grid.nx] = X[1:grid.ny,1:grid.nx]
    Spad = Array{Float64,2}(undef,grid.ny+1,grid.nx+1)
    Spad[1:grid.ny,1:grid.nx] = S[1:grid.ny,1:grid.nx]

    # println("Before Tpad Matrix")
    # display(Tpad)
    # println("Before Spad Matrix")
    # display(Spad)
    # println("Before Xpad Matrix")
    # display(Xpad)

    # Applying the boundary condition along top of the domain
    # -1 = insulating, 1 = constant temp
    if bctop == 1
        Tpad[1,2:grid.nx] = (2.0*bcval[3]) .- Tpad[2,2:grid.nx]
        Xpad[1,2:grid.nx] = (2.0*Xpad[2,2:grid.nx][1]) .- Xpad[2,2:grid.nx]
        Sbt = compute_S_from_T_X(Xpad[1,2:grid.nx][1],bcval[3],options)
        Spad[1,2:grid.nx] = (2.0*Sbt) .- Spad[2,2:grid.nx]
    elseif bctop == -1
        # Tpad[1,2:grid.nx] = Tpad[2,2:grid.nx] .- ((grid.yc[2]-grid.yc[1]) * bcval[3])
        # Xpad[1,2:grid.nx] .= 0.0
        # Sb = compute_S_from_T_X.(Xpad[1,2:grid.nx],Tpad[1,2:grid.nx],Ref(options))
        # Spad[1,2:grid.nx] = 2.0*Sb .- Spad[2,2:grid.nx]
    end

    # Applying the boundary condition along bottom of the domain
    # -1 = insulating, 1 = constant temp
    if bcbottom == 1
        Tpad[grid.ny+1,2:grid.nx] = 2.0*bcval[4] .- Tpad[grid.ny,2:grid.nx]
        Xpad[grid.ny+1,2:grid.nx] = (2.0*Xpad[grid.ny,2:grid.nx][1]) .- Xpad[grid.ny,2:grid.nx]
        Sbb = compute_S_from_T_X(Xpad[grid.ny+1,2:grid.nx][1],bcval[4],options)
        Spad[grid.ny+1,2:grid.nx] = (2.0*Sbb) .- Spad[grid.ny,2:grid.nx]
    elseif bcbottom == -1
        # Tpad[grid.ny+1,2:grid.nx] = Tpad[grid.ny,2:grid.nx] .+ ((grid.yc[grid.ny+1]-grid.yc[grid.ny]) * bcval[4])
        # Xpad[grid.ny+1,2:grid.nx] = 2.0*1.0 .- Xpad[grid
        # Spad[grid.ny+1,2:grid.nx] = compute_S_from_T_X.(Xpad[grid.ny+1,2:grid.nx],bcval[4],Ref(options))
    end

    # Applying the boundary condition along left of the domain
    # -1 = insulating, 1 = constant temp
    if bcleft == -1
        Tpad[:,1] = Tpad[:,2] # insulating
        Xpad[:,1] = Xpad[:,2]
        Spad[:,1] = Spad[:,2]
    elseif bcleft == 1
        println("assigning left boundary temperature ",bcval[1])
        Tpad[:,1] = 2.0*bcval[1] .- Tpad[:,2]
        Xpad[:,1] = Xpad[:,2]
        Spad[:,1] = Spad[:,2]
    end

    # Applying the boundary condition along right of the domain
    # -1 = insulating, 1 = constant temp
    if bcright == -1
        Tpad[:,grid.nx+1] = Tpad[:,grid.nx] # insulating
        Xpad[:,grid.nx+1] = Xpad[:,grid.nx]
        Spad[:,grid.nx+1] = Spad[:,grid.nx]
    elseif bcright == 1
        Tpad[:,grid.nx+1] = 2.0*bcval[2] .- Tpad[:,grid.nx]
        Xpad[:,grid.nx+1] = Xpad[:,grid.nx]
        Spad[:,grid.nx+1] = Spad[:,grid.nx]
    end

    # println("After Tpad Matrix")
    # display(Tpad)
    # println("After Spad Matrix")
    # display(Spad)
    # println("After Xpad Matrix")
    # display(Xpad)

    return Tpad,Xpad,Spad
end
#### end ####

#### start ####
##### Function to compute a subgird diffusion operation with entropy #####
function subgridSdiff(grid::CartesianGrid,markers::Markers,Slast::Matrix{Float64},dt::Float64,options::Dict;diffusivity::Float64=1.0)
    """
    Arguments:
        Slast -
        dt -
        options - allow the simulation of optional keyword arguments from a dictonary

    Returns:
        dSm -
    """
    Hfus = options["latent heat of fusion"] # J/kg
    Tm = options["Tm"] # K
    Cv = options["specific heat of ice"] # J/kg*K

    # Defining d a dimensionless numerical diffusion coefficient
    d = diffusivity

    # Creating a matrix for the subgrid entropy changes on the markers
    dS_subgrid_Sm = Array{Float64,2}(undef,1,markers.nmark)
    Sm_nodal = Array{Float64,2}(undef,1,markers.nmark)
    # cell centers -> markers for Slast -> S(nodal)
    cell_center_to_markers!(markers,grid,Slast,Sm_nodal)

    # Obtainig rho, Cp, k, and S on the markers
    rho = markers.scalarFields["rho"]
    T = markers.scalarFields["T"]
    S = markers.scalarFields["S"]
    kThermal = markers.scalarFields["kThermal"]
    X = markers.scalarFields["X"]

    Threads.@threads for i in 1:markers.nmark
        dx2 = (grid.x[markers.cell[1,i]+1] - grid.x[markers.cell[1,i]])^2
        dy2 = (grid.y[markers.cell[2,i]+1] - grid.y[markers.cell[2,i]])^2
        dS::Float64=0.0
        dT::Float64=0.0
        dX::Float64=0.0
               
        for iter in 1:3
             if iter == 1
                 # initial iteration - guess heat capacity equal to Cv
                 Cpm = Cv
             else
                 # after first iteration
                 if dT == 0 || dX == 0# if there is no change in melt fraction, use Cv
                     Cpm = Cv
                 else
                     Cpm = markers.scalars[T,i]*(dS/dT)
                 end
             end

            tdiff = markers.scalars[rho,i]*Cpm/markers.scalars[kThermal,i]/(2/dx2 + 2/dy2)

                        
            dS_subgrid_Sm[i] = (Sm_nodal[i] - markers.scalars[S,i])*( 1.0 - exp(-d*dt/tdiff) )
            Si = markers.scalars[S,i] + dS_subgrid_Sm[i] # new guess of marker entropy
            Ti,Xi = compute_T_X_from_S(Si,options) # temperature consistent with new marker entropy
            dS = Si - markers.scalars[S,i] # change in marker entropy
            dT = Ti - markers.scalars[T,i] # change in marker temperature
            dX = Xi - markers.scalars[X,i]           
        end
    end
    # update the marker entropy
    markers.scalars[S,1:markers.nmark] += dS_subgrid_Sm[1,:]
    rhoT = markers.scalars[markers.scalarFields["rho"],:] .* markers.scalars[markers.scalarFields["T"],:]
    dSm, = marker_to_stag(markers,grid,dS_subgrid_Sm,"center",extra_weight=rhoT)    
    #rhoT = markers.scalars[markers.scalarFields["rho"],:] .* markers.scalars[markers.scalarFields["T"],:]
    dSm, = marker_to_stag(markers,grid,dS_subgrid_Sm,"center")#,extra_weight=rhoT)
    dSm[isnan.(dSm)] .= 0.0
    return dSm
end
#### end ####