#### begin #### 
##### Equations for setting up Stefan conidtion #####
function get_lambda1(options::Dict)
    """
    Arguments:
        options - allow the simulation of optional keyword arguments from a dictonary
    """
    L = options["latent heat of fusion"] # J/kg
    c = options["specific heat of ice"] # J/kg*K
    dT = options["Tm"]-options["To"] # K 
    f(lambda1) = L*sqrt(pi)/(c*dT)-exp(-lambda1^2)/(lambda1*erf(lambda1))
    initial_guess = 0.1
    lambda1_solution = fzero(f,initial_guess)
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
#### end ####

#### begin ####
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
        T = exp((S/Cv))*Tm
        X = 0.0
    elseif S > (Hfus/Tm)
        T = exp(((S*Tm)-Hfus)/(Cv*Tm))*Tm
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

#### begin #### 
##### Numerical equations for conductive heat flux, entropy #####
function compute_q_cond(grid::CartesianGrid,T::Matrix{Float64},k_vx::Matrix{Float64},k_vy::Matrix{Float64}) 
    # Note - this function expects T to include ghost values on all sides of the domain.
    q_vx = zeros(grid.ny,grid.nx)
    q_vy = zeros(grid.ny,grid.nx)
    for j in 1:grid.nx
        for i in 2:grid.ny
            q_vx[i,j] = -k_vx[i,j] * (T[i,j+1]-T[i,j])/(grid.xc[j+1]-grid.xc[j])
        end
    end
    for j in 2:grid.nx
        for i in 1:grid.ny
            q_vy[i,j] = -k_vy[i,j] * (T[i+1,j]-T[i,j])/(grid.yc[i+1]-grid.yc[i])            
        end
    end
    return q_vx,q_vy
end

function compute_S_new(grid::CartesianGrid,Tlast::Matrix{Float64},rho::Matrix{Float64},H::Matrix{Float64},qx::Matrix{Float64},qy::Matrix{Float64},S_old::Matrix{Float64},dt::Float64)
    S = zeros(grid.ny+1,grid.nx+1)
    for j in 2:grid.nx
        for i in 2:grid.ny
            S[i,j] = (dt/( rho[i,j] * Tlast[i,j]) ) * 
                ( -( (qx[i,j]-qx[i,j-1] )/(grid.x[j]-grid.x[j-1]) +
                (qy[i,j]-qy[i-1,j])/(grid.y[i]-grid.y[i-1]) ) +
                H[i,j] ) +
            S_old[i,j]
        end
    end
    return S
end
#### end #### 

#### begin ###
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
    results = compute_T_X_from_S.(Snew, Ref(options))
    # Extracting T and X matrices from the results
    T = [result[1] for result in results]
    X = [result[2] for result in results]
    return T,X
end
#### end ####


#### begin ####
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
        Tpad[1,2:grid.nx] = 2.0*bcval[3] .- Tpad[2,2:grid.nx]
        # Xpad[1,2:grid.nx] = 2.0*0.0 .- Xpad[2,2:grid.nx]
        Xpad[1,2:grid.nx] = 2.0.*Xpad[2,2:grid.nx] .- Xpad[2,2:grid.nx]
        # Sb = compute_S_from_T_X.(0.0,bcval[3],Ref(options))
        Sb = compute_S_from_T_X.(Xpad[1,2:grid.nx],bcval[3],Ref(options))
        Spad[1,2:grid.nx] = (2.0*Sb) .- Spad[2,2:grid.nx]
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
        # Xpad[grid.ny+1,2:grid.nx] = (2.0*1.0) .- Xpad[grid.ny,2:grid.nx]
        Xpad[grid.ny+1,2:grid.nx] = (2.0.*Xpad[grid.ny,2:grid.nx]) .- Xpad[grid.ny,2:grid.nx]
        # Sb = compute_S_from_T_X.(1.0,bcval[4],Ref(options)) 
        Sb = compute_S_from_T_X.(Xpad[grid.ny+1,2:grid.nx],bcval[4],Ref(options)) 
        Spad[grid.ny+1,2:grid.nx] = (2.0*Sb) .- Spad[grid.ny,2:grid.nx]
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