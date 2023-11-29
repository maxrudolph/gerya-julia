function temp_of_P_S(X::Float64,S::Float64,options::Dict)
    """
    Arguments:
        X - melt fraction
        S - entropy 
        options - allow the simulation of optional keyword arguments from a dictonary
    """  
    Hfus = options["latent heat of fusion"] 
    Cv = options["specific heat"] 
    Tm = options["Tm"]
    Tref = options["Tref"]
    T = exp((S/Cv) + log(Tref) - (Hfus*X)/(Cv*Tm))
    return T
end

function entropy_of_P_T(X::Float64,T::Float64,options::Dict)
    """
    Arguments:
        X - melt fraction
        T - temperature 
        options - allow the simulation of optional keyword arguments from a dictonary
    """
    Hfus = options["latent heat of fusion"] # J/kg
    Cv = options["specific heat"] # J/kg*K
    Tm = options["Tm"] # K 
    Tref = options["Tref"] # K 
    S = Cv*(log(T) - log(Tref)) + (Hfus*X)/Tm
    return S
end

function compute_melt_fraction(T::Float64,S::Float64,options::Dict)
    Hfus = options["latent heat of fusion"] # J/kg
    Cv = options["specific heat"] # J/kg*K
    Tm = options["Tm"] # K 
    Tref = options["Tref"] # K 
    X = -Cv*(log(T) - log(Tref))*(Tm/Hfus) + (S*Tm)/Hfus
    return X
end 

function get_lambda1(options::Dict)
    """
    Arguments:
        options - allow the simulation of optional keyword arguments from a dictonary
    """
    L = options["latent heat of fusion"] # J/kg
    c = options["specific heat"] # J/kg*K
    dT = options["Tm"] - options["To"] # K 
    f(lambda1) = L * sqrt(pi) / (c * dT) - exp(-lambda1^2) / (lambda1 * erf(lambda1))
    initial_guess = 0.1
    lambda1_solution = fzero(f,initial_guess)
    return lambda1_solution
end

function get_t(lambda1::Float64,options::Dict)
    """
    Arguments:
        lambda1 - 
        options - allow the simulation of optional keyword arguments from a dictonary
    """
    kappa = options["thermal diffusivity"] # m^2/s
    ym = options["ym"]^2 # m 
    lambda = lambda1^2
    t = (ym)/(4*lambda*kappa) # seconds
    return t
end

function get_y(lambda::Float64,t::Float64,options::Dict)
    kappa = options["thermal diffusivity"] # m^2/s
    y = 2*lambda1*sqrt(kappa*t)
    return y
end 

function get_theta(y::Float64,t::Float64,lambda1::Float64)
    """
    Arguments:
        y - depth position in meters
        t - time in seconds
        lambda1 - 
    """
    kappa = options["thermal diffusivity"] # m^2/s
    eta = y/(2*sqrt(kappa*t))
    theta = erf(eta)/erf(lambda1)
    return theta
end

function stefan_initial_condition(theta::Float64,options::Dict)
    """
    Arguments:
        theta - 
        options - allow the simulation of optional keyword arguments from a dictonary
    """
    Tm = options["Tm"] # K
    To = options["To"] # K
    dT = Tm-To
    T = (theta*dT) + To
    return T
end

function compute_q_cond(grid::CartesianGrid,T::Matrix{Float64},k_vx::Matrix{Float64},k_vy::Matrix{Float64}) 
    # Note - this function expects T to include ghost values on all sides of the domain.
    q_vx = zeros(grid.ny,grid.nx)
    q_vy = zeros(grid.ny,grid.nx)
    for j in 1:grid.nx
        for i in 1:grid.ny
            q_vx[i,j] = -k_vx[i,j] * (T[i,j+1]-T[i,j])/(grid.xc[j+1]-grid.xc[j])
        end
    end
    for j in 1:grid.nx
        for i in 1:grid.ny
            q_vy[i,j] = -k_vy[i,j] * (T[i+1,j]-T[i,j])/(grid.yc[i+1]-grid.yc[i])            
        end
    end
    return q_vx,q_vy
end

function compute_S_new(grid::CartesianGrid,Tlast::Matrix{Float64},rho::Matrix{Float64},H::Matrix{Float64},qx::Matrix{Float64},qy::Matrix{Float64},S_old::Matrix{Float64},dt::Float64)
    S = zeros(grid.ny,grid.nx)
    # for j in 2:grid.nx
    #     for i in 2:grid.ny
    #         S[i,j] = begin (dt/(rho[i,j]*Tlast[i,j])) * (-((qx[i,j]-qx[i,j-1])/(grid.x[j]-grid.x[j-1]) 
    #                     + (qy[i,j]-qy[i-1,j])/(grid.y[i]-grid.y[i-1])) + H[i,j]) + S_old[i,j] end
    #     end
    # end
    # return S
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

function update_temp_from_entropy(grid::CartesianGrid,Xlast::Matrix{Float64},S_old::Matrix{Float64},options::Dict)
    T = zeros(grid.ny+1,grid.nx+1)
    for j in 1:grid.nx
        for i in 1:grid.ny
            T[i,j] = temp_of_P_S(Xlast[i,j],S_old[i,j],options)
        end
    end
    return T
end

function update_melt_fraction(grid::CartesianGrid,T_new::Matrix{Float64},S_new::Matrix{Float64},options::Dict)
    X = zeros(grid.ny,grid.nx)
    for j in 2:grid.nx
        for i in 2:grid.ny
            X[i,j] = compute_melt_fraction(T_new[i,j],S_new[i,j],options)
        end
    end
    return X
end
# function update_temp_from_entropy(grid::CartesianGrid,Xlast::Matrix{Float64},S_old::Matrix{Float64},options::Dict)
#     T = zeros(grid.ny+1,grid.nx+1)
#     for j in 2:grid.nx
#         for i in 2:grid.ny
#             T[i,j] = temp_of_P_S(Xlast[i,j],S_old[i,j],options)
#         end
#     end
#     return T
# end