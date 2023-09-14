# This file contains code to compute derived quantities from information stored on a structured grid.

# new code - compute stress and strain rate.
function compute_strainrate(grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64})
    # Inputs:
    # grid
    # vx and vy are velocities at the velocity nodes
    # Outputs:
    # exy at the basic nodes
    # exx and eyy at the cell centers

    exy = zeros(grid.ny,grid.nx)
    for i in 1:grid.ny
        for j in 1:grid.nx
            exy[i,j] = 0.5*( (vx[i+1,j]-vx[i,j])/(grid.yc[i+1]-grid.yc[i]) + (vy[i,j+1]-vy[i,j])/(grid.xc[j+1]-grid.xc[j]) )
        end
    end
    exx = zeros(grid.ny+1,grid.nx+1)
    eyy = zeros(grid.ny+1,grid.nx+1)
    for i in 2:grid.ny # first row is outside domain
        for j in 2:grid.nx # first column is outside domain
            exx[i,j] = (vx[i,j]-vx[i,j-1])/(grid.x[j]-grid.x[j-1])
            eyy[i,j] = (vy[i,j]-vy[i-1,j])/(grid.y[i]-grid.y[i-1])
        end
    end
    return exx,eyy,exy
end

function compute_stress(grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64},etan::Matrix{Float64},etas::Matrix{Float64})
    # compute stresses.
    # Inputs:
    # grid
    # vx and vy at the velocity nodes
    # etan is viscosity at the cell centers
    # etas is viscosity at the basic nodes
    # returns:
    # sxx and syy at the cell centers
    # sxy at the basic nodes
    exx,eyy,exy = compute_strainrate(grid,vx,vy)
    sxy = 2.0*etas .* exy
    sxx = 2.0*etan .* exx
    syy = 2.0*etan .* eyy
    return sxx,syy,sxy
end

function compute_shear_heating(grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64},etan::Matrix{Float64},etas::Matrix{Float64})
    # compute shear heating at the cell centers
    # inputs: 
    # grid - the Cartesian grid
    # vx, vy - velocity at the velocity nodes
    # etan - viscosity at the cell centers
    # etas - viscosity at the basic nodes.
    #
    # 1. compute strain-rate components
    exx,eyy,exy = compute_strainrate(grid,vx,vy)
    # 2. compute stress components
    sxy = 2.0*etas .* exy
    sxx = 2.0*etan .* exx
    syy = 2.0*etan .* eyy
    # 3. average sxy*exy over each cell
    sxyexy = zeros(grid.ny+1,grid.nx+1)
    for i in 2:grid.ny
        for j in 2:grid.nx
            sxyexy[i,j] = 0.25*(sxy[i,j]*exy[i,j]+sxy[i-1,j]*exy[i-1,j]+sxy[i,j-1]*exy[i,j-1]+sxy[i-1,j-1]*exy[i-1,j-1])
        end
    end
    # 4. compute shear heating
    shear_heating = sxx .* exx + syy .* eyy + 2.0*sxyexy
    return shear_heating    
end

function compute_adiabatic_heating(grid::CartesianGrid,rho_c::Matrix{Float64},T::Matrix{Float64},alpha_c::Matrix{Float64},gx::Float64,gy::Float64,vxc::Matrix{Float64},vyc::Matrix{Float64})
    #
    # Compute the adiabatic heating at cell centers
    # Inputs:
    # grid - the cartesian grid
    # rho_c - density at cell centers
    # T - temperature at cell centers
    # alpha_c - thermal expansivity at cell centers
    # gx - gravity acting in the +x direction
    # gy - gravity acting in the +y direction
    # vxc, vyc - velocities at the cell centers
    # Returns the adiabatic heating rate
    H_a = zeros(Float64,grid.ny+1,grid.nx+1)

    for j in 2:grid.nx+1
        for i in 2:grid.ny+1
            H_a[i,j] = T[i,j]*alpha_c[i,j]*rho_c[i,j]*(gx*vxc[i,j] + gy*vyc[i,j])
        end
    end
    return H_a
end

function replace_nan!(old_field::Matrix{Float64},new_field::Matrix{Float64})
    #
    # Replace NaN values in (new_field) with values from (old_field)
    # This is useful in case a node has no associated markers. When this happens,
    # we replace the NaN value with the value in old_field
    #
    local nanind = findall(isnan.(new_field))
    new_field[nanind] = old_field[nanind]
<<<<<<< HEAD
=======
end

function get_interface(grid::CartesianGrid,mat::Matrix{Float64},contour_level::Float64)
    for i in 1:grid.nx
        j=1;
        while j < contour_level
                
        
           j+=1 
        end
    end
>>>>>>> aac6bdb6701598d804a6533dbd969be3ac08d936
end