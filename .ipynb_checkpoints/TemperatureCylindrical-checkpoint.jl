# Define a function to form the energy equation left hand side and right hand side
# In cylindrical coordinates!!
# x -> r
# y -> z
function assemble_energy_equation_cylindrical(grid::CartesianGrid,rho_c::Matrix{Float64},Cp_c::Matrix{Float64},kThermal::Matrix{Float64},H::Matrix{Float64},Tlast::Matrix{Float64},dt::Float64,bctype,bcval)
    # Assemble the left hand side for the energy equation
    # Inputs:
    # grid - this is the CartesianGrid
    # rho_c - density at the cell centers
    # Cp_c - heat capacity at the cell centers
    # kThermal- thermal conductivity at the basic nodes
    # H - volumetric rate of internal heating at the cell centers
    # Tlast - the previous-timestep temperature.
    # dt - the time step
    # bctype - a vector containing flags for the boundary condition at the [left, right, top, bottom]
    # 1 = fixed temperature
    # 2 = fixed dT/dn
    # bcval - a vector containing the temperature or dT/dn values at the [left, right, top, bottom]
    # Returns:
    # L,R, the left hand side and right hand side of the energy equation
    
    bcleft  = bctype[1]#-1   # -1 = constant normal temperature gradient, 1 = constant temp
    bcright = bctype[2]#-1   #
    bctop   =  bctype[3]#1
    bcbottom  = bctype[4]#1
    # bcval should contain temperature or dT/dx values for left,right,top,bottom
    
    N = grid.nx*grid.ny
    row = zeros(Int64,5*N);
    col = zeros(Int64,5*N);
    val = zeros(Float64, 5*N);
    R = zeros(Float64,N,1);
    k = 1;
    #        j-1   j     j   j+1   j+1
    #   i-1   |----kC----|----vy----|
    #         |          |          |
    #   i     kA    T   kB    T     |
    #         |          |          |
    #   i     |----kD--(i,j)--vy----|
    #         |          |          |
    #   i+1   vx    c    vx   c     vx
    #         |          |          |
    #   i+1   |----vy----|----vy----|
    
    for j in 1:grid.nx
        dxc = j>1 ? grid.x[j] - grid.x[j-1] : grid.x[j+1] - grid.x[j]
        dxp = grid.xc[j+1] - grid.xc[j]
        dxm = j>1 ? grid.xc[j]-grid.xc[j-1] : grid.xc[j+1] - grid.xc[j]
        for i in 1:grid.ny
            dyc = i>1 ? grid.y[i] - grid.y[i-1] : grid.y[i+1] - grid.y[i]
            dyp = grid.yc[i+1] - grid.yc[i]
            dym = i>1 ? grid.yc[i]-grid.yc[i-1] : grid.yc[i+1] - grid.yc[i]
            
            this_row = node_index(i,j,grid.ny);
            if i==1 # ghost nodes along top
                row[k] = this_row
                col[k] = this_row
                val[k] = 1.0 #1.0/2.0
                k+=1
                
                row[k] = this_row
                col[k] = node_index(i+1,j,grid.ny)
                val[k] = bctop #bctop/2.0
                k+=1
                
                R[this_row] = bctop == 1 ? 2.0*bcval[3] : -dyp*bcval[3]
            elseif j==1 # ghost nodes along left side.
                row[k] = this_row
                col[k] = this_row
                val[k] = 1.0/2.0
                k+=1
                
                row[k] = this_row
                col[k] = node_index(i,j+1,grid.ny)
                val[k] = bcleft/2.0
                k+=1
                
                R[this_row] = bcval[1]
            else            
                # kA, kB, kC, kD
                kA = 0.5*(kThermal[i-1,j-1] + kThermal[i,j-1])
                kB = 0.5*(kThermal[i,j]     + kThermal[i-1,j])
                kC = 0.5*(kThermal[i-1,j-1] + kThermal[i-1,j])
                kD = 0.5*(kThermal[i,j-1]   + kThermal[i,j])
                #rho_c = 0.25*(rho[i-1,j-1] + rho[i,j-1] + rho[i-1,j] + rho[i,j])
                #Cp_c = 0.25*(Cp[i-1,j-1] + Cp[i,j-1] + Cp[i-1,j] + Cp[i,j])
                
                # diagonal entry
                row[k] = this_row;
                col[k] = this_row;
                val[k] = (rho_c[i,j]*Cp_c[i,j])/dt + kB*grid.x[j]/grid.xc[j]/dxp/dxc + kA*grid.x[j-1]/grid.xc[j]/dxm/dxc + kD/dyp/dyc + kC/dym/dyc;
                k+=1
                # right
                row[k] = this_row;
                col[k] = j==grid.nx ? node_index(i,j,grid.ny) : node_index(i,j+1,grid.ny);
                val[k] = j==grid.nx ? bcright*kB*grid.x[j]/grid.xc[j]/dxp/dxc : -kB/dxp/dxc*grid.x[j]/grid.xc[j];
                k+=1
                # left
                row[k] = this_row;
                col[k] = node_index(i,j-1,grid.ny);
                val[k] = -kA*grid.x[j-1]/grid.xc[j]/dxm/dxc;
                k+=1
                # down (+y)
                row[k] = this_row;
                col[k] = i==grid.ny ? node_index(i,j,grid.ny) : node_index(i+1,j,grid.ny);
                val[k] = i==grid.ny ? bcbottom*kD/dyp/dyc : -kD/dyp/dyc;
                k+=1
                # up (-y)
                row[k] = this_row;
                col[k] = node_index(i-1,j,grid.ny);
                val[k] = -kC/dym/dyc;
                k+=1
                
                R[this_row] = Tlast[i,j]*rho_c[i,j]*Cp_c[i,j]/dt + H[i,j];
                if j==grid.nx
                    R[this_row] += 2*bcval[2]*bcright*kB*grid.x[j]/grid.xc[j]/dxp/dxc
                end
                if i==grid.ny
                    R[this_row] += 2*bcval[4]*bcbottom*kD/dyp/dyc
                end
            end
        end
    end
    row = @views row[1:k-1]
    col = @views col[1:k-1]
    val = @views val[1:k-1]
    L = sparse(row,col,val)
    return L,R
end

function subgrid_temperature_relaxation_center!(markers::Markers,grid::CartesianGrid,Tlast::Matrix,Cp,kThermal,dt::Float64)
    # Perform the sub-grid scale temperature diffusion operation
    # Inputs:
    # markers - the markers
    # grid - the grid
    # Tlast - the previous timestep temperature solution at the cell centers
    # Cp (scalar) heat capacity
    # kthermal (scalar) thermal diffusivity
    # dt - the time step
    # Returns:
    # a matrix whose values are the change in temperature at the cell centers

    dsubgrid = 1.0; # subgrid temperature diffusivity
    dT_subgrid_m = Array{Float64,2}(undef,1,markers.nmark)
    # compuate the nodal temperature on the markers.    
    cell_center_to_markers!(markers,grid,Tlast,dT_subgrid_m)

    # compute the subgrid temperature changes on the markers
    rho = markers.scalarFields["rho"]
    T = markers.scalarFields["T"]
    Threads.@threads for i in 1:markers.nmark
        dx2 = (grid.x[markers.cell[1,i]+1] - grid.x[markers.cell[1,i]])^2
        dy2 = (grid.y[markers.cell[2,i]+1] - grid.y[markers.cell[2,i]])^2
        tdiff = markers.scalars[rho,i]*Cp/kThermal / (2/dx2 + 2/dy2)
        dT_subgrid_m[i] = (dT_subgrid_m[i]-markers.scalars[T,i])*( 1.0 - exp(-dsubgrid*dt/tdiff) )
    end
    # interpolate subgrid temperature changes back onto basic nodes.
    markers.scalars[T,1:markers.nmark] += dT_subgrid_m[1,:]
    # zero out nodal values for any cells without markers (nan values)
    dTm, = marker_to_stag(markers,grid,dT_subgrid_m,"center")
    dTm[isnan.(dTm)] .= 0.0
    return dTm
    
end



function compute_boundary_heat_flow(grid::CartesianGrid,Tc::Matrix{Float64},kThermal::Matrix{Float64})
    if size(Tc,1) != grid.ny+1 || size(Tc,2) != grid.nx + 1
        error("temperature array needs to contain ghost values")
    end    
    qsurf = (Tc[2,:] - Tc[1,:]) ./ (grid.xc[2]-grid.xc[1])
    Qsurf = 0.0
    for j in 2:grid.nx
        qsurf[j] *= 0.5*(kThermal[1,j]+kThermal[2,j])
        Qsurf += qsurf[j] * 2*pi*grid.xc[j]*(grid.x[j]-grid.x[j-1])
    end
    qbtm = (Tc[2,:] - Tc[1,:]) ./ (grid.xc[2]-grid.xc[1])
    Qbtm = 0.0
    for j in 2:grid.nx
        qbtm[j] *= 0.5*(kThermal[grid.ny+1,j]+kThermal[grid.ny,j])
        Qbtm += qbtm[j] * 2*pi*grid.xc[j]*(grid.x[j]-grid.x[j-1])
    end
    println("Surface heat flow, ",Qsurf," Bottom: ",Qbtm)
end