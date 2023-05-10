# Define a function to form the energy equation left hand side and right hand side
function assemble_energy_equation_center(grid::CartesianGrid,rho_c::Matrix{Float64},Cp_c::Matrix{Float64},kThermal::Matrix{Float64},H::Matrix{Float64},Tlast::Matrix{Float64},dt::Float64,bctype,bcval,Tsurf::Vector{Float64})
    # Assemble the left hand side for the energy equation
    # Inputs:
    # grid - this is the CartesianGrid
    # rho_c - density at the cell centers
    # Cp_c - heat capacity at the cell centers
    # kThermal- thermal conductivity at the basic nodes
    # H - volumetric rate of internal heating at the cell centers
    # Tlast - the previous-timestep temperature.
    # dt - the time step
    # bctype - a vector containing a flag (-1=prescribed gradient, +1=prescribed value) for [left,right,top,bottom]
    # bcval - a vector containing the temperature or normal gradient values at the [left, right, top, bottom]
    # Returns:
    # L,R, the left hand side and right hand side of the energy equation
    
    bcleft  = bctype[1]   # -1 = insulating, 1 = constant temp
    bcright = bctype[2]   #
    bctop   = bctype[3]
    bcbottom  = bctype[4]
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
            if i==1 # ghost nodes along top.
                row[k] = this_row
                col[k] = this_row
                val[k] = 1.0/2.0
                k+=1
                
                row[k] = this_row
                col[k] = node_index(i+1,j,grid.ny)
                val[k] = bctop/2.0
                k+=1
                
#                 R[this_row] = bcval[3]
                R[this_row] = Tsurf[j]
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
                val[k] = (rho_c[i,j]*Cp_c[i,j])/dt + kB/dxp/dxc + kA/dxm/dxc + kD/dyp/dyc + kC/dyp/dyc;
                k+=1
                # right
                row[k] = this_row;
                col[k] = j==grid.nx ? node_index(i,j,grid.ny) : node_index(i,j+1,grid.ny);
                val[k] = j==grid.nx ? bcright*kB/dxp/dxc : -kB/dxp/dxc;
                k+=1
                # left
                row[k] = this_row;
                col[k] = node_index(i,j-1,grid.ny);
                val[k] = -kA/dxm/dxc;
                k+=1
                # down (+y)
                row[k] = this_row;
                col[k] = i==grid.ny ? node_index(i,j,grid.ny) : node_index(i+1,j,grid.ny);
                val[k] = i==grid.ny ? bcbottom*kD/dyp/dyc : -kD/dyp/dyc;
                k+=1
                # up (-y)
                row[k] = this_row;
                col[k] = node_index(i-1,j,grid.ny);
                val[k] = -kC/dyp/dyc;
                k+=1
                
                R[this_row] = Tlast[i,j]*rho_c[i,j]*Cp_c[i,j]/dt + H[i,j];
                if j==grid.nx
                    R[this_row] += 2*bcval[2]*bcright*kB/dxp/dxc
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

function ghost_temperature_center(grid::CartesianGrid,T::Matrix{Float64},bctype,bcval,Tsurf::Vector{Float64})
    # Define a new grid that is (ny+1)x(nx+1) and insert the ghost temperature values.
    # bcval shoud be a vector containing the temperatures or temperature gradients 
    # along the left, right, top, and bottom (in that order)
    bcleft  = bctype[1]   # -1 = insulating, 1 = constant temp
    bcright = bctype[2]
    bctop   = bctype[3]
    bcbottom  = bctype[4]
    #bcval = [0.0,0.0,1000.0,1000.0] # left,right,top,bottom
    Tpad = Array{Float64,2}(undef,grid.ny+1,grid.nx+1)
    Tpad[1:grid.ny,1:grid.nx] = T[1:grid.ny,1:grid.nx]

    # enforce BCs along top and bottom.
    if bctop == 1
        Tpad[1,2:grid.nx] = 2.0*Tsurf[2:grid.nx] .- Tpad[2,2:grid.nx]
#         Tpad[1,2:grid.nx] = 2.0*bcval[3] .- Tpad[2,2:grid.nx]
    elseif bctop == -1
	    Tpad[1,2:grid.nx] = Tpad[2,2:grid.nx] .- ((grid.yc[2]-grid.yc[1]) * bcval[3])
    end
    # bottom
    if bcbottom == 1
        Tpad[grid.ny+1,2:grid.nx] = 2.0*bcval[4] .- Tpad[grid.ny,2:grid.nx]
    elseif bcbottom == -1
        Tpad[grid.ny+1,2:grid.nx] = Tpad[grid.ny,2:grid.nx] .+ ((grid.yc[grid.ny+1]-grid.yc[grid.ny]) * bcval[4])
    end

    # Left boundary
    if bcleft == -1
       Tpad[:,1] = Tpad[:,2] # insulating
    elseif bcleft == 1
        println("assigning left boundary temperature ",bcval[1])
       Tpad[:,1] = 2.0*bcval[1] .- Tpad[:,2]
    end

    # Right boundary
    if bcright == -1
        Tpad[:,grid.nx+1] = Tpad[:,grid.nx] # insulating
    elseif bcright == 1
        Tpad[:,grid.nx+1] = 2.0*bcval[2] .- Tpad[:,grid.nx]
    end    

    return Tpad
end

function temperature_to_basic_nodes(grid::CartesianGrid,Tc::Matrix{Float64})
    # Interpolate temperature from the cell centers to the basic nodes. 
    # Tc should be cell centered values including ghost values outside domain.
    if size(Tc,1) != grid.ny+1 || size(Tc,2) != grid.nx + 1
        error("temperature array needs to contain ghost values")
    end
    Tn = zeros(grid.ny,grid.nx)
    Threads.@threads for i in 1:grid.ny
        for j in 1:grid.nx
            Tn[i,j] = 0.25*(Tc[i,j] + Tc[i,j+1] + Tc[i+1,j] + Tc[i+1,j+1])
        end
    end
    return Tn
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
    # dTm, = marker_to_cell_center(markers,grid,dT_subgrid_m)
    dTm, = marker_to_stag(markers,grid,dT_subgrid_m,"center")
    dTm[isnan.(dTm)] .= 0.0
    return dTm
end