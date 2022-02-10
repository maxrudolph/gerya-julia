node_index(i::Int,j::Int,ny::Int) = ny*(j-1)+i
vxdof(i::Int,j::Int,ny::Int) = 3*(node_index(i,j,ny)-1)+1
vydof(i::Int,j::Int,ny::Int) = 3*(node_index(i,j,ny)-1)+2
pdof( i::Int,j::Int,ny::Int) = 3*(node_index(i,j,ny)-1)+3

struct BoundaryConditions
    # The intent here is that each boundary gets a flag
    # 0 = Free-slip
    # 1 = No-slip
    # other possibilities?
    top::Int
    bottom::Int
    left::Int
    right::Int
end

function form_stokes(grid::CartesianGrid,eta_s::Matrix,eta_n::Matrix,bc::BoundaryConditions,gx::Float64,gy::Float64)
    k::Int = 1 # index into dof arrays
    nx = grid.nx
    ny = grid.ny
    nn = nx*ny
    nnz = 2*11*nn + 5*nn # total number of nonzeros in matrix (not including BCs)
    row_index = zeros(Int64,nnz) # up to 5 nonzeros per row
    col_index = zeros(Int64,nnz) 
    value = zeros(Float64, nnz)
    kcont = 2*eta_s[1,1]/(grid.dx+grid.dy)# scaling factor for continuity equation
    kcont = 1e20/(grid.dx+grid.dy)*2
    kbond = 1.# scaling factor for dirichlet bc equations.

    R=zeros(3*nn,1)
    
    # loop over j
    for j in 1:nx
        # loop over i
        for i in 1:ny
            dxp = j<nx ? grid.x[j+1] - grid.x[j]   : grid.x[j]   - grid.x[j-1]
            dxm = j>1  ? grid.x[j]   - grid.x[j-1] : grid.x[j+1] - grid.x[j]
            dxc = 0.5*(dxp+dxm)
            dyp = i<ny ? grid.yc[i+1]- grid.yc[i]   : grid.yc[i]   -grid.yc[i-1]
            dym = i>1  ? grid.yc[i]  - grid.yc[i-1] : grid.yc[i+1] -grid.yc[i]
            dyc = 0.5*(dyp+dym)

            # discretize the x-stokes - note that numbering in comments refers to Gerya Figure 7.18a
            # and equation 7.22
            this_row = vxdof(i,j,ny)
            # Boundary cases first...            
            if j==1 || j == nx # left boundary or right boundary
                # vx = 0
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = kbond
                k+=1
                R[this_row] = 0.0 *kbond
            elseif i==1
                # dvx/dy = 0 (free slip)
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = -kbond
                k+=1
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j,ny)
                value[k] = kbond
                k+=1                
                R[this_row] = 0.0*kbond
            else
                # vx1
                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = 2*eta_n[i,j]/dxm/dxc
                k+=1
                # vx2
                row_index[k] = this_row
                col_index[k] = vxdof(i-1,j,ny)
                value[k] = eta_s[i-1,j]/dym/dyc
                k+=1
                # vx3
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = -2*eta_n[i,j+1]/dxp/dxc -2*eta_n[i,j]/dxm/dxc - eta_s[i,j]/dyp/dyc - eta_s[i-1,j]/dym/dyc
                if i == ny #vx4
                    # if i == nx, dvx/dy = 0 -> vx3 == vx4 (see Gerya fig 7.18a)
                    value[k] += eta_s[i,j]/dyp/dyc
                    k+=1
                else
                    k+=1
                    # vx4

                    # enforce dvx/dy = 0 (free slip)                
                    row_index[k] = this_row
                    col_index[k] = vxdof(i+1,j,ny)
                    value[k] = eta_s[i,j]/dyp/dyc
                    k+=1
                end

                # vx5
                row_index[k] = this_row
                col_index[k] = vxdof(i,j+1,ny)
                value[k] = 2*eta_n[i,j+1]/dxp/dxc
                k+=1
                # vy1
                row_index[k] = this_row
                col_index[k] = vydof(i-1,j,ny)
                value[k] = eta_s[i-1,j]/dxc/dyc
                k+=1
                # vy2
                row_index[k] = this_row
                col_index[k] = vydof(i,j,ny)
                value[k] = -eta_s[i,j]/dxc/dyc
                k+=1
                # vy3
                row_index[k] = this_row
                col_index[k] = vydof(i-1,j+1,ny)
                value[k] = -eta_s[i-1,j]/dxc/dyc
                k+=1
                # vy4
                row_index[k] = this_row
                col_index[k] = vydof(i,j+1,ny)
                value[k] = eta_s[i,j]/dxc/dyc
                k+=1
                # P1
                row_index[k] = this_row
                col_index[k] = pdof(i,j,ny)
                value[k] = kcont/dxc
                k+=1
                # P2
                row_index[k] = this_row
                col_index[k] = pdof(i,j+1,ny)
                value[k] = -kcont/dxc
                k+=1

                R[this_row] = -gx*(rho[i-1,j]+rho[i,j])/2.
            end
            # END X-STOKES
            
            # BEGIN Y-STOKES
            dxp = j < nx ? grid.xc[j+1] - grid.xc[j]   : grid.xc[j]  -grid.xc[j-1]
            dxm = j > 1  ? grid.xc[j]   - grid.xc[j-1] : grid.xc[j+1]-grid.xc[j]
            dxc = j > 1  ? grid.x[j]    - grid.x[j-1]  : grid.x[j+1] - grid.x[j]
            dyp = i < ny ? grid.y[i+1]  - grid.y[i]    : grid.y[i]   - grid.y[i-1]
            dym = i > 1  ? grid.y[i]    - grid.y[i-1]  : grid.y[i+1] - grid.y[i]
            dyc = i < ny ? grid.yc[i+1] - grid.yc[i]   : grid.yc[i]  - grid.yc[i-1]            
            
            this_row = vydof(i,j,ny)
            if i==1 || i == ny
                # top row / bottom row
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = kbond
                k+=1
                R[this_row] = 0.0*kbond
            elseif j==1
                # left boundary - free slip
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = kbond
                k+=1
                row_index[k] = this_row
                col_index[k] = vydof(i,j+1,ny)
                value[k] = -kbond
                k+=1
                R[this_row] = 0.0*kbond
            else            
                #vy1
                row_index[k] = this_row
                col_index[k] = vydof(i,j-1,ny)
                value[k] = eta_s[i,j-1]/dxm/dxc
                k+=1
                #vy2
                row_index[k] = this_row
                col_index[k] = vydof(i-1,j,ny)
                value[k] = 2*eta_n[i,j]/dym/dyc
                k+=1
                #vy3
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = -2*eta_n[i+1,j]/dyp/dyc -2*eta_n[i,j]/dym/dyc - eta_s[i,j]/dxp/dxc - eta_s[i,j-1]/dxm/dxc
                if j == nx
                   # free slip - vx5 = vx3.
                   value[k] += eta_s[i,j]/dxp/dxc
                end
                k+=1
                
                #vy4
                row_index[k] = this_row
                col_index[k] = vydof(i+1,j,ny)
                value[k] = 2*eta_n[i+1,j]/dyp/dyc
                k+=1
                #vy5
                if j<nx
                    row_index[k] = this_row
                    col_index[k] = vydof(i,j+1,ny)
                    value[k] = eta_s[i,j]/dxp/dxc
                    k+=1
                end
                #vx1
                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = eta_s[i,j-1]/dxc/dyc
                k+=1
                #vx2
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j-1,ny)
                value[k] = -eta_s[i+1,j-1]/dxc/dyc
                k+=1
                #vx3
                row_index[k] = this_row
                col_index[k] = vxdof(i,j,ny)
                value[k] = -eta_s[i,j]/dxc/dyc
                k+=1
                #vx4
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j,ny)
                value[k] = eta_s[i+1,j]/dxc/dyc
                k+=1
                #P1
                row_index[k] = this_row
                col_index[k] = pdof(i,j,ny)
                value[k] = kcont/dyc
                k+=1
                #P2
                row_index[k] = this_row
                col_index[k] = pdof(i+1,j,ny)
                value[k] = -kcont/dyc
                k+=1

                R[this_row] = -gy*(rho[i,j-1]+rho[i,j])/2.
            end
            # END Y-STOKES
            
            # discretize the continuity equation
            # dvx/dx + dvy/dy = 0
            this_row = pdof(i,j,ny)
            if i==1 || j == 1 || (i==2 && j == 2)
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] =  kbond
                k+=1  
                R[this_row] = 0.0
            else
                dxm = grid.x[j] - grid.x[j-1]
                dym = grid.y[i] - grid.y[i-1]
                
                row_index[k] = this_row
                col_index[k] = vxdof(i,j,ny)
                value[k] =  kcont/dxm
                k+=1

                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = -kcont/dxm
                k+=1

                row_index[k] = this_row
                col_index[k] = vydof(i,j,ny)
                value[k] =  kcont/dym
                k+=1

                row_index[k] = this_row
                col_index[k] = vydof(i-1,j,ny)
                value[k] = -kcont/dym
                k+=1
                
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = 0.0
                k+=1

                R[this_row] = 0.0
            end
            # END CONTINUITY
            
        end
    end
    @views row_index = row_index[1:(k-1)]
    @views col_index = col_index[1:(k-1)]
    @views value = value[1:(k-1)]

    L = sparse(row_index,col_index,value)
    return L,R    
end

function unpack(solution, grid::CartesianGrid)
    P = zeros(Float64,(grid.ny,grid.nx))
    vx = zeros(Float64,(grid.ny,grid.nx))
    vy = zeros(Float64,(grid.ny,grid.nx))
    ny = grid.ny
     for j in 1:grid.nx
        for i in 1:grid.ny                
            vx[i,j] = solution[vxdof(i,j,grid.ny)]
            vy[i,j] = solution[vydof(i,j,grid.ny)]
            P[i,j]  = solution[pdof(i,j,grid.ny)]            
        end
    end
    return vx,vy,P
end