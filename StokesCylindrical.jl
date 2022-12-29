function form_stokes_cylindrical(grid::CartesianGrid,eta_s::Matrix,eta_n::Matrix,eta_vx::Matrix,eta_vy::Matrix,rhoX::Matrix,rhoY::Matrix,bc::BoundaryConditions,gx::Float64,gy::Float64)
    # Form the Stokes system.    
    # For cylindrical axisymmetric coordinates.
    # note that vx and x corresponds to the radial velocity and radius
    # y and vy corresponds to the axial coordinate (up/down).
    #
    # Inputs:
    # grid - the cartesian grid
    # eta_s - viscosity at the basic nodes
    # eta_n - viscosity at the cell centers
    # eta_vx - viscosity at vx nodes
    # eta_vy - viscosity at vy nodes
    # rhoX - density at the vx nodes
    # rhoY - density at the vy nodes
    # bc - a vector describing the boundary conditions along the [left,right,top,bottom]
    # gx,gy - gravitational body force in the x and y direction
    # Outputs:
    # L,R - the left hand side (matrix) and right hand side (vector) of the stokes system
    
    k::Int64 = 1 # index into dof arrays
    nx = grid.nx
    ny = grid.ny
    nn = nx*ny
    nnz = 2*11*nn + 5*nn # total number of nonzeros in matrix (not including BCs)
    row_index = zeros(Int64,nnz) # up to 5 nonzeros per row
    col_index = zeros(Int64,nnz) 
    value = zeros(Float64, nnz)
    dx = grid.W/(grid.nx-1)
    dy = grid.H/(grid.ny-1)
    kcont = 2*eta_s[1,1]/(dx+dy)# scaling factor for continuity equation
    kcont = 1e20/(dx+dy)*2
    kbond = 1.# scaling factor for dirichlet bc equations.

    R=zeros(3*nn,1)
    
    # loop over j
    for j in 1:nx
        # loop over i
        for i in 1:ny
            # dxp is the dx in the +x direction, dxm is dx in the -x direction, dxc is the spacing between cell centers
            dxp = j<nx ? grid.x[j+1] - grid.x[j]   : grid.x[j]   - grid.x[j-1]
            dxm = j>1  ? grid.x[j]   - grid.x[j-1] : grid.x[j+1] - grid.x[j]
            dxc = 0.5*(dxp+dxm)
            # dyp and dym are spacing between vx nodes in the +y and -y directions
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
                # vx1 vx(i,j-1)
                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = 2*eta_n[i,j]/dxm/dxc - 2.0/grid.x[j]*eta_vx[i,j]/(dxp+dxm)
                k+=1
                # vx2 vx(i-1,j)
                row_index[k] = this_row
                col_index[k] = vxdof(i-1,j,ny)
                value[k] = eta_s[i-1,j]/dym/dyc
                k+=1
                # vx3 vx(i,j)
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = -2*eta_n[i,j+1]/dxp/dxc -2*eta_n[i,j]/dxm/dxc - eta_s[i,j]/dyp/dyc - eta_s[i-1,j]/dym/dyc - 2.0/grid.x[j]^2*eta_vx[i,j]
                if i == ny #vx4 vx(i+1,j)
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

                # vx5 vx(i,j+1)
                row_index[k] = this_row
                col_index[k] = vxdof(i,j+1,ny)
                value[k] = 2*eta_n[i,j+1]/dxp/dxc + 2.0/grid.x[j]*eta_vx[i,j]/(dxp+dxm)
                k+=1
                # vy1 vy(i-1,j)
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

                R[this_row] = -gx*rhoX[i,j]
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
                #vy1 vy(i,j-1)
                row_index[k] = this_row
                col_index[k] = vydof(i,j-1,ny)
                value[k] = eta_s[i,j-1]/dxm/dxc - eta_vy[i,j]/grid.xc[j]/(dxp+dxm)
                k+=1
                #vy2 vy(i-1,j)
                row_index[k] = this_row
                col_index[k] = vydof(i-1,j,ny)
                value[k] = 2*eta_n[i,j]/dym/dyc
                k+=1
                #vy3
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = -2*eta_n[i+1,j]/dyp/dyc -2*eta_n[i,j]/dym/dyc - eta_s[i,j]/dxp/dxc - eta_s[i,j-1]/dxm/dxc
                if j == nx
                   # free slip - vy5 = vx3.
                   value[k] += eta_s[i,j]/dxp/dxc + eta_vy[i,j]/grid.xc[j]/(dxp+dxm)
                end
                k+=1
                
                #vy4 vy(i+1,j)
                row_index[k] = this_row
                col_index[k] = vydof(i+1,j,ny)
                value[k] = 2*eta_n[i+1,j]/dyp/dyc
                k+=1
                #vy5 vy(i,j+1)
                if j<nx
                    row_index[k] = this_row
                    col_index[k] = vydof(i,j+1,ny)
                    value[k] = eta_s[i,j]/dxp/dxc + eta_vy[i,j]/grid.xc[j]/(dxp+dxm)
                    k+=1
                end
                #vx1 vx(i,j-1)
                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = eta_s[i,j-1]/dxc/dyc - eta_vy[i,j]/grid.xc[j]/dyc/2.0
                k+=1
                #vx2 vx(i+1,j-1)
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j-1,ny)
                value[k] = -eta_s[i+1,j-1]/dxc/dyc + eta_vy[i,j]/grid.xc[j]/dyc/2.0
                k+=1
                #vx3 vx(i,j)
                row_index[k] = this_row
                col_index[k] = vxdof(i,j,ny)
                value[k] = -eta_s[i,j]/dxc/dyc - eta_vy[i,j]/grid.xc[j]/dyc/2.0
                k+=1
                #vx4 vx(i+1,j)
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j,ny)
                value[k] = eta_s[i+1,j]/dxc/dyc + eta_vy[i,j]/grid.xc[j]/dyc/2.0
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

                R[this_row] = -gy*rhoY[i,j]
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
                xc = grid.xc[j]
                
                row_index[k] = this_row
                col_index[k] = vxdof(i,j,ny)
                value[k] =  kcont/dxm*(grid.x[j]/grid.xc[j])
                k+=1

                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = -kcont/dxm*(grid.x[j-1]/grid.xc[j])
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

function compute_velocity_divergence(grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64})
    # This function computes the divergence of the velocity field.
    # I wrote the function only to verify that the cylindrical version of the continuity equation is satisfied.
    # inputs - grid, velocities at the velocity nodes.
    nx = grid.nx
    ny = grid.ny
    divv = zeros(ny-1,nx-1)
    # compute cell-by-cell velocity divergence
    # 1/r d/dr(r vr)
    for i in 2:ny
        for j in 2:nx
            dvxdx = 1.0/grid.xc[j]/(grid.x[j]-grid.x[j-1])*(grid.x[j]*vx[i,j] - grid.x[j-1]*vx[i,j-1])
            dvydy = (vy[i,j]-vy[i-1,j])/(grid.y[i]-grid.y[i-1])
            divv[i-1,j-1] = dvxdx + dvydy
        end
    end
    return divv    
end