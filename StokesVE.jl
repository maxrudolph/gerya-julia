@inline node_index(i::Int64,j::Int64,ny::Int64) = ny*(j-1)+i
@inline vxdof(i::Int64,j::Int64,ny::Int64) = 3*(node_index(i,j,ny)-1)+1
@inline vydof(i::Int64,j::Int64,ny::Int64) = 3*(node_index(i,j,ny)-1)+2
@inline pdof( i::Int64,j::Int64,ny::Int64) = 3*(node_index(i,j,ny)-1)+3

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

function form_stokes(grid::CartesianGrid,eta_s::Matrix,eta_n::Matrix,mu_s::Matrix,mu_n::Matrix,rhoX::Matrix,rhoY::Matrix,bc::BoundaryConditions,gx::Float64,gy::Float64;dt::Float64=0.0)
    # Form the Stokes system.
    # Inputs:
    # grid - the cartesian grid
    # eta_s - viscosity at the basic nodes
    # eta_n - viscosity at the cell centers
    # rhoX - density at the vx nodes
    # rhoY - density at the vy nodes
    # bc - a vector describing the boundary conditions along the [left,right,top,bottom]
    # gx,gy - gravitational body force in the x and y direction
    # dt - the timestep, used in the free surface stabilization terms. dt=0.0 (default) 
    #         disables free surface stabilization.
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
                # add free surface stabilization
                drhodx = (rhoX[i,j+1]-rhoX[i,j-1])/2/dxc
                drhody = (rhoX[i+1,j]-rhoX[i-1,j])/2/dyc
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
                value[k] = -2*eta_n[i,j+1]/dxp/dxc -2*eta_n[i,j]/dxm/dxc - eta_s[i,j]/dyp/dyc - eta_s[i-1,j]/dym/dyc - drhodx*gx*dt
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
                value[k] = eta_s[i-1,j]/dxc/dyc- drhody*gx*dt/4
                k+=1
                # vy2
                row_index[k] = this_row
                col_index[k] = vydof(i,j,ny)
                value[k] = -eta_s[i,j]/dxc/dyc- drhody*gx*dt/4
                k+=1
                # vy3
                row_index[k] = this_row
                col_index[k] = vydof(i-1,j+1,ny)
                value[k] = -eta_s[i-1,j]/dxc/dyc- drhody*gx*dt/4
                k+=1
                # vy4
                row_index[k] = this_row
                col_index[k] = vydof(i,j+1,ny)
                value[k] = eta_s[i,j]/dxc/dyc- drhody*gx*dt/4
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
            
            # visco-elastic coefficient
            Z_n1::float64 = mu_n[i,j]*dt/(mu_n[i,j]*dt + eta_n[i,j]) # P1
            Z_n2::float64 = mu_n[i+1,j]*dt/(mu_n[i+1,j]*dt + eta_n[i+1,j]) # P1
            Z_s1::float64 = mu_s[i,j-1]*dt/(mu_s[i,j-1]*dt + eta_s[i,j-1]) # S1
            Z_s2::float64 = mu_s[i,j]*dt/(mu_s[i,j]*dt + eta_s[i,j]) # S2
            
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
                # add free surface stabilization
                drhodx = (rhoY[i,j+1]-rhoY[i,j-1])/2/dxc
                drhody = (rhoY[i+1,j]-rhoY[i-1,j])/2/dyc
                #vy1
                row_index[k] = this_row
                col_index[k] = vydof(i,j-1,ny)
                value[k] = eta_s[i,j-1]*Z_s1/dxm/dxc
                k+=1
                #vy2
                row_index[k] = this_row
                col_index[k] = vydof(i-1,j,ny)
                value[k] = eta_n[i,j]*Z_s2/dym/dyc
                k+=1
                #vy3
                row_index[k] = this_row
                col_index[k] = this_row
                value[k] = -eta_n[i+1,j]*Z_n2/dyp/dyc -eta_n[i,j]*Z_n1/dym/dyc - eta_s[i,j]*Z_s2/dxp/dxc - eta_s[i,j-1]*Z_s1/dxm/dxc - drhody*gy*dt
                if j == nx
                   # free slip - vx5 = vx3.
                   value[k] += eta_s[i,j]/dxp/dxc
                end
                k+=1
                
                #vy4
                row_index[k] = this_row
                col_index[k] = vydof(i+1,j,ny)
                value[k] = eta_n[i+1,j]*Z_n2/dyp/dyc
                k+=1
                #vy5
                if j<nx
                    row_index[k] = this_row
                    col_index[k] = vydof(i,j+1,ny)
                    value[k] = eta_s[i,j]*Z_s2/dxp/dxc
                    k+=1
                end
                #vx1
                row_index[k] = this_row
                col_index[k] = vxdof(i,j-1,ny)
                value[k] = eta_s[i,j-1]*Z_s1/dxc/dyc - drhodx*gy*dt/4 - eta_n[i,j]*Z_n1/dxc/dyc
                k+=1
                #vx2
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j-1,ny)
                value[k] = -eta_s[i,j-1]*Z_s1/dxc/dyc - drhodx*gy*dt/4 + eta_n[i+1,j]*Z_n2/dxc/dyc
                k+=1
                #vx3
                row_index[k] = this_row
                col_index[k] = vxdof(i,j,ny)
                value[k] = -eta_s[i,j]*Z_s2/dxc/dyc -drhodx*gy*dt/4 + eta_n[i,j]*Z_n1/dxc/dyc
                k+=1
                #vx4
                row_index[k] = this_row
                col_index[k] = vxdof(i+1,j,ny)
                value[k] = eta_s[i,j]*Z_s2/dxc/dyc - drhodx*gy*dt/4 - eta_n[i+1,j]*Z_n2/dxc/dyc
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

                # get old stress
                # tbd
                
                R[this_row] = -gy*rhoY[i,j] - (syy[i+1,j]*(1-Z_n2)-syy[i,j]*(1-Z_n1))/dyc - (sxy[i,j]*(1-Z_s2)-sxy[i,j-1]*(1-Z_s1))/dxc
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


function unpack(solution, grid::CartesianGrid; ghost::Bool=false)
    if ghost
        nx1 = grid.nx+1
        ny1 = grid.ny+1
        P = zeros(Float64,(ny1,nx1))
        vx = zeros(Float64,(ny1,nx1))
        vy = zeros(Float64,(ny1,nx1))
        ny = grid.ny
        for j in 1:grid.nx
            for i in 1:grid.ny                
                vx[i,j] = solution[vxdof(i,j,grid.ny)]
                vy[i,j] = solution[vydof(i,j,grid.ny)]
                P[i,j]  = solution[pdof(i,j,grid.ny)]            
            end
        end
        # right boundary
        j=nx1
        for i in 1:grid.ny
              vx[i,j] = 0.0
              vy[i,j] = vy[i,j-1];# free slip
        end
        i=ny1
        for j in 1:grid.nx
               vx[i,j] = vx[i-1,j];# free-slip along bottom
               vy[i,j] = 0.0
        end
    else
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
    end
    return vx,vy,P
end

function compute_timestep(grid::CartesianGrid,vxc::Matrix,vyc::Matrix;dtmax::Float64=Inf,cfl::Float64=0.5)
    # compute the maximum timestep based on cell-centered velocities in vxc and vyc and the cfl number.
    for i in 2:grid.ny
        for j in 2:grid.nx
            dx = grid.x[j]-grid.x[j-1]
            dy = grid.y[i]-grid.y[i-1]
            dtmax = min( dtmax , cfl*dx/abs(vxc[i,j]) , cfl*dy/abs(vyc[i,j]) )                            
        end
    end
    return dtmax
end