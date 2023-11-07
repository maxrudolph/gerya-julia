function compute_q_cond(grid::CartesianGrid,T::Matrix{Float64},kThermal::Matrix{Float64},stringtype::String) 
    if stringtype == "vx"
        q_cond = zeros(grid.ny+1,grid.nx)
        for j in 1:grid.nx
            for i in 1:grid.ny
                q_cond[i,j] = kThermal[i,j] * (T[i,j+1]-T[i,j])/(grid.xc[j+1]-grid.xc[j])
                # if j == 1
                #     q_cond[i,j] = kThermal[i,j]*(T[i,j]/grid.xc[j])
                # else
                #     q_cond[i,j] = kThermal[i,j]*((T[i,j]-T[i,j-1])/(grid.xc[j]-grid.xc[j-1]))
                # end
            end
        end
    elseif stringtype == "vy"
        q_cond = zeros(grid.ny,grid.nx+1)
        for j in 1:grid.nx+1
            for i in 1:grid.ny
                q_cond[i,j] = kThermal[i,j] * (T[i+1,j]-T[i,j])/(grid.yc[i+1]-grid.yc[i])
                # if i == 1
                #     q_cond[i,j] = kThermal[i,j]*(T[i,j]/grid.yc[i])
                # else
                #     q_cond[i,j] = kThermal_[i,j]*((T[i,j]-T[i-1,j])/(grid.yc[i]-grid.yc[i-1]))
                # end
            end
        end
    end
    return q_cond  
end

# Define a function to form the energy equation left hand side and right hand side
function assemble_energy_equation_center(grid::CartesianGrid,Tlast::Matrix{Float64},rho_c::Matrix{Float64},H::Matrix{Float64},qx::Matrix{Float64},qy::Matrix{Float64},Slast::Matrix{Float64},dt::Float64,Snew::Matrix{Float64},bctype,bcval)
    # Assemble the left hand side for the energy equation
    # Inputs:
    # grid - this is the CartesianGrid
    # rho_c - density at the cell centers
    # H - volumetric rate of internal heating at the cell centers
    # S_old - the previous-timestep entropy.
    # dt - the time step
    # bctype - a vector containing a flag (-1=prescribed gradient, +1=prescribed value) for [left,right,top,bottom]
    # bcval - a vector containing the temperature or normal gradient values at the [left, right, top, bottom]
    # Returns:
    # L,R, the left hand side and right hand side of the energy equation

    ### Setting up boundary conditions ###
    # insluating = -1
    # constant temperature = 1
    bcleft  = bctype[1]  
    bcright = bctype[2]   
    bctop   = bctype[3]
    bcbottom  = bctype[4]
    # bcval should contain temperature or dT/dx values for left,right,top,bottom

    """
               j-1  xc(j-1)   j   xc(j+1)  j+1
                |             |             |
       i-1  (i-1,j-1)-(B)--(i-1,j)--(B)-(i-1,j+1)-
                |             |             |
                |             |             |             A = q(x)
      yc(i-1)  (A)    (C)    (A)    (C)    (A)            B = q(y)
                |             |             |             C = T,S
                |             |             |
        i   -(i,j-1)--(B)---(i,j)---(B)--(i,j+1)--
                |             |             |
                |             |             |
      yc(i+1)  (A)    (C)    (A)    (C)    (A)
                |             |             |
                |             |             |
       i+1  (i+1,j-1)-(B)--(i+1,j)--(B)-(i+1,j+1)-
                |             |             |
    """
    
    N = grid.nx*grid.ny
    row = zeros(Int64,5*N);
    col = zeros(Int64,5*N);
    val = zeros(Float64, 5*N);
    R = zeros(Float64,N,1);
    k = 1;
    
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
                
                R[this_row] = bcval[3]
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
                
                """
                        j-1   j     j   j+1   j+1
                         |          |          |
                   i-1 --|----qC----|----vy----|--
                         |          |          |
                   i     qA    T   qB    T     |
                         |          |          |
                   i   --|----qD--(i,j)--vy----|--
                         |          |          |
                   i+1   vx    c    vx   c     vx
                         |          |          |
                   i+1 --|----vy----|----vy----|--
                         |          |          |
                """
                
                qA = 0.5*(qx[i-1,j-1] + qx[i,j-1])
                qB = 0.5*(qx[i,j]     + qx[i-1,j])
                qC = 0.5*(qy[i-1,j-1] + qy[i-1,j])
                qD = 0.5*(qy[i,j-1]   + qy[i,j])

                # diagonal entry
                row[k] = this_row;
                col[k] = this_row;
                val[k] = (rho_c[i,j]*Tlast[i,j]*(Snew[i,j]/dt)) + qB/dxp/dxc + qA/dxm/dxc + qD/dyp/dyc + qC/dyp/dyc;
                k+=1
                # right
                row[k] = this_row;
                col[k] = j==grid.nx ? node_index(i,j,grid.ny) : node_index(i,j+1,grid.ny);
                val[k] = j==grid.nx ? bcright*qB/dxp/dxc : -qB/dxp/dxc;
                k+=1
                # left
                row[k] = this_row;
                col[k] = node_index(i,j-1,grid.ny);
                val[k] = -qA/dxm/dxc;
                k+=1
                # down (+y)
                row[k] = this_row;
                col[k] = i==grid.ny ? node_index(i,j,grid.ny) : node_index(i+1,j,grid.ny);
                val[k] = i==grid.ny ? bcbottom*qD/dyp/dyc : -qD/dyp/dyc;
                k+=1
                # up (-y)
                row[k] = this_row;
                col[k] = node_index(i-1,j,grid.ny);
                val[k] = -qC/dyp/dyc;
                k+=1
                
                R[this_row] = H[i,j] + (rho_c[i,j]*Tlast[i,j]*(Slast/dt));
                if j==grid.nx
                    R[this_row] += 2*bcval[2]*bcright*qB/dxp/dxc
                end
                if i==grid.ny
                    R[this_row] += 2*bcval[4]*bcbottom*qD/dyp/dyc
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





function compute_S_new(grid::CartesianGrid,T::Matrix{Float64},rho::Matrix{Float64},H::Matrix{Float64},qx::Matrix{Float64},qy::Matrix{Float64},S::Matrix{Float64})
    
    S_new[i,j] = (dt/(rho[i,j]*Tnew[i,j])) * (-((qx[i,j]-qx[i,j-1])/(grid.xc[j]-grid.xc[j-1]) + (qy[i,j]-qy[i-1,j])/(grid.yc[i]-grid.yc[i-1])) + H[i,j]) + S_old[i,j]
end 