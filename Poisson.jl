module Poisson

export form_poisson_sphan

function form_poisson_sphan(grid::EquatorialGrid,rval::Float64)
    # pre-allocate arrays to store indices and values
    nn = grid.nr * grid.npsi
    row = zeros(Int64,5*nn) # up to 5 nonzeros per row
    col = zeros(Int64,5*nn)
    value = zeros(Float64, 5*nn)
    R = zeros(Float64,nn,1)
    k=1
    jmax = grid.periodic ? grid.npsi-1 : grid.npsi
    for j in 1:jmax
        for i in 1:grid.nr
            this_row = dof(i,j,grid)
            if i==1 || i == grid.nr # interior and exterior
                row[k] = this_row; col[k] = this_row; value[k] = 1.0; k+=1
                if i==1
                    R[this_row] = 2273.0;
                else
                    R[this_row] = 273.0
                end
            elseif !grid.periodic && (j==1 || j==grid.npsi)
                row[k] = this_row; col[k] = this_row; value[k] = 1.0; k+=1
                R[this_row] = 0.0;
            else
                rp = grid.rc[i+1]
                rm = grid.rc[i]
                r = grid.r[i]
                drp = i < grid.nr ? grid.r[i+1]-grid.r[i] : grid.r[i]-grid.r[i-1]
                drm = i > 1 ? grid.r[i]-grid.r[i-1] : grid.r[i+1]-grid.r[i]
                drc = 0.5*(drp+drm)                    
                dpsip = j < jmax ? grid.psi[j+1]-grid.psi[j]   : grid.psi[j]-grid.psi[j-1]
                dpsim = j > 1    ? grid.psi[j] - grid.psi[j-1] : grid.psi[j+1]-grid.psi[j]
                dpsic = 0.5*(dpsip+dpsim)

                # diagonal value
                row[k] = this_row; col[k] = this_row; value[k] = -1.0/r^2/dpsim/dpsic - 1.0/r^2/dpsip/dpsic-rp^2/r^2/drp/drc - rm^2/r^2/drm/drc; k+=1;
                # i-1,j (in)
                row[k] = this_row; col[k] = dof(i-1,j,grid); value[k] = rm^2/r^2/drm/drc; k+=1;
                # i+1,j (out)
                row[k] = this_row; col[k] = dof(i+1,j,grid); value[k] = rp^2/r^2/drp/drc; k+=1;
                # i,j-1 (left)
                row[k] = this_row; col[k] = dof(i,j-1,grid); value[k] = 1.0/r^2/dpsim/dpsic; k+=1;
                # i,j+1 (right)
                row[k] = this_row; col[k] = dof(i,j+1,grid); value[k] = 1.0/r^2/dpsip/dpsic; k+=1;
                R[this_row] = rval;
            end
        end
    end
    if grid.periodic
        j = grid.npsi
        for i in 1:grid.nr
            this_row = dof(i,j,grid)
            row[k] = this_row; col[k] = this_row;      value[k] =  1.0; k+=1;
            row[k] = this_row; col[k] = dof(i,1,grid); value[k] = -1.0; k+=1;
        end
    end
    @views row = row[1:(k-1)]
    @views col = col[1:(k-1)]
    @views value = value[1:(k-1)]
    L = sparse(row,col,value)
    return L,R
end

end