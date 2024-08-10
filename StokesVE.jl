using Base: throw_invalid_char
function form_vep_stokes()
    # looping over j
    for j in 1:nx
        # looping over i
        for i in 1:ny

            dxx = (grid.x[j+1] - grid[j-1])
            dyy = (grid.y[i] - grid.y[i-1])

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
                # Computing Z term from equation 13.3 (See Gerya)
                Z_n_f = dt * mu_n(i,j+1) / (dt * mu_n(i,j+1) + eta_vp_n(i,j+1))
                Z_n = dt * mu_n(i,j) / (dt * mu_n(i,j) + eta_vp_n(i,j))
                Z_s = dt * mu_s(i,j) / (dt * mu_s(i,j) + eta_vp_s(i,j))
                Z_s_b = dt * mu_s(i-1,j) / (dt * mu_s(i-1,j) + eta_vp_s(i-1,j))

                # Computing left-hand side terms
                L_term1 = 4/dxx
                L_term2 = (eta_vp_n(i,j+1) * exx_dot(i,j+1) * Z_n_f) - (eta_vp_n(i,j) * exx_dot(i,j) * Z_n)
                L_term3 = 2/dyy
                L_term4 = (eta_vp_s(i,j) * exx_dot(i,j) * Z_s) - (eta_vp_s(i-1,j) * exy_dot(i-1,j) * Z_s_b)
                L_term5 = -2 * (P(i,j+1) - P(i,j)) / dxx
                L_term6 = -gx * dt
                L_term7 = (vx(i,j) * (rho_h(i,j+1) - rho_h(i,j-1))/dxx) +
                    (vy(i-1,j) + vy(i-1,j+1) + vy(i-1,j)) * (rho_h(i+1,j) -
                    rho_h(i-1,j))/(2*(grid.y[i+1] + grid.y[i] - grid.y[i-1] - grid.y[i-2]))

                # Computing right-hand side terms
                R_term1 = -rho_h(i,j) * gx
                R_term2 = -2/dxx
                R_term3 = (sigma_xx_o(i,j+1) * eta_vp_n(i,j+1)) / (mu_n(i,j+1) * dt + eta_vp_n(i,j+1))
                R_term4 = -(sigma_xx_o(i,j) * eta_vp_n(i,j)) / (mu_n(i,j) * dt + eta_vp_n(i,j))
                R_term5 = -1/dyy
                R_term6 = (sigma_xy_o(i,j) * eta_vp_s(i,j)) / (mu_s(i,j) * dt + eta_vp_s(i,j))
                R_term7 = -(sigma_xy_o(i-1,j) * eta_vp_s(i-1,j)) / (mu_s(i-1,j) * dt + eta_vp_s(i-1,j))

                # Left-hand side
                
            end
        end
    end
    @views row_index = row_index[1:(k-1)]
    @views col_index = col_index[1:(k-1)]
    @views value = value[1:(k-1)]
    L = sparse(row_index,col_index,value)
    return L,R
end
