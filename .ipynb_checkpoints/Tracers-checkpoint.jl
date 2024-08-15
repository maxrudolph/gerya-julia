mutable struct Tracers    
    x::Array{Float64,2}
    cell::Array{Int64,2}
    #scalarFields::Dict
    #scalars::Array{Float64,2}
    #integerFields::Dict
    #integers::Array{Int16,3} # note - this could be changed if larger numbers need to be stored...

    nmark::Int64
    max_mark::Int64
    
    function Tracers(grid::CartesianGrid, nmx::Integer=5,nmy::Integer=5,random::Bool=false)
        # scalarFieldNames,integerFieldNames;
        N = nmx*nmy*(grid.nx-1)*(grid.ny-1) # total number of markers
        Nmax = Int(ceil( N*1.20 ))
        mdx = grid.W/nmx/(grid.nx-1)
        mdy = grid.H/nmy/(grid.ny-1)
        
        #n_fields = length(scalarFieldNames)
        #scalarFields = Dict()
        #ind=1
        #for field in scalarFieldNames
        #   scalarFields[field] = ind
        #   ind += 1
        #end
        
        #n_ifields = length(integerFieldNames)
        #integerFields = Dict()
        #ind=1
        #for field in integerFieldNames
        #    integerFields[field] = ind
        #    ind += 1
        #end
        
        x = Array{Float64,2}(undef,2,Nmax)
        cell = Array{Int64,2}(undef,2,Nmax)
        
        #scalars = Array{Float64,2}(undef,n_fields,Nmax)
        #integers = Array{Int16,4}(undef,n_ifields,Nmax)
        
        k=1
        for i in 1:(grid.ny-1)
            for j in 1:(grid.nx-1)
                for ii in 1:nmy
                     for jj in 1:nmx
                        x[1,k] = mdx/2. + mdx*(jj-1) + mdx*nmx*(j-1) + ( random ? (rand()-0.5)*mdx : 0.0 )
                        x[2,k] = mdy/2. + mdy*(ii-1) + mdy*nmy*(i-1) + ( random ? (rand()-0.5)*mdy : 0.0 )
                        cell[1,k] = j
                        cell[2,k] = i
                        k+=1
                    end
                end
            end
        end

        new(x,cell,k-1,Nmax)
    end
end

function velocity_to_markers(m::Markers,grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64};vxc=nothing,vyc=nothing,continuity_weight::Float64=0.0)
    # This function expects the velocities to be defined at the cell centers. vxc and vyc should each have
    # an 'extra' column and row corresponding to the ghost degrees of freedom that are needed to interpolate
    # velocities along the bottom and left of the domain.
    mvx,mvy = velocity_to_points(m.x,m.cell,grid,vx,vy;continuity_weight=continuity_weight,N=m.nmark)
    return mvx,mvy
end
function velocity_to_points(x::Matrix{Float64},cell::Matrix{Int64},grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64}; continuity_weight::Float64=0.0,N::Int64=-1,vxc=nothing,vyc=nothing)
    # compute the velocity at points, using the continuity-based velocity interpolation
    # compute velocity at cell centers
    # compute the velocity at the markers from the velocity nodes:
    mvx,mvy = velocity_nodes_to_points(x,cell,grid,vx,vy,N=N)
    if continuity_weight == 0.0
        return mvx,mvy
    end

    # compute velocity at cell centers
    if vxc == nothing || vyc == nothing
        vxc,vyc = velocity_to_centers(grid,vx,vy);
    end
    # compute the velocity at the markers from the cell centers:
    mvxc,mvyc = velocity_center_to_points(x,cell,grid,vxc,vyc,N=N)

    mvx = continuity_weight*(mvxc) + (1.0-continuity_weight)*(mvx)
    mvy = continuity_weight*(mvyc) + (1.0-continuity_weight)*(mvy)
    return mvx,mvy
end
function move_markers_rk4!(markers::Markers,grid::CartesianGrid,vx::Matrix{Float64},vy::Matrix{Float64},dt::Float64; continuity_weight::Float64=1.0/3.0)
    # This function implements the 4th-order Runge-Kutta scheme for advection of markers. It expects
    # vx and vy are the velocities at the velocity nodes
    # dt is the timestep
    if continuity_weight != 0.0
        vxc,vyc = velocity_to_centers(grid,vx,vy)
    else
        vxc=nothing
        vyc=nothing
    end
    # 1. compute velocity at point A
    vxA::Vector{Float64}, vyA::Vector{Float64} = velocity_to_markers(markers,grid,vx,vy,continuity_weight=continuity_weight,vxc=vxc,vyc=vyc)
    # 2. compute xB=xA + vA*dt/2
    xB = Array{Float64,2}(undef,2,markers.nmark)
    for i in 1:markers.nmark
        xB[1,i] = markers.x[1,i] + dt/2*vxA[i]
        xB[2,i] = markers.x[2,i] + dt/2*vyA[i]
    end
    # 3. locate xB and compute vxB
    cell::Matrix{Int64} = copy(markers.cell[:,1:markers.nmark])
    Threads.@threads for i in 1:markers.nmark
        cell[1,i] = find_cell(xB[1,i], grid.x, grid.nx, guess=cell[1,i])
        cell[2,i] = find_cell(xB[2,i], grid.y, grid.ny, guess=cell[2,i])
    end
    vxB, vyB = velocity_to_points(xB,cell,grid,vx,vy,continuity_weight=continuity_weight,vxc=vxc,vyc=vyc,N=markers.nmark)
    # 4. compute xC = xA+vB*dt/2
    xC = Array{Float64,2}(undef,2,markers.nmark)
    for i in 1:markers.nmark
        xC[1,i] = markers.x[1,i] + dt/2*vxB[i]
        xC[2,i] = markers.x[2,i] + dt/2*vyB[i]
    end
    # 5. locate cells for xC and compute vC
    Threads.@threads for i in 1:markers.nmark
        cell[1,i] = find_cell(xC[1,i], grid.x, grid.nx, guess=cell[1,i])
        cell[2,i] = find_cell(xC[2,i], grid.y, grid.ny, guess=cell[2,i])
    end
    vxC, vyC = velocity_to_points(xC,cell,grid,vx,vy,continuity_weight=continuity_weight,vxc=vxc,vyc=vyc,N=markers.nmark)
    # 6. compute xD = xA + vC*dt
    xD = Array{Float64,2}(undef,2,markers.nmark)
    for i in 1:markers.nmark
        xD[1,i] = markers.x[1,i] + dt*vxC[i]
        xD[2,i] = markers.x[2,i] + dt*vyC[i]
    end
    # 7. locate cells for xD and compute vD
    Threads.@threads for i in 1:markers.nmark
        cell[1,i] = find_cell(xD[1,i], grid.x, grid.nx, guess=cell[1,i])
        cell[2,i] = find_cell(xD[2,i], grid.y, grid.ny, guess=cell[2,i])
    end
    vxD, vyD = velocity_to_points(xD,cell,grid,vx,vy,continuity_weight=continuity_weight,vxc=vxc,vyc=vyc,N=markers.nmark)
    # 8. Compute v_eff = 1/6*(vA+2*vB+2*vC+vD) and move markers by v_eff*dt
    Threads.@threads for i in 1:markers.nmark
        markers.x[1,i] += dt/6.0*(vxA[i] + 2*vxB[i] + 2*vxC[i] + vxD[i]) 
        markers.x[2,i] += dt/6.0*(vyA[i] + 2*vyB[i] + 2*vyC[i] + vyD[i])
    end
    # 9. relocate markers in their cells.
    find_cells!(markers,grid)
end
