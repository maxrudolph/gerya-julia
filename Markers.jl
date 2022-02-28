mutable struct Markers    
    x::Array{Float64,2}
    cell::Array{Int,2}
    scalarFields::Dict
    scalars::Array{Float64,2}
    integerFields::Dict
    integers::Array{Int16,2} # note - this could be changed if larger numbers need to be stored...

    nmark::Int64
    
    function Markers(grid::CartesianGrid,scalarFieldNames,integerFieldNames; nmx::Integer=5,nmy::Integer=5,random::Bool=false)
        N = nmx*nmy*grid.nx*grid.ny
        mdx = grid.W/nmx/(grid.nx-1)
        mdy = grid.H/nmy/(grid.ny-1)
        
        n_fields = length(scalarFieldNames)
        scalarFields = Dict()
        ind=1
        for field in scalarFieldNames
           scalarFields[field] = ind
           ind += 1
        end
        
        n_ifields = length(integerFieldNames)
        integerFields = Dict()
        ind=1
        for field in integerFieldNames
            integerFields[field] = ind
            ind += 1
        end
        
        x = Array{Float64,2}(undef,2,N)
        cell = Array{Int,2}(undef,2,N)
        
        scalars = Array{Float64,2}(undef,n_fields,N)
        integers = Array{Int16,2}(undef,n_ifields,N)
        
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

        new(x,cell,scalarFields,scalars,integerFields,integers,k-1)
    end
end

function find_cell(x::Float64,gridx::Vector{Float64},nx::Int ; guess::Integer=nothing)
    # find the cell in the array gridx that contains the current marker
    # first, see whether the initial guess is correct.
    lower::Int = 1
    upper::Int = nx
    if guess != nothing && guess >= 1 && guess < nx
        if x >= gridx[guess] && x < gridx[guess+1]
            return guess            
        elseif guess < nx-1 && x>=gridx[guess+1] && x<gridx[guess+2]
            return guess+1
        elseif guess > 1 && x < gridx[guess] && x >= gridx[guess-1]
            return guess-1
        else
            if x>=gridx[guess+1]
                lower = guess
            else
                upper = guess+1
            end
        end
     end
    # locate cell using bisection on lower,upper
    while upper-lower > 1 
        midpoint::Int = lower + floor((upper-lower)/2)
        if x >= gridx[midpoint]
            lower = midpoint
        else
            upper = midpoint
        end
    end
    return lower
end

function find_cells!(markers::Markers,grid::CartesianGrid)
   for i in 1:markers.nmark
        markers.cell[1,i] = find_cell(markers.x[1,i] , grid.x, grid.nx, guess=markers.cell[1,i])
        markers.cell[2,i] = find_cell(markers.x[2,i] , grid.y, grid.ny, guess=markers.cell[2,i])
    end
end

function marker_to_basic_node0(m::Markers,grid::CartesianGrid,markerfield::Array{Float64,1})
     # move quantities from the markers to the basic nodes.
     # currently moves rho and eta.
     # returns rho, eta, each as a ny-by-nx matrix

     weights = zeros(Float64,grid.ny,grid.nx)
     field = zeros(Float64,grid.ny,grid.nx)
     # loop over the markers
     for i in 1:m.nmark
        # calculate weights for four surrounding basic nodes
          cellx::Int = m.cell[1,i]
          celly::Int = m.cell[2,i]
          wx = (m.x[1,i] - grid.x[cellx])/(grid.x[cellx+1]-grid.x[cellx]) # mdx/dx
          wy = (m.x[2,i] - grid.y[celly])/(grid.y[celly+1]-grid.y[celly])
          #i,j
          wt_i_j=(1.0-wx)*(1.0-wy)
          #i+1,j        
          wt_i1_j = (1.0-wx)*(wy)
          #i,j+1
          wt_i_j1 = (wx)*(1.0-wy)
          #i+1,j+1
          wt_i1_j1 = (wx)*(wy)


          field[celly,cellx] += wt_i_j*markerfield[i]
          field[celly+1,cellx] += wt_i1_j*markerfield[i]
          field[celly,cellx+1] += wt_i_j1*markerfield[i]
          field[celly+1,cellx+1] += wt_i1_j1*markerfield[i]

          weights[celly,cellx] += wt_i_j
          weights[celly+1,cellx] += wt_i1_j
          weights[celly,cellx+1] += wt_i_j1
          weights[celly+1,cellx+1] += wt_i1_j1       
     end

     field = field ./ weights

     return field
 end

function marker_to_basic_node(m::Markers,grid::CartesianGrid,markerfield::Array{Float64})
    # move quantities from the markers to the basic nodes.
    # currently moves rho and eta.
    # returns rho, eta, each as a ny-by-nx matrix
    if size(markerfield,1) > 1 && size(markerfield,2) == 1
        markerfield = reshape(markerfield,1,size(markerfield,1))
        nfield=1
    else
        nfield = size(markerfield,1)
    end
    
    # loop over the markers
    row = Vector{Int64}(undef,4*m.nmark)
    col = Vector{Int64}(undef,4*m.nmark)
    val_wt = Vector{Float64}(undef,4*m.nmark)
    val_field = Array{Float64,2}(undef,nfield,4*m.nmark)
    
    # idea - create an N-x-nmark sparse matrix. Compute the weights by taking the rowsum    
    Threads.@threads for i in 1:m.nmark
       # calculate weights for four surrounding basic nodes
         local cellx::Int64 = m.cell[1,i]
         local celly::Int64 = m.cell[2,i]
         local wx::Float64 = (m.x[1,i] - grid.x[cellx])/(grid.x[cellx+1]-grid.x[cellx]) # mdx/dx
         local wy::Float64 = (m.x[2,i] - grid.y[celly])/(grid.y[celly+1]-grid.y[celly])
         #i,j
         local wt_i_j::Float64=(1.0-wx)*(1.0-wy)
         #i+1,j        
         local wt_i1_j::Float64 = (1.0-wx)*(wy)
         #i,j+1
         local wt_i_j1::Float64 = (wx)*(1.0-wy)
         #i+1,j+1
         local wt_i1_j1::Float64 = (wx)*(wy)

         local ind::Int64 = 4*(i-1) + 1
        row[ind] = ind
        col[ind] = node_index(celly,cellx,grid.ny)
        val_wt[ind] = wt_i_j

        row[ind+1] = ind
        col[ind+1] = node_index(celly+1,cellx,grid.ny)
        val_wt[ind+1] = wt_i1_j

        row[ind+2] = ind
        col[ind+2] = node_index(celly,cellx+1,grid.ny)
        val_wt[ind+2] = wt_i_j1

        row[ind+3] = ind
        col[ind+3] = node_index(celly+1,cellx+1,grid.ny)
        val_wt[ind+3] = wt_i1_j1
        for k in 1:nfield
            val_field[k,ind] = wt_i_j*markerfield[k,i]
            val_field[k,ind+1] = wt_i1_j*markerfield[k,i]
            val_field[k,ind+2] = wt_i_j1*markerfield[k,i]
            val_field[k,ind+3] = wt_i1_j1*markerfield[k,i]    
        end
    end
    fields = [ sum( sparse(row,col,@views val_field[k,:]),dims=1) for k in 1:nfield ];
    weights = sparse(row,col,val_wt)

    wtsum = sum(weights,dims=1)
    for k in 1:nfield
        fields[k] = reshape(fields[k] ./ wtsum,grid.ny,grid.nx)
    end
    if nfield == 1
        return fields[1]
    else
        return fields
    end
end

function marker_to_cell_center(m::Markers,grid::CartesianGrid,fieldnames)
    # move quantities from the markers to the basic nodes.
    # currently moves rho and eta.
    # returns rho, eta, each as a ny-by-nx matrix
    nfields = length(fieldnames)
    # markerfields will be indices into the 'scalars' array
    markerfields = [m.scalarFields[tmp] for tmp in fieldnames]
    
    weights = zeros(Float64,grid.ny+1,grid.nx+1)
    field = zeros(Float64,nfields,grid.ny+1,grid.nx+1)
    # loop over the markers
    for i in 1:m.nmark
       # calculate weights for four surrounding basic nodes
         cellx::Int =  m.cell[1,i]
         cellx += cellx < grid.nx && m.x[1,i] >= grid.xc[cellx+1] ? 1 : 0
         celly::Int = m.cell[2,i]
         celly += celly < grid.ny && m.x[2,i] >= grid.yc[celly+1] ? 1 : 0
        
         wx = (m.x[1,i] - grid.xc[cellx])/(grid.xc[cellx+1]-grid.xc[cellx]) # mdx/dx
         wy = (m.x[2,i] - grid.yc[celly])/(grid.yc[celly+1]-grid.yc[celly])
         #i,j
         wt_i_j=(1.0-wx)*(1.0-wy)
         #i+1,j        
         wt_i1_j = (1.0-wx)*(wy)
         #i,j+1
         wt_i_j1 = (wx)*(1.0-wy)
         #i+1,j+1
         wt_i1_j1 = (wx)*(wy)
        
         for k in 1:nfields
             kp = markerfields[k]
             field[k,celly,cellx] += wt_i_j*m.scalars[kp,i]
             field[k,celly+1,cellx] += wt_i1_j*m.scalars[kp,i]
             field[k,celly,cellx+1] += wt_i_j1*m.scalars[kp,i]
             field[k,celly+1,cellx+1] += wt_i1_j1*m.scalars[kp,i]
         end
         weights[celly,cellx] += wt_i_j
         weights[celly+1,cellx] += wt_i1_j
         weights[celly,cellx+1] += wt_i_j1
         weights[celly+1,cellx+1] += wt_i1_j1       
    end
    for k in 1:nfields
        field[k,:,:] = field[k,:,:] ./ weights
    end
    
    return field
end

function marker_to_cell_center(m::Markers,grid::CartesianGrid,markerfield::Array{Float64})
    # markerfields will be indices into the 'scalars' array
    
    if size(markerfield,1) > 1 && size(markerfield,2) == 1
        markerfield = reshape(markerfield,1,size(markerfield,1))
        nfield=1
    else
        nfield = size(markerfield,1)
    end
    
    # loop over the markers
    row = Vector{Int64}(undef,4*m.nmark)
    col = Vector{Int64}(undef,4*m.nmark)
    val_wt = Vector{Float64}(undef,4*m.nmark)
    val_field = Array{Float64,2}(undef,nfield,4*m.nmark)
    
    # loop over the markers
    Threads.@threads for i in 1:m.nmark
       # calculate weights for four surrounding basic nodes
         cellx::Int =  m.cell[1,i]
         cellx += cellx < grid.nx && m.x[1,i] >= grid.xc[cellx+1] ? 1 : 0
         celly::Int = m.cell[2,i]
         celly += celly < grid.ny && m.x[2,i] >= grid.yc[celly+1] ? 1 : 0
        
         local wx::Float64 = (m.x[1,i] - grid.xc[cellx])/(grid.xc[cellx+1]-grid.xc[cellx]) # mdx/dx
         local wy::Float64 = (m.x[2,i] - grid.yc[celly])/(grid.yc[celly+1]-grid.yc[celly])
         #i,j
         local wt_i_j::Float64=(1.0-wx)*(1.0-wy)
         #i+1,j        
         local wt_i1_j::Float64 = (1.0-wx)*(wy)
         #i,j+1
         local wt_i_j1::Float64 = (wx)*(1.0-wy)
         #i+1,j+1
         local wt_i1_j1::Float64 = (wx)*(wy)
        
         local ind::Int64 = 4*(i-1) + 1
         local cell_ind::Int64 = celly + (cellx-1)*(grid.ny+1)
        
         row[ind] = ind
        col[ind] = cell_ind; #node_index(celly,cellx,grid.ny)
        val_wt[ind] = wt_i_j

        row[ind+1] = ind
        col[ind+1] = cell_ind + 1 #node_index(celly+1,cellx,grid.ny)
        val_wt[ind+1] = wt_i1_j

        row[ind+2] = ind
        col[ind+2] = cell_ind+grid.ny+1 # node_index(celly,cellx+1,grid.ny)
        val_wt[ind+2] = wt_i_j1

        row[ind+3] = ind
        col[ind+3] = cell_ind+grid.ny+2 # node_index(celly+1,cellx+1,grid.ny)
        val_wt[ind+3] = wt_i1_j1
        for k in 1:nfield
            val_field[k,ind] = wt_i_j*markerfield[k,i]
            val_field[k,ind+1] = wt_i1_j*markerfield[k,i]
            val_field[k,ind+2] = wt_i_j1*markerfield[k,i]
            val_field[k,ind+3] = wt_i1_j1*markerfield[k,i]    
        end                   
    end
    
    fields = [ sum( sparse(row,col,@views val_field[k,:]),dims=1) for k in 1:nfield ];
    weights = sparse(row,col,val_wt)

    wtsum = sum(weights,dims=1)
    for k in 1:nfield
        fields[k] = reshape(fields[k] ./ wtsum,grid.ny+1,grid.nx+1)
    end
    if nfield == 1
        return fields[1]
    else
        return fields
    end
end


function marker_to_cell_center0(m::Markers,grid::CartesianGrid,markerfield::Array{Float64,1})
    # markerfields will be indices into the 'scalars' array
    weights = zeros(Float64,grid.ny+1,grid.nx+1)
    field = zeros(Float64,grid.ny+1,grid.nx+1)
    # loop over the markers
    for i in 1:m.nmark
       # calculate weights for four surrounding basic nodes
         cellx::Int =  m.cell[1,i]
         cellx += cellx < grid.nx && m.x[1,i] >= grid.xc[cellx+1] ? 1 : 0
         celly::Int = m.cell[2,i]
         celly += celly < grid.ny && m.x[2,i] >= grid.yc[celly+1] ? 1 : 0
        
         wx = (m.x[1,i] - grid.xc[cellx])/(grid.xc[cellx+1]-grid.xc[cellx]) # mdx/dx
         wy = (m.x[2,i] - grid.yc[celly])/(grid.yc[celly+1]-grid.yc[celly])
         #i,j
         wt_i_j=(1.0-wx)*(1.0-wy)
         #i+1,j        
         wt_i1_j = (1.0-wx)*(wy)
         #i,j+1
         wt_i_j1 = (wx)*(1.0-wy)
         #i+1,j+1
         wt_i1_j1 = (wx)*(wy)
        
         field[celly,cellx] += wt_i_j*markerfield[i]
         field[celly+1,cellx] += wt_i1_j*markerfield[i]
         field[celly,cellx+1] += wt_i_j1*markerfield[i]
         field[celly+1,cellx+1] += wt_i1_j1*markerfield[i]

         weights[celly,cellx] += wt_i_j
         weights[celly+1,cellx] += wt_i1_j
         weights[celly,cellx+1] += wt_i_j1
         weights[celly+1,cellx+1] += wt_i1_j1       
    end
    field = field ./ weights    
    return field
end

function marker_to_basic_node(m::Markers,grid::CartesianGrid,fieldnames)
    # move quantities from the markers to the basic nodes.
    # This is a convenience function that allows the field names to be passed as a list of strings.
    nfields = length(fieldnames)
    # markerfields will be indices into the 'scalars' array
    markerfields = [m.scalarFields[tmp] for tmp in fieldnames]
    
    fields = marker_to_basic_node(m,grid, m.scalars[markerfields,:] )
    return fields
end

function basic_node_to_markers!(m::Markers,grid::CartesianGrid,field::Matrix,mfield::String)
    k = m.scalarFields[mfield]
    Threads.@threads for i in 1:m.nmark
        cellx = m.cell[1,i]
        celly = m.cell[2,i]
        wx::Float64 = (m.x[1,i] - grid.x[cellx])/(grid.x[cellx+1]-grid.x[cellx]) # mdx/dx
        wy::Float64 = (m.x[2,i] - grid.y[celly])/(grid.y[celly+1]-grid.y[celly])
        
        m.scalars[k,i] = (1.0-wx)*(1.0-wy)*field[celly,cellx] +
            + (wx)*(1.0-wy)*field[celly,cellx+1] +
            + (1.0-wx)*(wy)*field[celly+1,cellx] +
            + (wx)*(wy)*field[celly+1,cellx+1]
    end
end

function basic_node_to_markers!(m::Markers,grid::CartesianGrid,field::Matrix,mfield::Array{Float64,1})
    Threads.@threads for i in 1:m.nmark
        cellx = m.cell[1,i]
        celly = m.cell[2,i]
        wx::Float64 = (m.x[1,i] - grid.x[cellx])/(grid.x[cellx+1]-grid.x[cellx]) # mdx/dx
        wy::Float64 = (m.x[2,i] - grid.y[celly])/(grid.y[celly+1]-grid.y[celly])
        
        mfield[i] = (1.0-wx)*(1.0-wy)*field[celly,cellx] +
            + (wx)*(1.0-wy)*field[celly,cellx+1] +
            + (1.0-wx)*(wy)*field[celly+1,cellx] +
            + (wx)*(wy)*field[celly+1,cellx+1]
    end
end

function cell_center_to_markers!(m::Markers,grid::CartesianGrid,field::Matrix,mfield::Array{Float64,1})
    if size(field,1) == grid.nx+1
        cellx_max = grid.nx
    else
        cellx_max = grid.nx-1
    end
    if size(field,2) == grid.ny+1
        celly_max = grid.ny
    else
        celly_max = grid.ny-1
    end
    
    Threads.@threads for i in 1:m.nmark
        cellx = m.cell[1,i]
        celly = m.cell[2,i]
        
        cellx += cellx < cellx_max && m.x[1,i] >= grid.xc[cellx+1] ? 1 : 0
        celly::Int = m.cell[2,i]
        celly += celly < celly_max && m.x[2,i] >= grid.yc[celly+1] ? 1 : 0
        
        wx::Float64 = (m.x[1,i] - grid.xc[cellx])/(grid.xc[cellx+1]-grid.xc[cellx]) # mdx/dx
        wy::Float64 = (m.x[2,i] - grid.yc[celly])/(grid.yc[celly+1]-grid.yc[celly])
        
        mfield[i] = (1.0-wx)*(1.0-wy)*field[celly,cellx] +
            + (wx)*(1.0-wy)*field[celly,cellx+1] +
            + (1.0-wx)*(wy)*field[celly+1,cellx] +
            + (wx)*(wy)*field[celly+1,cellx+1]
    end
end

function cell_center_change_to_markers!(m::Markers,grid::CartesianGrid,field::Matrix,mfield::String)
    if size(field,1) == grid.nx+1
        cellx_max = grid.nx
    else
        cellx_max = grid.nx-1
    end
    if size(field,2) == grid.ny+1
        celly_max = grid.ny
    else
        celly_max = grid.ny-1
    end
    k = m.scalarFields[mfield]
    Threads.@threads for i in 1:m.nmark
        cellx = m.cell[1,i]
        celly = m.cell[2,i]
        
        cellx += cellx < cellx_max && m.x[1,i] >= grid.xc[cellx+1] ? 1 : 0
         celly::Int = m.cell[2,i]
         celly += celly < celly_max && m.x[2,i] >= grid.yc[celly+1] ? 1 : 0
        
        wx::Float64 = (m.x[1,i] - grid.xc[cellx])/(grid.xc[cellx+1]-grid.xc[cellx]) # mdx/dx
        wy::Float64 = (m.x[2,i] - grid.yc[celly])/(grid.yc[celly+1]-grid.yc[celly])
        
        m.scalars[k,i] += (1.0-wx)*(1.0-wy)*field[celly,cellx] +
            + (wx)*(1.0-wy)*field[celly,cellx+1] +
            + (1.0-wx)*(wy)*field[celly+1,cellx] +
            + (wx)*(wy)*field[celly+1,cellx+1]
    end
end

function basic_node_change_to_markers!(m::Markers,grid::CartesianGrid,field::Matrix,mfield::String)
    k = m.scalarFields[mfield]
    Threads.@threads for i in 1:m.nmark
        cellx = m.cell[1,i]
        celly = m.cell[2,i]
        wx::Float64 = (m.x[1,i] - grid.x[cellx])/(grid.x[cellx+1]-grid.x[cellx]) # mdx/dx
        wy::Float64 = (m.x[2,i] - grid.y[celly])/(grid.y[celly+1]-grid.y[celly])
        
        m.scalars[k,i] += (1.0-wx)*(1.0-wy)*field[celly,cellx] +
            + (wx)*(1.0-wy)*field[celly,cellx+1] +
            + (1.0-wx)*(wy)*field[celly+1,cellx] +
            + (wx)*(wy)*field[celly+1,cellx+1]
    end
end

# function basic_node_to_markers!(m::Markers,grid::CartesianGrid,field::Matrix)
#     Threads.@threads for i in 1:m.nmark
#         cellx = m.cell[1,i]
#         celly = m.cell[2,i]
#         wx::Float64 = (m.x[1,i] - grid.x[cellx])/(grid.x[cellx+1]-grid.x[cellx]) # mdx/dx
#         wy::Float64 = (m.x[2,i] - grid.y[celly])/(grid.y[celly+1]-grid.y[celly])
        
#         m.rho[i] = (1.0-wx)*(1.0-wy)*field[celly,cellx] +
#             + (wx)*(1.0-wy)*field[celly,cellx+1] +
#             + (1.0-wx)*(wy)*field[celly+1,cellx] +
#             + (wx)*(wy)*field[celly+1,cellx+1]
#     end
# end

function viscosity_to_cell_centers(grid::CartesianGrid,etas::Matrix)
    # compute the harmonic average of the viscosities at the nodal points
   etan = zeros(grid.ny,grid.nx)
    for i in 2:grid.ny
        for j in 2:grid.nx
            etan[i,j] = 1/( (1/etas[i-1,j-1] + 1/etas[i-1,j] + 1/etas[i,j-1] + 1/etas[i,j])/4. )
        end
    end
    return etan
end

function velocity_to_centers(grid::CartesianGrid,vx::Matrix,vy::Matrix)
    # compute vx and vy at cell centers
     vxc = zeros(grid.ny+1,grid.nx+1);
     vyc = zeros(grid.ny+1,grid.nx+1);
     # cell centers are offset in (-) direction from basic nodes.
     #               |
     # (center)     vx[i,j]
     #               |
     # ---vy[i,j]---(i,j)       
     for i in 2:grid.ny        # interior...
        for j in 2:grid.nx
            # left
            vxm = vx[i,j-1] # this will produce vx=0 along the left boundary
            vxp = vx[i,j]
            # top
            vym = vy[i-1,j] # vy=0 along the top boundary
            vyp = vy[i,j]
            vxc[i,j] = 0.5*(vxp+vxm)
            vyc[i,j] = 0.5*(vyp+vym)            
        end
    end
    # vx - top
    vxc[1,2:grid.nx] = vxc[2,2:grid.nx]
    # bottom
    vxc[grid.ny+1,2:grid.nx] = vxc[grid.ny,2:grid.nx]
    # left
    vxc[:,1] = -vxc[:,2]
    # right
    vxc[:,grid.nx+1] = - vxc[:,grid.nx]

    # vy - left
    vyc[2:grid.ny,1] = vyc[2:grid.ny,2]
    # vy - right
    vyc[2:grid.ny,grid.nx+1] = vyc[2:grid.ny,grid.nx]
    # vy - top
    vyc[1,:] = -vyc[2,:]
    # vy - bottom
    vyc[grid.ny+1,:] = -vyc[grid.ny,:]        
    
    return vxc,vyc
end

function velocity_to_basic_nodes(grid::CartesianGrid,vxc::Matrix,vyc::Matrix)
    # this gets the velocity in a format suitable for visualization.
    # NOTE - performs a transpose on the grid!!!
    vn = Array{Float64,3}(undef,2,grid.nx,grid.ny)
    for i in 1:grid.ny
        for j in 1:grid.nx
            vn[1,j,i] = 0.25*(vxc[i,j]+vxc[i+1,j]+vxc[i,j+1]+vxc[i+1,j+1])
            vn[2,j,i] = 0.25*(vyc[i,j]+vyc[i+1,j]+vyc[i,j+1]+vyc[i+1,j+1])
        end
    end
    return vn
end

function velocity_to_points(x::Matrix,cell::Matrix,grid::CartesianGrid,vxc::Matrix,vyc::Matrix;N::Int=-1)
    # compute velocity at N points given in x (2-by-N)
    # cell should contain the cells in which the points are located (2-by-N)
    # this routine assumes that the velocity (vxc and vyc) is defined at the cell centers
    if N==-1
        N = size(x,2)
    end
    mvx = Array{Float64,1}(undef,N) # velocities at specified locations
    mvy = Array{Float64,1}(undef,N) 
    Threads.@threads for i in 1:N
        cellx::Int = x[1,i] < grid.xc[cell[1,i]+1] ? cell[1,i] : cell[1,i] + 1
        celly::Int = x[2,i] < grid.yc[cell[2,i]+1] ? cell[2,i] : cell[2,i] + 1
        mdx::Float64 = (x[1,i] - grid.xc[cellx])/(grid.xc[cellx+1]-grid.xc[cellx])
        mdy::Float64 = (x[2,i] - grid.yc[celly])/(grid.yc[celly+1]-grid.yc[celly])
        mvx[i] = (1-mdx)*(1-mdy)*vxc[celly,cellx] +
            + (mdx)*(1-mdy)*vxc[celly,cellx+1] +
            + (1-mdx)*(mdy)*vxc[celly+1,cellx] +
            + (mdx)*(mdy)*vxc[celly+1,cellx+1]
        mvy[i] = (1-mdx)*(1-mdy)*vyc[celly,cellx] +
            + (mdx)*(1-mdy)*vyc[celly,cellx+1] +
            + (1-mdx)*(mdy)*vyc[celly+1,cellx] +
            + (mdx)*(mdy)*vyc[celly+1,cellx+1]
    end    
    return mvx,mvy
end

function velocity_to_markers(m::Markers,grid::CartesianGrid,vxc::Matrix,vyc::Matrix)
    # This function expects the velocities to be defined at the cell centers. vxc and vyc should each have
    # an 'extra' column and row corresponding to the ghost degrees of freedom that are needed to interpolate
    # velocities along the bottom and left of the domain.
    mvx,mvy = velocity_to_points(m.x,m.cell,grid,vxc,vyc;N=m.nmark)
    return mvx,mvy
end

function move_markers!(markers::Markers,grid::CartesianGrid,vxc::Matrix,vyc::Matrix,dt::Float64)
    # move the markers using the 1st-order algorithm (forward Euler)
    mvx,mvy = velocity_to_markers(markers,grid,vxc,vyc)
    # determine the maximal timestep
    vxmax = maximum(abs.(mvx))
    vymax = maximum(abs.(mvy))
    Threads.@threads for i in 1:markers.nmark
        markers.x[1,i] += dt*mvx[i]
        markers.x[2,i] += dt*mvy[i]
    end    
    find_cells!(markers,grid)
    return dt
end

function move_markers_rk2!(markers::Markers,grid::CartesianGrid,vxc::Matrix,vyc::Matrix,dt::Float64)
    # move the markers using the 2nd-order Runge-Kutta algorithm.
    # compute velocities for each marker at current position
    mvx,mvy = velocity_to_markers(markers,grid,vxc,vyc)
    # compute marker location at xA, xB
    xB = Array{Float64,2}(undef,2,markers.nmark)
    Threads.@threads for i in 1:markers.nmark
        xB[1,i] = markers.x[1,i] + dt/2*mvx[i]
        xB[2,i] = markers.x[2,i] + dt/2*mvy[i]
    end
    # re-locate markers, which may now be in a different cell.
    cell = copy(markers.cell)
    Threads.@threads for i in 1:markers.nmark
        cell[1,i] = find_cell(xB[1,i], grid.x, grid.nx, guess=cell[1,i])
        cell[2,i] = find_cell(xB[2,i], grid.y, grid.ny, guess=cell[2,i])
    end
    # compute velocity at xB
    mvx,mvy = velocity_to_points(xB,cell,grid,vxc,vyc)
    # Move the markers using the velocity at xB.
    Threads.@threads for i in 1:markers.nmark
        markers.x[1,i] += dt*mvx[i]
        markers.x[2,i] += dt*mvy[i]
    end
    # re-locate markers in their new cells.
    find_cells!(markers,grid)
end



