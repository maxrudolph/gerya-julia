# Unit tests for marker removal (add_remove_markers!) and stag_to_points.
# Run from the repository root:  julia --project test/test_markers.jl
# (needs Test, Random and NearestNeighbors)
using Test, Random, NearestNeighbors
include(joinpath(@__DIR__,"..","Grid.jl"))
include(joinpath(@__DIR__,"..","Markers.jl"))

# Independent reference: bilinear interpolation of a staggered field at point (xp,yp).
# stagx/stagy = -1 means the field lives at cell centers (xc/yc, including ghost nodes),
# 0 means it lives at basic nodes (x/y). The bracketing nodes are found by search,
# not from the marker's cell, so this does not share any logic with stag_to_points.
function reference_interp(grid,field,stagx,stagy,xp,yp)
    xn = stagx == -1 ? grid.xc : grid.x
    yn = stagy == -1 ? grid.yc : grid.y
    j = clamp(searchsortedlast(xn,xp),1,size(field,2)-1)
    i = clamp(searchsortedlast(yn,yp),1,size(field,1)-1)
    wx = (xp-xn[j])/(xn[j+1]-xn[j])
    wy = (yp-yn[i])/(yn[i+1]-yn[i])
    return (1-wx)*(1-wy)*field[i,j] + wx*(1-wy)*field[i,j+1] +
           (1-wx)*wy*field[i+1,j] + wx*wy*field[i+1,j+1]
end

function points_and_cells(grid,x)
    cell = similar(x,Int64)
    for k in 1:size(x,2)
        cell[1,k] = find_cell(x[1,k],grid.x,grid.nx,guess=1)
        cell[2,k] = find_cell(x[2,k],grid.y,grid.ny,guess=1)
    end
    return cell
end

# Build a marker set with a unique id per marker, then put `extra[c]` additional markers
# into cell c = (row,col). Returns the markers and the starting per-cell count.
function markers_with_overfull_cells(grid,extra;nmx=5,nmy=5)
    m = Markers(grid,["T","id"],["material"];nmx=nmx,nmy=nmy)
    nextra = sum(values(extra))
    nextra + m.nmark > m.max_mark && add_markers!(m, nextra + m.nmark - m.max_mark)
    for ((i,j),n) in extra
        for _ in 1:n
            k = m.nmark += 1
            m.x[1,k] = grid.x[j] + rand()*(grid.x[j+1]-grid.x[j])
            m.x[2,k] = grid.y[i] + rand()*(grid.y[i+1]-grid.y[i])
            m.cell[1,k] = j; m.cell[2,k] = i
        end
    end
    m.scalars[m.scalarFields["id"],1:m.nmark] = 1:m.nmark
    m.scalars[m.scalarFields["T"],1:m.nmark] .= 1.0
    m.integers[m.integerFields["material"],1:m.nmark] .= 1
    return m, nmx*nmy
end

@testset "Markers.jl" begin

@testset "stag_to_points" begin
    Random.seed!(1234)
    # nx != ny in both orientations (where the swapped size checks matter),
    # and nx == ny as a control (where they happen to agree).
    for (nx,ny) in ((7,5),(5,7),(6,6),(12,4))
        grid = CartesianGrid(1.3*(nx-1),0.9*(ny-1),nx,ny)
        dx = grid.x[2]-grid.x[1]; dy = grid.y[2]-grid.y[1]
        # Random points everywhere, plus points in the outer half of the last cell in x and in y
        # (between the last interior cell center and the wall) and in the far corner. In those
        # half-cells the correct stencil uses the ghost-node values.
        npts = 200
        x = [grid.W .* rand(1,npts); grid.H .* rand(1,npts)]
        xedge = [grid.W .- 0.5*dx .* rand(1,npts); grid.H .* rand(1,npts)]
        yedge = [grid.W .* rand(1,npts); grid.H .- 0.5*dy .* rand(1,npts)]
        corner = [grid.W .- 0.5*dx .* rand(1,npts); grid.H .- 0.5*dy .* rand(1,npts)]
        walls = [grid.W 0.0 grid.W 0.5*grid.W; 0.5*grid.H grid.H grid.H 0.0]

        # A field that is not bilinear, so interpolating from the wrong pair of cell
        # centers (linear extrapolation) gives a different answer.
        f(xx,yy) = sin(1.7*xx) + xx^2*yy - 0.3*cos(2.3*yy)
        @testset "cell centers, nx=$nx ny=$ny" begin
            Tc = [f(grid.xc[j],grid.yc[i]) for i in 1:grid.ny+1, j in 1:grid.nx+1]
            for (label,pts) in (("interior",x),("right half-cell",xedge),
                                ("bottom half-cell",yedge),("corner",corner),("walls",walls))
                vals = stag_to_points(pts,points_and_cells(grid,pts),grid,Tc,-1,-1)
                ref = [reference_interp(grid,Tc,-1,-1,pts[1,k],pts[2,k]) for k in 1:size(pts,2)]
                @test vals[1,:] ≈ ref atol=1e-12
            end
            # interpolation reproduces the field at every cell center inside the domain
            pts = [repeat(grid.xc[2:grid.nx],inner=grid.ny-1)'; repeat(grid.yc[2:grid.ny],outer=grid.nx-1)']
            vals = stag_to_points(pts,points_and_cells(grid,pts),grid,Tc,-1,-1)
            @test vals[1,:] ≈ [Tc[i,j] for j in 2:grid.nx for i in 2:grid.ny] atol=1e-12
        end
        @testset "vx-type field (basic nodes in x, centers in y), nx=$nx ny=$ny" begin
            vx = [f(grid.x[j],grid.yc[i]) for i in 1:grid.ny+1, j in 1:grid.nx]
            for pts in (x,yedge,corner)
                vals = stag_to_points(pts,points_and_cells(grid,pts),grid,vx,0,-1)
                ref = [reference_interp(grid,vx,0,-1,pts[1,k],pts[2,k]) for k in 1:size(pts,2)]
                @test vals[1,:] ≈ ref atol=1e-12
            end
        end
    end
end

@testset "add_remove_markers! removal" begin
    Random.seed!(4321)
    grid = CartesianGrid(6.0,4.0,7,5)   # 6 x 4 cells, 25 markers per cell
    min_markers, target, max_markers = 12, 25, 40
    Tc = zeros(grid.ny+1,grid.nx+1)
    # Overfull cells chosen so that "remove target_markers" (the old bug) and
    # "remove per_cell - target_markers" give different results:
    #   41 -> 25 (old: 16), 60 -> 25 (old: 35, still overfull), 100 -> 25 (old: 75)
    # and a cell at exactly max_markers, which must be left alone.
    extra = Dict((1,1)=>16, (2,3)=>35, (4,6)=>75, (3,2)=>15)
    m, n0 = markers_with_overfull_cells(grid,extra)
    pc0 = markers_per_cell(grid,m)
    @test pc0[1,1] == 41 && pc0[2,3] == 60 && pc0[4,6] == 100 && pc0[3,2] == 40

    ids_before = Dict(c => sort(m.scalars[m.scalarFields["id"],
                      findall(k -> m.cell[2,k]==c[1] && m.cell[1,k]==c[2], 1:m.nmark)])
                      for c in Tuple.(CartesianIndices(pc0)))
    nmark0 = m.nmark
    newm = add_remove_markers!(m,grid,Tc,min_markers,target,max_markers)
    pc = markers_per_cell(grid,m)

    @testset "overfull cells are reduced to the target" begin
        @test pc[1,1] == target
        @test pc[2,3] == target
        @test pc[4,6] == target
        @test all(pc .<= max_markers)
    end
    @testset "only the excess is removed" begin
        @test pc[3,2] == 40                               # at max_markers: not overfull
        others = trues(size(pc)); others[1,1] = others[2,3] = others[4,6] = others[3,2] = false
        @test all(pc[others] .== n0)                      # untouched cells keep their markers
        @test m.nmark == nmark0 - (41-target) - (60-target) - (100-target)
        @test !any(newm[1:m.nmark])                       # nothing was added
    end
    @testset "survivors are original markers with their own data" begin
        ids = m.scalars[m.scalarFields["id"],1:m.nmark]
        @test allunique(ids)
        for c in Tuple.(CartesianIndices(pc))
            in_c = findall(k -> m.cell[2,k]==c[1] && m.cell[1,k]==c[2], 1:m.nmark)
            # the surviving ids in each cell are a subset of that cell's original ids,
            # and for cells that were not overfull they are exactly the original ids
            @test issubset(ids[in_c], ids_before[c])
            pc0[c...] <= max_markers && @test sort(ids[in_c]) == ids_before[c]
            # positions still lie in the cell the marker is assigned to
            @test all(grid.x[c[2]] .<= m.x[1,in_c] .<= grid.x[c[2]+1])
            @test all(grid.y[c[1]] .<= m.x[2,in_c] .<= grid.y[c[1]+1])
        end
    end
    @testset "removal is a random subset, not always the same markers" begin
        # repeat the removal from an identical starting state: the surviving set should vary
        survivors = Set{Vector{Float64}}()
        for trial in 1:5
            m2, _ = markers_with_overfull_cells(grid,Dict((2,3)=>35))
            add_remove_markers!(m2,grid,Tc,min_markers,target,max_markers)
            in_c = findall(k -> m2.cell[2,k]==2 && m2.cell[1,k]==3, 1:m2.nmark)
            push!(survivors, sort(m2.scalars[m2.scalarFields["id"],in_c]))
        end
        @test length(survivors) > 1
    end
end

end # Markers.jl
