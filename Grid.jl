struct BoundaryConditions
    # The intent here is that each boundary gets a flag
    # -1 = Free-slip d/dn = 0
    #  1 = prescribed velocity v=const
    # If the flag 1 is passed, the bc value should be set to 0.0 (i.e., dvx/dx = 0.0).
    # other possibilities?
    vx_bc_type::Vector{Int64}
    vx_bc_value::Vector{Float64}
    vy_bc_type::Vector{Int64}
    vy_bc_value::Vector{Float64}
    function BoundaryConditions(vx_bc_type,vx_bc_value,vy_bc_type,vy_bc_value)
        # vx_bc_type [left,right,top,bottom]
        # vx_bc_value [left,right,top,bottom]
        # same for vy [left,right,top,bottom]        
        new(vx_bc_type,vx_bc_value,vy_bc_type,vy_bc_value)
    end
    function BoundaryConditions()
        # free slip everywhere
        new([1,1,-1,-1],[0.0,0.0,0.0,0.0],[-1,-1,1,1],[0.0,0.0,0.0,0.0])
    end
end

struct CartesianGrid
    x::Array{Float64,1}
    y::Array{Float64,1}
    xc::Array{Float64,1}
    yc::Array{Float64,1}
    nx::Int64
    ny::Int64
    W::Float64
    H::Float64
    yperiodic::Bool
    xperiodic::Bool
    #dx::Float64
    #dy::Float64
    function CartesianGrid(W::Float64,H::Float64,nx::Int64,ny::Int64;xstart::Float64=0.0,ystart::Float64=0.0)
        dx = W/(nx-1)
        dy = H/(ny-1)
        new(xstart.+LinRange(0,W,nx),ystart.+LinRange(0,H,ny),
            xstart.+LinRange(0-dx/2,W+dx/2,nx+1),ystart.+LinRange(0-dy/2,H+dy/2,ny+1),
            nx,ny,W,H,false,false)
    end
end

# struct CylindricalGrid
#     r::Array{Float64,1}
#     z::Array{Float64,1}
#     rc::Array{Float64,1}
#     zc::Array{Float64,1}
#     nr::Int64
#     nz::Int64
#     W::Float64
#     H::Float64
#     #dx::Float64
#     #dy::Float64
#     function CylindricalGrid(W::Float64,H::Float64,nr::Int,nz::Int)
#         dr = W/(nr-1)
#         dz = H/(nz-1)
#         new(LinRange(0,W,nr),LinRange(0,H,nz),
#             LinRange(0-dr/2,W+dz/2,nr+1),LinRange(0-dz/2,H+dz/2,nz+1),
#             nr,nz,W,H)#,W/(nx-1),H/(ny-1))
#     end
# end

struct ScalarField
    #
    # This is a data structure to hold data defined on a (staggered) grid.
    # stagx and stagy are meant to indicate the offset (-1 or 0) in the x and y directions
    # Name is a field name, to be used in visualization
    # data contains the actual data.
    # old_data contains the previous data (last timestep), which is used to infill NaN values if there are no local markers.
    #
    stagx::Int64
    stagy::Int64
    name::String
    data::Matrix{Float64}
    old_data::Matrix{Float64}
    function ScalarField(stagx::Int64,stagy::Int64,name::String,data::Matrix{Float64})
        new(stagx,stagy,name,data)
    end
end

struct VectorField
    #
    # This is a data structure to hold data defined on a (staggered) grid.
    # stagx and stagy are meant to indicate the offset (-1 or 0) in the x and y directions
    # Name is a field name, to be used in visualization
    # data contains the actual data.
    #
    stagx::Int64
    stagy::Int64
    components::Int64
    name::String
    data::Array{Float64,3}
    function GridField(stagx::Int64,stagy::Int64,components,name::String,data)
        new(stagx,stagy,components,name,data)
    end
end

struct EquatorialGrid
    #
    # This data structure represents an equatorial slice
    # r is the radial coordinate
    # psi is the azimuthal coordinate
    # nr is the number of nodes in the radial direction
    # npsi is the number of nodes in the azimuthal direction
    # periodic indicates whether the grid is periodic (i.e. a complete spherical slice or a chunk)
    #
    r::Vector{Float64}
    psi::Vector{Float64}
    rc::Vector{Float64}
    psic::Vector{Float64}
    nr::Int64
    npsi::Int64
    periodic::Bool
    function EquatorialGrid(psi_max::Float64,npsi::Int64,r_min::Float64,r_max::Float64,nr::Int64; periodic::Bool=false)
        #
        # Construct a new equatorial grid with psi_min <= psi <= psi_max and r_min <= r <= r_max
        # if psi_max is 'close' to 2*pi, this is a complete spherical slice and periodic is set to true.
        #
        if psi_max ≈ 2.0*pi || periodic
            periodic = true
        else
            periodic = false
        end
        psi_min = 0.0
        dpsi = (psi_max-psi_min)/(npsi-1)
        dr   = (r_max-r_min)/(nr-1)
        psi  = LinRange(psi_min,psi_max,npsi)
        r    = LinRange(r_min,r_max,nr)
        rc   = LinRange(r_min-dr/2.,r_max+dr/2.,nr+1)
        psic = LinRange(psi_min-dpsi/2.,psi_max+dpsi/2.,npsi+1)
        new(r,psi,rc,psic,nr,npsi,periodic)
    end
end