struct CartesianGrid
    x::Array{Float64,1}
    y::Array{Float64,1}
    xc::Array{Float64,1}
    yc::Array{Float64,1}
    nx::Int64
    ny::Int64
    W::Float64
    H::Float64
    #dx::Float64
    #dy::Float64
    function CartesianGrid(W::Float64,H::Float64,nx::Int,ny::Int)
        dx = W/(nx-1)
        dy = H/(ny-1)
        new(LinRange(0,W,nx),LinRange(0,H,ny),
            LinRange(0-dx/2,W+dx/2,nx+1),LinRange(0-dy/2,H+dy/2,ny+1),
            nx,ny,W,H)#,W/(nx-1),H/(ny-1))
    end
end

struct EquatorialGrid
    r::Vector
    psi::Vector
    rc::Vector
    psic::Vector
    nr::Int
    npsi::Int
    dr::Float64
    dpsi::Float64
    periodic::Bool
    function EquatorialGrid(psi_max,npsi,r_min,r_max,nr; periodic=false)
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
        new(r,psi,rc,psic,nr,npsi,dr,dpsi,periodic)
    end
end