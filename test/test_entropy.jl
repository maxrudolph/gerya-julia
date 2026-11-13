# Unit tests for the entropy solver (TemperatureEntropy.jl): boundary entropy at a fixed-temperature
# boundary, subgrid entropy diffusion, and the Stefan-problem similarity variable get_theta.
# Run from the repository root:  julia --project test/test_entropy.jl
# (needs Test, Random, NearestNeighbors, SpecialFunctions and PyPlot)
using Test, Random, NearestNeighbors, SpecialFunctions
using PyPlot   # TemperatureEntropy.jl uses PyPlot's L"..." strings in its plotting routines
include(joinpath(@__DIR__,"..","Grid.jl"))
include(joinpath(@__DIR__,"..","Markers.jl"))
include(joinpath(@__DIR__,"..","TemperatureEntropy.jl"))

opts = Dict("latent heat of fusion"=>334e3, "Tm"=>273.0, "specific heat of ice"=>2096.0)
Hfus, Tm, Cv = opts["latent heat of fusion"], opts["Tm"], opts["specific heat of ice"]
S_of(T,X) = compute_S_from_T_X(X,T,opts)

# Subgrid diffusion relaxes a marker's entropy toward the value interpolated from the nodes.
# Reference: one marker exchanging heat with a node held at Tn, rho*T dS/dt = k*(2/dx^2+2/dy^2)*(Tn - T),
# integrated with RK4. Heat flow is driven by the temperature difference only.
const rho, kThermal = 920.0, 2.1
function reference_entropy(S0,Sn,dt,c)
    Tn, = compute_T_X_from_S(Sn,opts)
    f(s) = kThermal*c*(Tn - compute_T_X_from_S(s,opts)[1])/(rho*compute_T_X_from_S(s,opts)[1])
    s = S0; n = 4000; h = dt/n
    for _ in 1:n
        k1 = f(s); k2 = f(s+h/2*k1); k3 = f(s+h/2*k2); k4 = f(s+h*k3)
        s += h/6*(k1+2k2+2k3+k4)
    end
    return s
end

# One subgrid step for markers that all start at entropy S0, with every node at entropy Sn.
function subgrid_step(S0,Sn,dt)
    grid = CartesianGrid(1.0,1.0,3,3)
    m = Markers(grid,["rho","T","S","kThermal","X"],["mat"];nmx=1,nmy=1)
    T0,X0 = compute_T_X_from_S(S0,opts)
    for (field,value) in (("rho",rho),("kThermal",kThermal),("S",S0),("T",T0),("X",X0))
        m.scalars[m.scalarFields[field],1:m.nmark] .= value
    end
    subgridSdiff(grid,m,fill(Sn,grid.ny+1,grid.nx+1),dt,opts)
    return m.scalars[m.scalarFields["S"],1:m.nmark]
end

@testset "TemperatureEntropy.jl" begin

@testset "boundary entropy at T = Tm is set per column" begin
    grid = CartesianGrid(8.0,4.0,9,5)
    X = zeros(grid.ny+1,grid.nx+1)
    X[grid.ny,:] = [0.0, 0.0, 1.0, 1.0, 0.5, 1.0, 0.0, 1.0, 1.0, 0.0]   # ice, water and mush along the bottom row
    T = fill(Tm,grid.ny+1,grid.nx+1)
    S = [S_of(T[i,j],X[i,j]) for i in 1:grid.ny+1, j in 1:grid.nx+1]
    Tp,Xp,Sp = ghost_nodes_center_TXS(grid,T,X,S,[-1,-1,1,1],[0.0,0.0,100.0,Tm],opts)
    cols = 2:grid.nx
    # entropy on the boundary (average of ghost and interior) equals that of the adjacent cell in each column
    @test (Sp[grid.ny+1,cols] .+ Sp[grid.ny,cols]) ./ 2 ≈ Hfus/Tm .* X[grid.ny,cols]
    @test Xp[grid.ny+1,cols] == X[grid.ny,cols]   # zero-gradient melt fraction, stays within [0,1]
end

@testset "subgrid entropy diffusion" begin
    c = 2/0.5^2 + 2/0.5^2                # 2/dx^2 + 2/dy^2 for the 3x3 grid in subgrid_step
    tdiff = rho*Cv/(kThermal*c)          # diffusion time of a solid marker
    @testset "mushy marker, mushy node: same temperature, no change" begin
        S0, Sn = S_of(Tm,0.3), S_of(Tm,0.6)
        for r in (0.1, 1.0, 100.0)
            @test all(subgrid_step(S0,Sn,r*tdiff) .== S0)
        end
    end
    @testset "mushy marker, colder solid node: the marker freezes" begin
        S0, Sn = S_of(Tm,0.3), S_of(250.0,0.0)
        # for short steps the entropy change follows the reference (heat flow set by Tm - 250 K)
        for r in (0.01, 0.1)
            dS = subgrid_step(S0,Sn,r*tdiff)[1] - S0
            dS_ref = reference_entropy(S0,Sn,r*tdiff,c) - S0
            @test dS < 0
            @test isapprox(dS, dS_ref, rtol=0.03)
        end
        # the change grows with the step and approaches the nodal entropy
        steps = [subgrid_step(S0,Sn,r*tdiff)[1] for r in (0.01,0.1,1.0,10.0,100.0)]
        @test issorted(steps, rev=true)
        @test steps[end] ≈ Sn rtol=1e-3
    end
    @testset "solid marker, colder solid node: relaxes with Cv" begin
        S0, Sn = S_of(250.0,0.0), S_of(200.0,0.0)
        dt = tdiff
        @test subgrid_step(S0,Sn,dt) ≈ fill(S0 + (Sn-S0)*(1-exp(-dt/tdiff)),4)
    end
end

@testset "get_theta takes kappa as an argument" begin
    kappa = 1e-6; lambda1 = 0.4; t = 1e10
    @test !isdefined(Main,:options)          # no global options is needed
    @test get_theta(0.0,t,lambda1,kappa) == 0.0
    y = get_y(lambda1,t,Dict("thermal diffusivity"=>kappa))   # position of the melting front
    @test get_theta(y,t,lambda1,kappa) ≈ 1.0
end

end # TemperatureEntropy.jl
