function visualization(grid::CartesianGrid,rho::Matrix,eta::Matrix,vn::Array{Float64},time ; filename="test.vts")
    # write the visualization output from the regular grid as a .vts file.
    vtk_grid(filename, grid.x, grid.y) do vtk
        vtk["rho"] = transpose(rho)
        vtk["viscosity"] = transpose(eta)
        # add a fake third dimension to the velocity vectors
        v3 = Array{Float64,3}(undef,3,grid.nx,grid.ny)
        v3[1:2,:,:] = vn
        v3[3,:,:] .= 0.0
        vtk["Velocity"] = v3
        vtk["TIME"] = time
    end
end

function draw_grid(grid::CartesianGrid ; nodes::Bool=false)
   for y in grid.y
        plot([grid.x[1],grid.x[end]],[y,y],"k")
    end
    for x in grid.x
        plot([x,x],[grid.y[1],grid.y[end]],"k") 
    end
    if nodes
        for i in 1:grid.ny
            for j in 1:grid.nx
                plot(grid.y[i],grid.xc[j],"gs")
                plot(grid.yc[i],grid.x[j],"bo")
                plot(grid.yc[i],grid.xc[j],"ro")
            end
        end
    end
end

function plots(markers::Markers,grid::CartesianGrid,rho::Matrix,time::Float ; draw_grid::Bool = false)
    @views mx = markers.x[:,1:markers.nmark]
    @views mrho = markers.rho[1:markers.nmark]
    figure()
    subplot(1,2,1)
    #pcolor(grid.x,grid.y,rho)
    scatter(mx[1,:],mx[2,:],c=mrho,s=0.1)
    colorbar()
    gca().invert_yaxis()
    title(time)
    gca().set_aspect("equal")
    if draw_grid
        draw_grid(grid)
    end
    subplot(1,2,2)
    contourf(grid.x,grid.y,rho,)
    colorbar()
    gca().invert_yaxis()
    title(time)
    gca().set_aspect("equal")
end