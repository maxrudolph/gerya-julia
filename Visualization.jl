
function visualization(grid::CartesianGrid,rho::Matrix,eta::Matrix,vn::Array{Float64},pressure::Matrix,temperature::Matrix,time ; filename="test.vts")
    # write the visualization output from the regular grid as a .vts file.
    vtk_grid(filename, grid.x, grid.y) do vtk
        vtk["rho"] = transpose(rho)
        vtk["viscosity"] = transpose(eta)
        # add a fake third dimension to the velocity vectors
        v3 = Array{Float64,3}(undef,3,grid.nx,grid.ny)
        v3[1:2,:,:] = vn
        v3[3,:,:] .= 0.0
        vtk["Velocity"] = v3
        vtk["Temperature"] = transpose(temperature)
        vtk["pressure"] = transpose(pressure[2:end,2:end])
        vtk["TimeValue"] = time
    end
end

function visualization(grid::CartesianGrid,output_fields::Dict,time ; filename="test.vts")
    # write the visualization output from the regular grid as a .vts file.
    vtk_grid(filename, grid.x, grid.y) do vtk
        for key in keys(output_fields)
            if length(size(output_fields[key])) == 3
                # this is a 3D array - assume that it's (vector) velocity.
                # add a fake third dimension to the velocity vectors
                v3 = Array{Float64,3}(undef,3,grid.nx,grid.ny)
                v3[1:2,:,:] = output_fields[key]
                v3[3,:,:] .= 0.0
                vtk[key] = v3
            else
                vtk[key] = transpose(output_fields[key])
            end
        end
        vtk["TimeValue",VTKFieldData()] = time
    end
end

function visualization(markers::Markers,time; filename="markers.vtp")  
    p3 = Array{Float64,2}(undef,3,markers.nmark)
    p3[1:2,:] = markers.x[1:2,1:markers.nmark]
    p3[3,:] .= 0.0
      
    @time polys = [MeshCell(PolyData.Polys(),i:i) for i in 1:markers.nmark]
    vtk_grid(filename,p3,polys) do vtk    
    #vtk_grid(filename,p3) do vtk
        for key in keys(markers.scalarFields)
            vtk[key] = markers.scalars[markers.scalarFields[key],1:markers.nmark]
        end
        for key in keys(markers.integerFields)
            vtk[key] = Float64.(markers.integers[markers.integerFields[key],1:markers.nmark])
        end
       vtk["TimeValue"] = time
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

function plots(markers::Markers,grid::CartesianGrid,rho::Matrix,time::Float64 ; draw_grid::Bool = false)
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