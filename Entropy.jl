using Interpolations
using HDF5

fname = "SeaFreeze_table_with_density.hdf5"
# Opening file
fid = h5open(fname, "r");
# Reading Data from file
entropy = fid["Entropy"];
pressure = fid["Pressure"];
temperature = fid["Temperature"];
density = fid["Density"];
# Accessing Entropy contents
ed = read(entropy, "Data");
# Accessing Pressure content that I want
pd = read(pressure, "Data");
# Accessing Temperature content that I want
td = read(temperature, "Data");
# Accessing Density content that I want
rhod = read(density, "Data");
# Close file
close(fid)

entropy = ed;
pressure = pd*1e6;
temperature = td;
density = rhod;
# linear interpolation for entropy
interp_linear_entropy = linear_interpolation((pressure,temperature),entropy)

# linear interpolation for density
interp_linear_density = linear_interpolation((pressure,temperature),density)

function update_entropy(markers::Markers)
    for i in 1:markers.nmark
        markers.scalars[markers.scalarFields["En"],i] = interp_linear_entropy(markers.scalars[markers.scalarFields["hp"],i],markers.scalars[markers.scalarFields["T"],i])
    end
end

function analytic_hydrostatic_pressure(grid::CartesianGrid,density::Array,gy::Float64)
    """
    Arguments:
    depth -- the depth in (m)
    density -- the density of the fluid in (kg/^3)
    gravity -- the acceleration due to gravity in (m/s^2)
    
    Returns: 
    pressure -- the hydostatic pressure in (Pa)
    """
    
    pressure = zeros(Float64,grid.ny+1)
    hice = 1.5e4
    hice = initial_ice_depth(1.0)
    hsurf = initial_surface_depth(1.0)

    for i in 1:grid.ny+1
        if grid.yc[i] < hsurf
            pressure[i] = density[1]*gy*grid.yc[i]
        elseif hsurf < grid.yc[i] < hice
            Pair = density[1]*gy*hsurf
            pressure[i] = Pair + density[2]*gy*(grid.yc[i]-hsurf)
        else
            Pair = density[1]*gy*hsurf
            Pice = density[2]*gy*(hice-hsurf)
            pressure[i] = Pair + Pice + density[3]*gy*(grid.yc[i]-hice)
        end
    end
    return pressure
end

function numerical_hydrostatic_pressure(grid::CartesianGrid,rho_vy::Matrix{Float64},gy::Float64)
    """
    Arguments:
    rho_vy -- the density on the y-velocity node in(kg/m^3)
    gy -- the acceleration due to gravity in the y-direction in (m/s^2)
    
    Returns: 
    hydro_press -- the hydostatic pressure in (Pa)
    """
    
    hydro_press = zeros(Float64,(grid.ny+1,grid.nx+1)) 
    for j in 1:grid.nx+1
        for i in 1:grid.ny+1
            if i == 1
                hydro_press[i,j] = gy*rho_vy[1,j]*grid.yc[i]
            else
                hydro_press[i,j] = hydro_press[i-1,j] + gy*rho_vy[i-1,j]*(grid.yc[i]-grid.yc[i-1])
            end
        end 
    end 
    return hydro_press
end

