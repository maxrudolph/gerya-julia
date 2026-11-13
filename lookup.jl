# Material thermal properties from lookup table using 2d interpolation func

using Interpolations
using CSV
using DataFrames

struct lookup
    values
    temperatures::Array{Float64,1}
    pressures::Array{Float64,1}
    property_lookup
    
    function lookup(filename::String,index::Int64)
        # column names in table order - index selects one of these.
        # (a Dict must not be used here because its iteration order is not guaranteed)
        prop = ["rho,kg/m3", "alpha,1/K", "cp,J/K/kg"]
        # "vp,km/s", "vs,km/s"
        # Read the header and data efficiently
        dataframe = CSV.File(filename; delim=' ', skipto=14, silencewarnings=true,ntasks=1) |> DataFrame
    
        # Assume the header is in the first skipped line
        header = CSV.File(filename; delim=' ', skipto=13, limit=1, silencewarnings=true,ntasks=1) |> DataFrame |> first |> collect
        rename!(dataframe, header)
    
        # Determine if loop_over_P
        loop_over_P = dataframe[1, "T(K)"] == dataframe[2, "T(K)"]
    
        # Extract unique values and compute deltas
        temperatures = unique(dataframe[!, "T(K)"])
        pressures = unique(dataframe[!, "P(bar)"]) .* 1e5  # Convert bar to Pa
        
    
        n_temperatures = length(temperatures)
        n_pressures = length(pressures)
        values = zeros(n_temperatures,n_pressures,length(prop))
        
        # Reshape densities once
        for (i,column) in enumerate(prop)
            values[:,:,i] = loop_over_P ?
                reshape(dataframe[!, column], n_pressures, n_temperatures)' :
                reshape(dataframe[!, column], n_temperatures, n_pressures)
        end
        property_lookup = linear_interpolation((temperatures,pressures),values[:,:,index],extrapolation_bc=Line())
        return new(values, temperatures, pressures, property_lookup )
    end
end
