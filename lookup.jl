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
        prop = Dict()
        prop["1"] = "rho,kg/m3"
        prop["2"] = "alpha,1/K"
        prop["3"] = "cp,J/K/kg"
        # prop["4"] = "vp,km/s"
        # prop["5"] = "vs,km/s"
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
        i=1
        for item in prop
            values[:,:,i] = loop_over_P ?
                reshape(dataframe[!, item.second], n_pressures, n_temperatures)' :
                reshape(dataframe[!, item.second], n_temperatures, n_pressures)
            i += 1
        end
        property_lookup = linear_interpolation((temperatures,pressures),values[:,:,index],extrapolation_bc=Line())
        return new(values, temperatures, pressures, property_lookup )
    end
end
