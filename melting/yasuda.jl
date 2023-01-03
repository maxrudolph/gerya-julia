# Routines related to the Yasuda et al. 1994 model for melting of eclogite
# This was used by Leitch and Davies (2001) to calculate melt fraction in an 
# ascending eclogitic plume head
# the P/T curves were digitized using a web work-alike of DataTheif from Figure 6 of Yasuda et al. (1994)

using Interpolations
using CSV

struct yasuda
    solidus_P::Array{Float64,1}
    solidus_T::Array{Float64,1}
    liquidus_P::Array{Float64,1}
    liquidus_T::Array{Float64,1}

    function yasuda()
        solidus = CSV.read("yasuda_1994_solidus.csv",DataFrame)
        liquidus = CSV.read("yasuda_1994_liquidus.csv",DataFrame)
        
        return new(solidus[:,1],solidus[:,2],liquidus[:,1],liquidus[:,2])
    end
end

