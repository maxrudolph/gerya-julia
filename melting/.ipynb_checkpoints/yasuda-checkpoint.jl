# Routines related to the Yasuda et al. 1994 model for melting of eclogite
# This was used by Leitch and Davies (2001) to calculate melt fraction in an 
# ascending eclogitic plume head
# the P/T curves were digitized using a web work-alike of DataTheif from Figure 6 of Yasuda et al. (1994)

using Interpolations
using CSV
using DataFrames

struct yasuda
    solidus_P::Array{Float64,1}
    solidus_T::Array{Float64,1}
    liquidus_P::Array{Float64,1}
    liquidus_T::Array{Float64,1}
    solidus
    liquidus
    pressure
    function yasuda()
        solidus = CSV.read("melting/yasuda_1994_solidus.csv",DataFrame);
        liquidus = CSV.read("melting/yasuda_1994_liquidus.csv",DataFrame);
        
        # Pressure-to-depth conversion from Yasuda et al. 1994 Figure 6:
        P = [3.00751879699248,6.466165413533833,9.962406015037592,13.308270676691729,17.105263157894733,21.428571428571423]
        depth = [1e5,2e5,3e5,4e5,5e5,6e5]
        
        return new(solidus[:,1],solidus[:,2],liquidus[:,1],liquidus[:,2],linear_interpolation(solidus[:,1],solidus[:,2],extrapolation_bc=Line()), linear_interpolation(liquidus[:,1],liquidus[:,2],extrapolation_bc=Line()),linear_interpolation(depth,P,extrapolation_bc=Line()) )
    end
end



function melt_fraction(model::yasuda,depth::Float64,T::Float64)
    P = model.pressure(depth)
    Tliq = model.liquidus(P)
    Tsol = model.solidus(P)
    if T<Tsol
        return 0.0
    elseif T>Tliq
        return 1.0
    else
        # assume linear variation of melt fraction between solidus and liquidus per Leitch and Davies
        theta = (T-Tsol)/(Tliq-Tsol)
        return theta        
    end
end
