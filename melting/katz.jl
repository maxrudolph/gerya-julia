# Routines related to the Yasuda et al. 1994 model for melting of eclogite
# This was used by Leitch and Davies (2001) to calculate melt fraction in an 
# ascending eclogitic plume head
# the P/T curves were digitized using a web work-alike of DataTheif from Figure 6 of Yasuda et al. (1994)

using Interpolations
using MAT

struct katz
    melt_frac
    temperatures::Array{Float64,1}
    pressures::Array{Float64,1}
    melt_lookup
    function katz()    
        filename = matopen("melting/melt_table_0.mat")
        meta = read(filename,"melt_table")
        pressure_melt = unique(meta["P"])
        temperature_melt = unique(meta["T"])
        melt_frac = meta["F"]
        
        melt_lookup = linear_interpolation((pressure_melt,temperature_melt),melt_frac,extrapolation_bc=Line())
        return new(melt_frac, temperature_melt, pressure_melt, melt_lookup )
    end
end