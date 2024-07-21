function fitting_data(topography::Vector{Any},times::Vector{Any},itime::Int64,sub_plots::String)
    ampfit = []
    for i in 1:itime-1
        amp = maximum(topography[i])-minimum(topography[i])
        append!(ampfit,amp)
    end
    x = convert(Array{Float64},times)
    y = convert(Array{Float64},ampfit)
    fit = fitexp(x,y)
    fitted_time = fit.b

    figure()
    for i in 1:itime-1
        plot(times[i]/3.15e7,maximum(topography[i])-minimum(topography[i]),".")
    end
    plot(fit.x,fit.y,"r-",label="fitted data")
    gca().set_ylabel(L"Amplitude\,(m)")
    gca().set_xlabel(L"Time\,(yrs)")
    #legend(fancybox="True",shadow="True")
    savefig(sub_plots*"/fitted_data.png",dpi=300)
    close()
    return fitted_time/3.15e7
end

function hdf5_file(options::Dict,t_hs::Float64,t_rel::Float64,t_thic::Float64,t_fit::Float64,top_dir::String)
    hice = options["hice"]
    lambda = options["wavelength"]
    lambda = vcat(map(x -> x', lambda)...)
    hice = vcat(map(x -> x', hice)...)
    h5open(top_dir*"/data.hdf5","w") do file
        println("Creating HDF5 File")
        println("Saving Data into a HDF5 File")
        # Creating Groups for Data
        g = create_group(file, "Model Run")
        # Storing Data Inside the Group
        g["Wavelength"] = lambda[:]
        g["Ice Shell Thickness"] = hice[:]
        g["Viscous Relaxation Time(Half-Space)"] = t_hs
        g["Viscous Relaxation Time(Model)"] = t_rel
        g["Fitted Viscous Relaxation Time"] = t_fit
        g["Thickening Time"] = t_thic
        # Apply an Attribute to Groups
        attrs(g)["Description"] = "This group contains only a 4 dataset"
        println("Finished Saving Data into a HDF5 File")
    end
end
