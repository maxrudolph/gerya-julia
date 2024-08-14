function fitting_amp_data(amplitude::Vector{Any},times::Vector{Any},itime::Int64,sub_plots::String)
    ampfit = []
    for i in 1:itime-1
        amp = amplitude[i]
        append!(ampfit,amp)
    end
    x = convert(Array{Float64},times)
    y = convert(Array{Float64},ampfit)
    fit = fitexp(x/3.15e7,y,options=Options(fine=100))
    fitted_time_amp = fit.b

    figure()
    for i in 1:itime-1
        plot(times[i]/3.15e7,amplitude[i],".")
    end
    plot(fit.x,fit.y,"r-",label="fitted data")
    title(L"Amplitude\,\,Over\,\,Time")
    gca().set_ylabel(L"Amplitude\,(m)")
    gca().set_xlabel(L"Time\,(yrs)")
    legend()
    savefig(sub_plots*"/ampfitted_data.png",dpi=300)
    close()
    return fitted_time_amp/3.15e7
end

function fitting_thickingd_data(thickness::Vector{Any},times::Vector{Any},itime::Int64,sub_plots::String)
    hicefit = []
    for i in 1:itime-1
        hice = thickness[i]
        append!(hicefit,hice)
    end
    x = convert(Array{Float64},times)
    y = convert(Array{Float64},hicefit)
    fit = fitexp(x/3.15e7/1e3,y,options=Options(fine=100))
    fitted_time_hice = fit.b

    figure()
    for i in 1:itime-1
        plot(times[i]/3.15e7/1e3,thickness[i],".")
    end
    plot(fit.x,fit.y,"r-",label="fitted data")
    title(L"Ice\,\,Shell\,\,Thickness\,\,Over\,\,Time")
    gca().set_ylabel(L"Ice Shell Thickness\,(m)")
    gca().set_xlabel(L"Time\,(kyr)")
    legend()
    savefig(sub_plots*"/hicefitted_data.png",dpi=300)
    close()
    return fitted_time_hice/3.15e7
end

function hdf5_file(options::Dict,t_hs::Float64,t_rel::Float64,t_thic::Float64,t_rel_fit::Float64,t_thick_fit::Float64,top_dir::String)
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
        g["Fitted Viscous Relaxation Time"] = t_rel_fit
        g["Thickening Time"] = t_thic
        g["Fitted Thickening Time"] = t_thick_fit
        # Apply an Attribute to Groups
        attrs(g)["Description"] = "This group contains only a 4 dataset"
        println("Finished Saving Data into a HDF5 File")
    end
end
