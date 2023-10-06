function fitting_data(topography::Vector{Any},times::Vector{Any},itime::Int64,sub_dir_plots::String)
    ampfit = []
    for i in 1:itime-1
        amp = maximum(topography[i])-minimum(topography[i])
        append!(ampfit,amp)
    end
    x = convert(Array{Float64}, times)
    y = convert(Array{Float64}, ampfit)
    fit = fitexp(x,y)
        # ,options=Options(besttol=100))
    # # figure()
    # for i in 1:itime-1
    #     plot(times[i]/3.15e7,maximum(topography[i])-minimum(topography[i]),".") 
    # end
    # plot(fit.x,fit.y,"r-",label="fitted data")
    # gca().set_ylabel(L"Amplitude\,(m)")
    # gca().set_xlabel(L"Time\,(years)")
    # legend(fancybox="True",shadow="True")
    # savefig(sub_dir_plots*"/Fitted_Data.pdf")
    # close()
    fitted_time = fit.b
    return fitted_time/3.15e7
end