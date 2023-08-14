using EasyFit
using PyPlot

function fitting_data(topography::Vector{Any},times::Vector{Any},max_step::Int64)
    ampfit = []
    for i in 1:max_step
        amp = maximum(topography[i])-minimum(topography[i])
        append!(ampfit,amp)
    end
    x = convert(Array{Float64}, times)
    y = convert(Array{Float64}, ampfit)
    fit = fitexp(x/3.15e7,y)
    figure()
    for i in 1:max_step
        plot(times[i]/3.15e7,maximum(topography[i])-minimum(topography[i]),".") 
    end
    plot(fit.x,fit.y,"r-",label="fitted data")
    gca().set_ylabel(L"Amplitude\,(m)")
    gca().set_xlabel(L"Time\,(years)")
    legend(fancybox="True",shadow="True")
    show()
end