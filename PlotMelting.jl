# Usage: julia PlotMelting.jl [output_directory output_directory2...] start_time end_time
# Give start time and end time in millions of years
# output will be written in output_directory.
using CSV, DataFrames, FastRunningMedian
using PyPlot
using Interpolations

eclogite_fraction = 0.15
# Parse arguments:
output_dirs = String[]

if length(ARGS) > 1
   nfiles = length(ARGS) - 2
   for i in 1:nfiles
       push!(output_dirs,ARGS[i])
   end
   start_time = parse( Float64, ARGS[nfiles+1] )
   end_time = parse( Float64, ARGS[nfiles+2] )
end

# parse temperatures and lithospheric thicknesses
line_colors = []
line_styles = []
for i in 1:nfiles
    parts = split(output_dirs[i],"_") # split on underscore
    temp = parse(Float64,parts[2])
    thickness = parse(Float64,parts[3])
    if thickness > 75000.0
       push!( line_styles, "--")
    else
       push!( line_styles, "-")
    end
    if temp < 350.0
       push!( line_colors, [241,193,107]/255 )
    elseif temp < 450.0
        push!( line_colors, [194,111,42]/255 )
    else
	push!( line_colors, [255,33,12]/255 )
    end
end
println(line_styles)
println(line_colors)
#break


# Create Figure Windows
fig1,ax1 = plt.subplots(1,1, figsize=(4,2.5) )
fig2,ax2 = plt.subplots(1,1, figsize=(4,2.5) )
fig5,ax5 = plt.subplots(1,1, figsize=(4,2.5) )
fig4,ax4 = plt.subplots(1,1, figsize=(4,2.5) )
fig6,ax6 = plt.subplots(1,1, figsize=(4,2.5) )


for ifile in 1:nfiles
    # Read statistics file:

    filename = output_dirs[ifile] * "/statistics.txt"
    println("Reading from: ",filename)
    melt_output = CSV.read(filename,DataFrame,comment="#");
    
    time = melt_output[:,2]/3.15e7/1e6
    melt_km3yr = melt_output[:,3] .* (3.15e7/1e9)*eclogite_fraction
    
    ind = findfirst( melt_km3yr .> 0.0 .&& time .> 2.0 )
    
    mask = (time .> start_time) .&& (time .<= end_time) 
    mask[1:ind] .= false

    # melt output is in units of m^3/s
    time = time[mask]
    time = time .- time[1]
    melt_km3yr = melt_output[mask,3] .* (3.15e7/1e9)*eclogite_fraction
    melt_km3yr = running_median(melt_km3yr, 1)
    # carbon should be measured in units of kg/yr
    # (integrated wt ppm (ppm*m^3)) * (kg/m^3) * (1 tonne)/(1000 kg)
    carbon = melt_output[mask,4] .* (1e-6) .* 3300.0 /1e3

    ax1.plot(time,melt_km3yr,label=output_dirs[ifile],color=line_colors[ifile],linestyle=line_styles[ifile])
    ax1.set_yscale("log")
    ax1.set_xlabel("Time (Myr)")
    ax1.set_ylabel("Melt Production (km\$^3\$/yr)")

    # make a figure in each output folder with the cumulative amount of melting.
    n = length(melt_km3yr)
    cumsum = zeros(n,1)
    
    for i=2:n
        cumsum[i] = cumsum[i-1]+melt_km3yr[i]*(time[i]-time[i-1])
    end
    fig3,ax3 = plt.subplots(1,1)
    ax3.plot(time,cumsum,color=line_colors[ifile],linestyle=line_styles[ifile])
    ax3.set_xlim([0.0,5.0])
    ax3.set_xlabel("Time (Myr)")
    ax3.set_ylabel("Cumulative Melt (10\$^6\$ km\$^3\$)")
    fig3.savefig(output_dirs[ifile] * "/cumulative_melt.png")

    ax4.plot(time,cumsum,label=output_dirs[ifile],color=line_colors[ifile],linestyle=line_styles[ifile])
    ax4.set_xlim([0.0,5.0])

    # cumulative amount of carbon
    n = length(carbon)
    cumsum = zeros(n)
    for i=2:n
        cumsum[i] = cumsum[i-1]+carbon[i]
    end
    ax6.plot(time,cumsum,label=output_dirs[ifile],color=line_colors[ifile],linestyle=line_styles[ifile])
    ax6.set_xlim([0.0,5.0])

    # put the carbon release onto a uniformly-spaced time vector
    n1 = 200
    time_plot = LinRange(0.0,10.0,n1)
    interpolant = linear_interpolation(time,cumsum)
    carbon_plot = interpolant(time_plot)
    dCdt = zeros(n1)
    for i in 1:n1-1
    	dCdt[i+1] = (carbon_plot[i+1]-carbon_plot[i])/(time_plot[i+1]-time_plot[i])/1e6
    end
    
    ax5.plot(time_plot,dCdt,label=output_dirs[ifile],color=line_colors[ifile],linestyle=line_styles[ifile])
    ax5.set_yscale("log")
    ax5.set_xlabel("Time (Myr)")
    ax5.set_ylabel("Carbon release (Tonne/yr)")

end
ax1.set_ylim([1e-2,5e1])
ax1.set_xlim([0.0,5.0])
ax1.legend()
fig1.savefig("./melt_vs_time.eps",bbox_inches="tight")

ax4.legend()
ax4.set_xlabel("Time (Myr)")
ax4.set_ylabel("Cumulative Melt (10\$^6\$ km\$^3\$)")
fig4.savefig("cumulative_melt.eps",bbox_inches="tight")

ax5.legend()
ax5.set_xlabel("Time (Myr)")
ax5.set_ylabel("CO$_2$ (Tonne/yr)")
fig5.savefig("carbon_vs_time.eps",bbox_inches="tight")

ax6.legend()
ax6.set_xlim([0.0,10.0])
ax6.set_xlabel("Time (Myr)")
ax6.set_ylabel("CO$_2$ (Tonne)")
fig6.savefig("cumulative_carbon.eps",bbox_inches="tight")
