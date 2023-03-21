# Usage: julia PlotMelting.jl [output_directory output_directory2...] start_time end_time
# Give start time and end time in millions of years
# output will be written in output_directory.
using CSV, DataFrames, FastRunningMedian
using PyPlot

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
# Create Figure Windows
fig1,ax1 = plt.subplots(1,1, figsize=(4,2.5) )
fig2,ax2 = plt.subplots(1,1, figsize=(4,2.5) )
fig4,ax4 = plt.subplots(1,1, figsize=(4,2.5) )
for i in 1:nfiles
    # Read statistics file:

    filename = output_dirs[i] * "/statistics.txt"
    println("Reading from: ",filename)
    melt_output = CSV.read(filename,DataFrame,comment="#");
    
    time = melt_output[:,2]/3.15e7/1e6
    melt_km3yr = melt_output[:,3] .* (3.15e7/1e9)*eclogite_fraction
    ind = findfirst( melt_km3yr .> 0.0 .&& time .> 2.0 )
    
    mask = (time .> start_time) .&& (time .<= end_time) 
    mask[1:ind] .= false

    # melt output is in units of m^3/s
    time=time[mask]
    time = time .- time[1]
    melt_km3yr = melt_output[mask,3] .* (3.15e7/1e9)*eclogite_fraction
    melt_km3yr = running_median(melt_km3yr, 1)

    ax1.plot(time,melt_km3yr,label=output_dirs[i])
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
    ax3.plot(time,cumsum)
    ax3.set_xlim([0.0,5.0])
    ax3.set_xlabel("Time (Myr)")
    ax3.set_ylabel("Cumulative Melt (10\$^6\$ km\$^3\$)")
    fig3.savefig(output_dirs[i] * "/cumulative_melt.png")

    ax4.plot(time,cumsum,label=output_dirs[i])
    ax4.set_xlim([0.0,5.0])

end
ax1.set_ylim([1e-2,5e1])
ax1.set_xlim([0.0,5.0])
ax1.legend()
fig1.savefig("./melt_vs_time.eps",bbox_inches="tight")

ax4.legend()
ax4.set_xlabel("Time (Myr)")
ax4.set_ylabel("Cumulative Melt (10\$^6\$ km\$^3\$)")
fig4.savefig("cumulative_melt.png")

