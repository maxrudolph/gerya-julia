using CSV, DataFrames, FastRunningMedian
using PyPlot

output_directory = "plume_test_300_high"

filename = output_directory *"/statistics.txt"
println(filename)
melt_output = CSV.read(filename,DataFrame,comment="#");
time = melt_output[:,2]/3.15e7/1e6

mask = time .> 20.0

# melt output is in units of m^3/s
time=time[mask]
melt_km3yr = melt_output[mask,3] .* (3.15e7/1e9)
melt_km3yr = running_median(melt_km3yr, 1)

# determine an appropriate plot range
#mask = melt_km3yr .>0 
figure(figsize=(4,2.5))
plot(time.-282.0,melt_km3yr)
gca().set_yscale("log")
gca().set_xlim([0,20])
gca().set_ylim([1e-2,1e1])
gca().set_xlabel("Time (Myr)")
gca().set_ylabel("Melt Production (-)")
#show()
savefig(output_directory * "/melt_vs_time.eps",bbox_inches="tight")

n = length(melt_km3yr)
cumsum = zeros(n,1)

for i=2:n
    cumsum[i] = cumsum[i-1]+melt_km3yr[i]*(time[i]-time[i-1])
end
figure()
plot(time,cumsum)
#gca().set_xlim([120,130])
savefig(output_directory * "/cumulative_melt.png")
