using CSV, DataFrames, FastRunningMedian
using PyPlot

output_directory = "plume_test_550"

filename = output_directory *"/statistics.txt"
println(filename)
melt_output = CSV.read(filename,DataFrame,comment="#");
time_out = melt_output[:,2]/3.15e7/1e6
mask = time_out .>= 20

# melt output is in units of m^3/s

melt_km3yr = melt_output[mask,3] .* (3.15e7/1e9)
melt_km3yr = running_median(melt_km3yr, 1)
figure()
plot(time_out[mask],melt_km3yr)
gca().set_yscale("log")
#show()
savefig(output_directory * "/melt_vs_time.png")

n = length(melt_km3yr)
cumsum = zeros(n,1)

for i=2:n
    cumsum[i] = cumsum[i-1]+melt_km3yr[i]*(time_out[i]-time_out[i-1])
end
figure()
plot(time_out[mask],cumsum)
#gca().set_xlim([120,130])
savefig(output_directory * "/cumulative_melt.png")
