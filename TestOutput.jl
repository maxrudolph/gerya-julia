# if length(ARGS) > 9
#     error("specify proper input arguments for range function for ice shell thickness and wavelength topogaphy for ocean-ice interface")
# else
#     ice_start = parse(Float64,ARGS[1])
#     ice_stop = parse(Float64,ARGS[2])
#     ice_length = parse(Int64,ARGS[3])
#     wavelength_start = parse(Float64,ARGS[4])
#     wavelength_stop = parse(Float64,ARGS[5])
#     wavelength_length = parse(Int64,ARGS[6])
#     percent_amplitude_start = parse(Float64,ARGS[7])
#     percent_amplitude_stop = parse(Float64,ARGS[8])
#     percent_amplitude_length = parse(Int64,ARGS[9])
#     println("Ice Shell Thickness Range(start:difference between two subsequent objects:stop): ",range(ice_start,ice_stop,ice_length))
#     println("Wavelength Topography for Ocean-Ice Interface Range(start:difference between two subsequent objects:stop): ",range(wavelength_start,wavelength_stop,wavelength_length))
#     println("Amplitude Percentage Range(start:difference between two subsequent objects:stop): ",range(percent_amplitude_start,percent_amplitude_stop,percent_amplitude_length))
# end

if length(ARGS) > 4
    error("specify proper input arguments for range function for ice shell thickness and wavelength topogaphy for ocean-ice interface")
else
    ice_shell_thickness = parse(Float64,ARGS[1])
    wavelength = parse(Float64,ARGS[2])
    percent_amplitude = parse(Float64,ARGS[3])
    top_dir = ARGS[4]
    println("Model run for ice shell thickness of $ice_shell_thickness , wavelength of $wavelength , amplitude percentage of $percent_amplitude")
end

options = Dict()

# include("Outputs.jl")
include("Outputs_sbatch.jl")

function modelrun()
    sub_plots,sub_data = mk_sub_dir(top_dir)
    options["visualization file path"] = sub_data
    amp_decimal = percent_amplitude/100
    options["wavelength"] = wavelength*1e3
    options["ice thickness"] = ice_shell_thickness*1e3
    options["amplitude"] = amp_decimal*options["ice thickness"]
    options["surface depth"] = options["amplitude"]
    io = open(top_dir*"/output.txt","w")
    println(io,options)
    close(io)
end

# function modelrun()
#     top_dir = mk_modelrun_dir()
#     prompt_and_get_model_description(top_dir)
#     nlambda = wavelength_length
#     nhice = ice_length
#     namp = percent_amplitude_length
#     local lambda = range(wavelength_start,wavelength_stop,nlambda)
#     local hice =  range(ice_start,ice_stop,nhice)
#     local per_amp = range(percent_amplitude_start,percent_amplitude_stop,namp)

#     # Defining cases
#     max_ice_shell_thickness = maximum(hice)
#     min_ice_shell_thickness = minimum(hice)
#     max_wavelength = maximum(lambda)
#     min_wavelength = minimum(lambda)
#     cases = [
#         (max_ice_shell_thickness,max_wavelength), # Case 1
#         (min_ice_shell_thickness,min_wavelength), # Case 2
#         (min_ice_shell_thickness,max_wavelength), # Case 3 ice shell thickness < wavelength
#         (max_ice_shell_thickness,min_wavelength), # Case 4 ice shell thickness > wavelength
#     ]

#     irun = 1
#     for k in 1:namp
#         t_halfspace = zeros(nlambda,nhice)
#         t_rel = zeros(nlambda,nhice)
#         t_tic = zeros(nlambda,nhice)
#         amplitude_percentage = per_amp[k]
#         amp_decimal = amplitude_percentage/100
#         sub_dir = mk_sub_dir(top_dir,amplitude_percentage)
#         # for i in 1:nlambda
#         #     for j in 1:nhice
#                 for (ice_shell_thickness,wavelength) in cases
#                     options["wavelength"] = wavelength*1e3
#                     options["ice thickness"] = ice_shell_thickness*1e3
#                     options["amplitude"] = amp_decimal*options["ice thickness"]
#                     options["surface depth"] = options["amplitude"]
#                     # if options["wavelength"] >= options["ice thickness"]
#                         println("Starting model execution for model run $irun...")
#                         sub_dir_by_run,sub_dir_plots,sub_dir_data = mk_output_dir(sub_dir,irun)
#                         data_table_info(ice_start,ice_stop,nhice,wavelength_start,wavelength_stop,nlambda,sub_dir,amplitude_percentage)
#                         println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)"," , ","Using Amplitde Percentage: $amplitude_percentage%")
#                         io = open(sub_dir_by_run*"/output.txt","w")
#                         model_runtime(sub_dir_by_run,"Start")
#                         model_runtime(sub_dir_by_run,"End")
#                         close(io)
#                         println("Model ran successfully for model run $irun. Outputs saved to output.txt")
#                         irun += 1
#                     # end
#                 end
#         #     end
#         # end
#         if namp != 1
#             irun = 1
#         end
#     end
# end

try
    modelrun();
catch e
    println("Model encountered an error. Error details saved to error_log.txt")
    open("error_log.txt","w")do out
        redirect_stdout(out) do
            showerror(stdout,e,catch_backtrace())
        end
    end
    @error "Something went wrong" exception=(e,catch_backtrace())
end
