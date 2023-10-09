if length(ARGS) > 9
    error("specify proper input arguments for range function for ice shell thickness and wavelength topogaphy for ocean-ice interface")
else
    ice_start = parse(Float64,ARGS[1])
    ice_stop = parse(Float64,ARGS[2])
    ice_length = parse(Int64,ARGS[3])
    wavelength_start = parse(Float64,ARGS[4])
    wavelength_stop = parse(Float64,ARGS[5])
    wavelength_length = parse(Int64,ARGS[6])   
    percent_amplitude_start = parse(Float64,ARGS[7])
    percent_amplitude_stop = parse(Float64,ARGS[8])
    percent_amplitude_length = parse(Int64,ARGS[9])
    println("Ice Shell Thickness Range(start:difference between two subsequent objects:stop): ",range(ice_start,ice_stop,ice_length))
    println("Wavelength Topography for Ocean-Ice Interface Range(start:difference between two subsequent objects:stop): ",range(wavelength_start,wavelength_stop,wavelength_length))
    println("Amplitude Percentage Range(start:difference between two subsequent objects:stop): ",range(percent_amplitude_start,percent_amplitude_stop,percent_amplitude_length))
end


using Printf
using Dates
using HDF5

include("Outputs.jl")


function modelrun()
    options = Dict()   
    top_dir = mk_modelrun_dir()
    nlambda = wavelength_length
    nhice = ice_length
    namp = percent_amplitude_length
    local lambda = range(wavelength_start,wavelength_stop,nlambda)
    local hice =  range(ice_start,ice_stop,nhice)
    local per_amp = range(percent_amplitude_start,percent_amplitude_stop,namp)
    irun = 1
    for k in 1:namp
        t_halfspace = zeros(nlambda,nhice)
        t_rel = zeros(nlambda,nhice)
        t_tic = zeros(nlambda,nhice)
        amplitude_percentage = per_amp[k]
        amp_decimal = amplitude_percentage/100
        sub_dir = mk_sub_dir(top_dir,amplitude_percentage)
        for i in 1:nlambda
            for j in 1:nhice
                options["wavelength"] = lambda[i]*1e3
                options["ice thickness"] = hice[j]*1e3
                options["amplitude"] = amp_decimal*options["ice thickness"]
                options["surface depth"] = options["amplitude"] 
                # if options["wavelength"] >= options["ice thickness"]
                    println("Starting model execution for model run $irun...")
                    sub_dir_by_run,sub_dir_plots,sub_dir_data = mk_output_dir(sub_dir,irun)
                    test_data_info(ice_start,ice_stop,nhice,wavelength_start,wavelength_stop,nlambda,amplitude_percentage)
                    # println("Using Wavelength: ",options["wavelength"]/1e3,"(km)"," , ","Using Ice Shell Thickness: ",options["ice thickness"]/1e3,"(km)"," , ","Using Amplitde Percentage: $amplitude_percentage%")
                    open(sub_dir_by_run*"/output.txt", "w") do out
                        redirect_stdout(out) do
                            model_runtime(sub_dir_by_run,"Start")
                            print("Hello World!")
                            model_runtime(sub_dir_by_run,"End")
                            amplitude = rand(Float64,10)
                            time_thickening = rand(Float64,10)
                            amplitude = copy(amplitude)
                            time_thickening = copy(time_thickening)
                            data_thickening_time(amplitude,time_thickening,sub_dir_by_run)
                        end
                    end
                    println("Model ran successfully for model run $irun. Outputs saved to output.txt")
                    irun += 1
                # end
            end
        end
        lambda = vcat(map(x->x',lambda)...)
        hice = vcat(map(x->x',hice)...)
        h5open(sub_dir*"/Data.hdf5", "w") do file
            println("Creating HDF5 File")
            println("Saving Data into a HDF5 File")
            # Creating Groups for Data
            g = create_group(file, "Model Run")
            # Storing Data Inside the Group
            g["Wavelength"] = lambda[:]
            g["Ice Shell Thickness"] = hice[:]
            # Apply an Attribute to Groups
            attrs(g)["Description"] = "This group contains only a 4 dataset"
            println("Finished Saving Data into a HDF5 File")
        end
        if namp != 1
            irun = 1
        end
    end
end


modelrun()
