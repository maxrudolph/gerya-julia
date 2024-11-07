# Import necessary module
using InteractiveUtils
using Printf
using Dates

function mk_modelrun_dir()
    day_time_stamp = Dates.format(now(),"mm-dd-YYYY_HH:MM:SS")
    dir_name = joinpath(@__DIR__,"Model_Outputs","ModelRun$day_time_stamp")
    mkpath(dir_name)
    return dir_name
end

function mk_output_dir(top_dir::String,irun::Int64)
    sub_dir = mkdir(joinpath(top_dir,"Run_$irun"))
    sub_dir_plots = mkdir(joinpath(sub_dir,"Plots"));
    sub_dir_data = mkdir(joinpath(sub_dir,"Data"));
    return sub_dir,sub_dir_plots,sub_dir_data
end

function mk_sub_dir(top_dir::String,amplitude_percentage::Float64)
    sub_dir = mkdir(joinpath(top_dir,"Amplitude_$amplitude_percentage%"))
    return sub_dir
end

function model_runtime(sub_dir::String,time_type::String)
    if time_type == "Start"
        daystamp = Dates.format(now(),"mm/dd/YYYY")
        touch(sub_dir*"/Model_Runtime.txt")
        # Open file in write mode
        start_model = open(sub_dir*"/Model_Runtime.txt","w")
        # Writing to a file using write() method
        write(start_model,"Date: $daystamp\n")
        write(start_model,"Model Runtime:\n")
        start_time = Dates.format(now(),"HH:MM:SS")
        write(start_model,"Start - $start_time\n")
        # Closing the file in order to write the content from the disk to file
        close(start_model)
    elseif time_type == "End"
        # Open file in append mode
        end_model =  open(sub_dir*"/Model_Runtime.txt","a")
        # Writing to a file using write() method
        end_time = Dates.format(now(),"HH:MM:SS");
        write(end_model,"Finished - $end_time")
        close(end_model)
    end
end

function data_to_hdf5_file(lambda::Any,hice::Any,t_halfspace::Matrix{Float64},t_rel::Matrix{Float64},t_tic::Matrix{Float64},top_dir::String)
    lambda = vcat(map(x->x',lambda)...)
    hice = vcat(map(x->x',hice)...)
    h5open(top_dir*"/data.hdf5", "w") do file
        println("Creating HDF5 File")
        println("Saving Data into a HDF5 File")
        # Creating Groups for Data
        g = create_group(file, "Model Run")
        # Storing Data Inside the Group
        g["Wavelength"] = lambda[:]
        g["Ice Shell Thickness"] = hice[:]
        g["Viscous Relaxation Time(Half-Space)"] = t_halfspace
        g["Viscous Relaxation Time(Model)"] = t_rel
        g["Thickening Time"] = t_tic
        # Apply an Attribute to Groups
        attrs(g)["Description"] = "This group contains only a 4 dataset"
        println("Finished Saving Data into a HDF5 File")
    end
end

function data_table_info(ice_start::Float64,ice_stop::Float64,ice_length::Int64,wavelength_start::Float64,wavelength_stop::Float64,wavelength_length::Int64,sub_dir::String,amplitude_percentage::Float64)
    lambda = range(wavelength_start,wavelength_stop,wavelength_length)
    hice =  range(ice_start,ice_stop,ice_length)  
    irun_mat = zeros(wavelength_length,ice_length)
    irun = 1
    open(sub_dir*"/DataInfo.txt","w") do out
        redirect_stdout(out) do
            println("------------------Paramaters Used For Each Model Run------------------")
            println("Model used a ice shell thickness for a range of $ice_start(km) to $ice_stop(km) for a length of $ice_length.")
            println("Model used a wavelength of topography for the ocean-ice interface for a range of $wavelength_start(km) to $wavelength_stop(km) for a length of $wavelength_length.")
            println("Model used a amplitude percentage of $amplitude_percentage%.")
            println("Table:")
            println("Model Run\t| Ice Thickness(km)\t| Wavelength(km)\t| Amplitude Percentage")
            for i in 1:wavelength_length 
                for j in 1:ice_length
                    # if lambda[i] >= hice[j]
                        print(irun)
                        print("\t\t\t\t")
                        # print(hice[j])
                        @printf(out,"%4.3g",hice[j])
                        print("\t\t\t\t")
                        # print(lambda[i])
                        @printf(out,"%4.3g",lambda[i])
                        print("\t\t\t\t")
                        print(amplitude_percentage)
                        print("\n")
                        irun_mat[i,j] = irun
                        irun += 1
                    # end
                end
            end
            println("------------------Info------------------")
            println("Data that are in a matrix will correspond to how the following table matrix is setup below.")     
            println("Number of rows: ",wavelength_length)
            println("Number of columns: ",ice_length)
            println("Each value represent the model run number, look at the model run data to see what parameters where used in that run as follows:")
            for row in eachrow(irun_mat)
                println(row)
            end
        end
    end
end

### The following code will prompt the user to write a short description about the purpose of the model runs ###
## starts here ##
function prompt_and_get_model_description(top_dir::String)
    while true
        println("Please write a short decription about the purpose of the model run(s):")
        description = readline()
        println("You entered the following description:")
        println(description)
        println("Is this correct? (y/n)")
        confirmation = readline()

        if lowercase(confirmation) == "y"
            file_name = top_dir*"/model_runs_info.txt"
            open(file_name,"w") do file
                write(file,description)
            end
            break  # Exit the loop if description is confirmed
        elseif lowercase(confirmation) == "n"
            println("Lets try again.")
        else 
            println("Invalid input. Please enter 'y' for yes or 'n' for no.")
        end
    end
end
## ends here ##