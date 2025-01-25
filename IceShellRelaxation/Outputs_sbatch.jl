# Import necessary module
using Printf

function mk_main_dir(hice::Float64,lambda::Float64,amp::Float64,g::Float64)
    dir_name = joinpath(@__DIR__,"Model_Outputs","h_$hice"*"_lambda_$lambda"*"_amp_$amp"*"_g_$g")
    mkpath(dir_name)
    return dir_name
end

function mk_sub_dir(main_dir::String)
    sub_dir_plots = mkdir(joinpath(main_dir,"Plots"))
    sub_dir_data = mkdir(joinpath(main_dir,"Visualization"))
    return sub_dir_plots,sub_dir_data
end

function data_table_info(ice_start::Float64,ice_stop::Float64,ice_length::Int64,wavelength_start::Float64,wavelength_stop::Float64,wavelength_length::Int64,dir::String,amplitude_percentage::Float64)
    lambda = range(wavelength_start,wavelength_stop,wavelength_length)
    hice =  range(ice_start,ice_stop,ice_length)  
    irun_mat = zeros(wavelength_length,ice_length)
    irun = 1
    open(dir*"/DataInfo.txt","w") do out
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
