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

function mk_output_dir(sub_dir::String,irun::Int64)
    sub_dir_by_run = mkdir(joinpath(sub_dir,"Run_$irun"))
    sub_dir_plots = mkdir(joinpath(sub_dir_by_run,"Plots"));
    sub_dir_data = mkdir(joinpath(sub_dir_by_run,"Data"));
    return sub_dir_by_run,sub_dir_plots,sub_dir_data
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

########## Output for Fitted Data ##########
# function topo_fitting_data(topography::Vector{Any},times::Vector{Any},itime::Int64,sub_dir_plots::String)
#     topofit = []
#     for i in 1:itime-1
#         amp = maximum(topography[i])-minimum(topography[i])
#         append!(topofit,amp)
#     end
#     x = convert(Array{Float64}, times)
#     y = convert(Array{Float64}, topofit)
#     fit = fitexp(x/3.15e7,y,options=Options(besttol=100))
#     # figure()
#     for i in 1:itime-1
#         plot(times[i]/3.15e7,maximum(topography[i])-minimum(topography[i]),".") 
#     end
#     plot(fit.x,fit.y,"r-",label="fitted data")
#     gca().set_ylabel(L"Amplitude\,(m)")
#     gca().set_xlabel(L"Time\,(years)")
#     legend(fancybox="True",shadow="True")
#     savefig(sub_dir_plots*"/Fitted_Topo_Data.pdf")
#     close()
#     fitted_time = fit.b
#     return fitted_time
# end

# function amp_fitting_data(amplitude::Vector{Any},times::Vector{Any},itime::Int64,sub_dir_plots::String)
#     ampfit = []
#     for i in 1:itime-1
#         amp = amplitude[i]
#         append!(ampfit,amp)
#     end
#     xx = convert(Array{Float64}, times)
#     yy = convert(Array{Float64}, ampfit)
#     fitt = fitexp(xx/3.15e7,yy,options=Options(besttol=100))
#     # figure()
#     for i in 1:itime-1
#         plot(times[i]/3.15e7,ampfit[i],".") 
#     end
#     plot(fitt.x,fitt.y,"r-",label="fitted data")
#     gca().set_ylabel(L"Amplitude\,(m)")
#     gca().set_xlabel(L"Time\,(years)")
#     legend(fancybox="True",shadow="True")
#     savefig(sub_dir_plots*"/Fitted_Amp_Data.pdf")
#     close()
#     fitted_time = fitt.b
#     fitted_amp = last(fitt.y)
#     return fitted_time,fitted_amp
# end

########## Output for Plots ##########
function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},times::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
        # topography::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
    x_time = @sprintf("%.3g",time/3.15e7/1e3)
    # Inital Model Schematic Profile
    figure() 
    pcolor(grid.xc/1000,grid.yc/1000,i_mat)
    colorbar(cmap="viridis")
    plot(grid.xc[2:end-1]/1000,i_interface_1/1000,"m",label="air-ice interface")
    plot(grid.xc[2:end-1]/1000,i_interface_2/1000,"r",label="air-ice interface")
    title(L"Initial\,\,Model\,\,Schematic")
    gca().invert_yaxis()
    gca().set_aspect("equal")
    gca().set_ylabel(L"Height\,(km)")
    gca().set_xlabel(L"Width\,(km)")
    # Legend is at the bottom
    legend(loc="upper center", bbox_to_anchor=(0.5,-0.15),fancybox="True",shadow="True",ncol=5)
    savefig(sub_dir_plots*"/Initial_Model_Schematic.pdf")
    close()

    # # Initial ice shell thickness profile along the horizontal direction
    # figure() 
    # plot(grid.xc/1000,(i_interface_2-i_interface_1)/1000)
    # title(L"Initial\,\,Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction")
    # gca().invert_yaxis()
    # gca().set_ylabel(L"thickness\,(km)")
    # gca().set_xlabel(L"x\,(km)")
    # savefig(sub_dir_plots*"/Initial_Ice_Thickness.pdf")
    # close()

    # Model Schematic Profile
    figure()
    pcolor(grid.xc/1000,grid.yc/1000,mat)
    colorbar(cmap="viridis")
    plot(grid.xc[2:end-1]/1000,interface_1/1000,"m",label="air-ice interface")
    plot(grid.xc[2:end-1]/1000,interface_2/1000,"r",label="ocean-ice interface")
    gca().invert_yaxis()
    gca().set_aspect("equal")
    gca().set_ylabel(L"Height\,(km)")
    gca().set_xlabel(L"Width\,(km)")
    title(L"Final\,\,Model\,\,Schematic")
    # Legend is at the bottom
    legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5) 
    savefig(sub_dir_plots*"/Final_Model_Schematic.pdf")
    close()

    # # Ice shell thickness profile along the horizontal direction
    # figure()
    # plot(grid.xc/1000,(interface_2-interface_1)/1000)
    # title(L"Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction\,\,at\,\,%$x_time\,kyr")
    # gca().invert_yaxis()
    # gca().set_ylabel(L"thickness\,(km)")
    # gca().set_xlabel(L"x\,(km)")
    # savefig(sub_dir_plots*"/Final_Ice_Thickness.pdf")
    # close()

    # Profile of ice-water inferface topograpgy over time
    # figure()
    # for i in 1:5:itime-1
    #     plot(grid.xc[2:end-1],topography[i])
    # end
    # title(L"Profile\,\,of\,\,Ice-Water\,\,Interface\,\,Topography\,\,Over\,\,Time")
    # gca().invert_yaxis()
    # #     gca().set_xlim([0.0,1e4])
    # #     gca().set_ylim([1.8e4,2.2e4])
    # # Legend is at the bottom
    # # legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5)
    # savefig(sub_dir_plots*"/Ice_Water_Inferface_Topography.pdf")            
    # close()

    # # Profile of maximum topograpgy over time 
    # figure()
    # plot(times/3.15e7/1e3,max_topo/1000)
    # title(L"Maximum\,\,Topography\,\,After\,\,%$x_time\,ka")
    # gca().set_ylabel(L"Max.\,topography\,(km)")
    # gca().set_xlabel(L"Time\,(ka)")
    # show()
end


# function amp_fitting_data(amplitude::Vector{Any},times::Vector{Any},itime::Int64,sub_dir_plots::String)
#     ampfit = []
#     for i in 1:itime-1
#         amp = amplitude[i]
#         append!(ampfit,amp)
#     end
#     xx = convert(Array{Float64}, times)
#     yy = convert(Array{Float64}, ampfit)
#     fitt = fitexp(xx/3.15e7,yy,options=Options(besttol=100))
#     # figure()
#     for i in 1:itime-1
#         plot(times[i]/3.15e7,ampfit[i],".") 
#     end
#     plot(fitt.x,fitt.y,"r-",label="fitted data")
#     gca().set_ylabel(L"Amplitude\,(m)")
#     gca().set_xlabel(L"Time\,(years)")
#     legend(fancybox="True",shadow="True")
#     savefig(sub_dir_plots*"/Fitted_Amp_Data.pdf")
#     close()
#     fitted_time = fitt.b
#     fitted_amp = last(fitt.y)
#     # fitted_time is in yr
#     # fitted_amp is in meters
#     return fitted_time,fitted_amp
# end

########## Output for Plots ##########
# function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},time::Float64,itime::Int64,sub_dir_plots::String)
# # function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},times::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
#         # topography::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
#     x_time = @sprintf("%.3g",time/3.15e7/1e3)
#     # Inital Model Schematic Profile
#     figure() 
#     pcolor(grid.xc/1000,grid.yc/1000,i_mat)
#     colorbar(cmap="viridis")
#     # plot(grid.xc/1000,i_interface_1/1000,"m",label="air-ice interface")
#     # plot(grid.xc/1000,i_interface_2/1000,"r",label="ocean-ice interface")
#     plot(grid.xc[2:end-1]/1000,i_interface_1/1000,"m",label="air-ice interface")
#     plot(grid.xc[2:end-1]/1000,i_interface_2/1000,"r",label="ocean-ice interface")
#     title(L"Initial\,\,Model\,\,Schematic")
#     gca().invert_yaxis()
#     gca().set_aspect("equal")
#     gca().set_ylabel(L"Height\,(km)")
#     gca().set_xlabel(L"Width\,(km)")
#     # Legend is at the bottom
#     legend(loc="upper center", bbox_to_anchor=(0.5,-0.15),fancybox="True",shadow="True",ncol=5)
#     savefig(sub_dir_plots*"/Initial_Model_Schematic.pdf")
#     close()

#     # # Initial ice shell thickness profile along the horizontal direction
#     # figure() 
#     # plot(grid.xc/1000,(i_interface_2-i_interface_1)/1000)
#     # title(L"Initial\,\,Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction")
#     # gca().invert_yaxis()
#     # gca().set_ylabel(L"thickness\,(km)")
#     # gca().set_xlabel(L"x\,(km)")
#     # savefig(sub_dir_plots*"/Initial_Ice_Thickness.pdf")
#     # close()

#     # Model Schematic Profile
#     figure()
#     pcolor(grid.xc/1000,grid.yc/1000,mat)
#     colorbar(cmap="viridis")
#     # plot(grid.xc/1000,interface_1/1000,"m",label="air-ice interface")
#     # plot(grid.xc/1000,interface_2/1000,"r",label="ocean-ice interface")
#     plot(grid.xc[2:end-1]/1000,interface_1/1000,"m",label="air-ice interface")
#     plot(grid.xc[2:end-1]/1000,interface_2/1000,"r",label="ocean-ice interface")
#     gca().invert_yaxis()
#     gca().set_aspect("equal")
#     gca().set_ylabel(L"Height\,(km)")
#     gca().set_xlabel(L"Width\,(km)")
#     title(L"Final\,\,Model\,\,Schematic")
#     # Legend is at the bottom
#     legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5) 
#     savefig(sub_dir_plots*"/Final_Model_Schematic.pdf")
#     close()

#     # # Ice shell thickness profile along the horizontal direction
#     # figure()
#     # plot(grid.xc/1000,(interface_2-interface_1)/1000)
#     # title(L"Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction\,\,at\,\,%$x_time\,kyr")
#     # gca().invert_yaxis()
#     # gca().set_ylabel(L"thickness\,(km)")
#     # gca().set_xlabel(L"x\,(km)")
#     # savefig(sub_dir_plots*"/Final_Ice_Thickness.pdf")
#     # close()

#     # Profile of ice-water inferface topograpgy over time
#     # figure()
#     # for i in 1:5:itime-1
#     #     plot(grid.xc[2:end-1],topography[i])
#     # end
#     # title(L"Profile\,\,of\,\,Ice-Water\,\,Interface\,\,Topography\,\,Over\,\,Time")
#     # gca().invert_yaxis()
#     # #     gca().set_xlim([0.0,1e4])
#     # #     gca().set_ylim([1.8e4,2.2e4])
#     # # Legend is at the bottom
#     # # legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5)
#     # savefig(sub_dir_plots*"/Ice_Water_Inferface_Topography.pdf")            
#     # close()

#     # # Profile of maximum topograpgy over time 
#     # figure()
#     # plot(times/3.15e7/1e3,max_topo/1000)
#     # title(L"Maximum\,\,Topography\,\,After\,\,%$x_time\,ka")
#     # gca().set_ylabel(L"Max.\,topography\,(km)")
#     # gca().set_xlabel(L"Time\,(ka)")
#     # show()
# end
    


function thickening_data_to_hdf5_file(lambda::Any,hice::Any,t_tic::Matrix{Float64},sub_dir::String)
    lambda = vcat(map(x->x',lambda)...)
    hice = vcat(map(x->x',hice)...)
    h5open(sub_dir*"/Thickening_Data.hdf5", "w") do file
        println("Creating HDF5 File")
        println("Saving Data into a HDF5 File")
        # Creating Groups for Data
        g = create_group(file, "Model Run")
        # Storing Data Inside the Group
        g["Wavelength"] = lambda[:]
        g["Ice Shell Thickness"] = hice[:]
        g["Thickening Time"] = t_tic
        # Apply an Attribute to Groups
        attrs(g)["Description"] = "This group contains only a 3 dataset"
        println("Finished Saving Data into a HDF5 File")
    end
end


function data_to_hdf5_file(lambda::Any,hice::Any,t_halfspace::Matrix{Float64},t_rel::Matrix{Float64},t_tic::Matrix{Float64},sub_dir::String)
    lambda = vcat(map(x->x',lambda)...)
    hice = vcat(map(x->x',hice)...)
    # t_halfspace = vcat(map(x->x',t_halfspace)...)
    # t_vis = vcat(map(x->x',t_vis)...)
    # t_tic = vcat(map(x->x',t_tic)...)
    h5open(sub_dir*"/Data.hdf5", "w") do file
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
        attrs(g)["Description"] = "This group contains only a 5 dataset"
        println("Finished Saving Data into a HDF5 File")
    end
end
    
function data_info(ice_start::Float64,ice_stop::Float64,ice_length::Int64,wavelength_start::Float64,wavelength_stop::Float64,wavelength_length::Int64,sub_dir::String,amplitude_percentage::Float64)
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

function test_data_info(ice_start::Float64,ice_stop::Float64,ice_length::Int64,wavelength_start::Float64,wavelength_stop::Float64,wavelength_length::Int64,amplitude_percentage::Float64)
    lambda = range(wavelength_start,wavelength_stop,wavelength_length)
    hice =  range(ice_start,ice_stop,ice_length)  
    irun_mat = zeros(wavelength_length,ice_length)
    irun = 1
    open("DataInfo.txt","w") do out
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