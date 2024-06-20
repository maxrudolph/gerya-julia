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
#     # fitted_time is in yr
#     # fitted_amp is in meters
#     return fitted_time,fitted_amp
# end

########## Output for Plots ##########
function get_plots(grid::CartesianGrid,S::Matrix{Float64},T::Matrix{Float64},X::Matrix{Float64},time_type::String)
    if time_type == "initial"
        figure()
        title(L"Initial\,Melt\,Fraction")
        cs = contour(grid.xc/1e3,grid.yc/1e3,X,[0.5],colors="red")
        clabel(cs,inline=true,fontsize=8,fmt="ocean-ice interface")
        pcolor(grid.xc/1e3,grid.yc/1e3,X)
        colorbar(cmap="viridis")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        tight_layout()
        show()
        
        figure()
        title(L"Initial\,Temperature")
        cs = contour(grid.xc/1e3,grid.yc/1e3,X,[0.5],colors="red")
        clabel(cs,inline=true,fontsize=8,fmt="ocean-ice interface")
        pcolor(grid.xc/1e3,grid.yc/1e3,T)
        colorbar(cmap="viridis")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        tight_layout()
        show()
        
        figure()
        title(L"Initial\,Entropy")
        cs = contour(grid.xc/1e3,grid.yc/1e3,X,[0.5],colors="red")
        clabel(cs,inline=true,fontsize=8,fmt="ocean-ice interface")
        pcolor(grid.xc/1e3,grid.yc/1e3,S)
        colorbar(cmap="viridis")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        tight_layout()
        show()
    elseif time_type == "final"
        figure()
        title(L"Final\,Melt\,Fraction")
        cs = contour(grid.xc/1e3,grid.yc/1e3,X,[0.5],colors="red")
        clabel(cs,inline=true,fontsize=8,fmt="ocean-ice interface")
        pcolor(grid.xc/1e3,grid.yc/1e3,X)
        colorbar(cmap="viridis")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        tight_layout()
        show() 
        
        figure()
        title(L"Final\,Temperature")
        cs = contour(grid.xc/1e3,grid.yc/1e3,X,[0.5],colors="red")
        clabel(cs,inline=true,fontsize=8,fmt="ocean-ice interface")
        pcolor(grid.xc/1e3,grid.yc/1e3,T)
        colorbar(cmap="viridis")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        tight_layout()
        show() 
        
        figure()
        title(L"Final\,Entropy")
        cs = contour(grid.xc/1e3,grid.yc/1e3,X,[0.5],colors="red")
        clabel(cs,inline=true,fontsize=8,fmt="ocean-ice interface")
        pcolor(grid.xc/1e3,grid.yc/1e3,S)
        colorbar(cmap="viridis")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        tight_layout()
        show()
    end
end

function thickness_over_time(grid::CartesianGrid,thickness_array::Vector{Any},time_plot::Vector{Any},itime::Int64)
    figure()
    for i in 1:ceil(Int,itime/10):itime-1
        plot(grid.xc/1e3,thickness_array[i]/1e3,label=(L"At",@sprintf("%.3g",time_plot[i]/3.15e7/1e6),L"Myr"))
    end
    title(L"Profile\,of\,Ocean-Ice\,Inferface\,Topograpgy\,Over\,Time")
    gca().invert_yaxis()
    # Legend is at the bottom
    legend(loc="upper center",bbox_to_anchor=(0.5,-0.15),ncol=5)
    show()
end

function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},time::Float64,itime::Int64,sub_dir_plots::String)
# function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},times::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
        # topography::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
    x_time = @sprintf("%.3g",time/3.15e7/1e3)
    # Inital Model Schematic Profile
    figure() 
    pcolor(grid.xc/1000,grid.yc/1000,i_mat)
    colorbar(cmap="viridis")
    # plot(grid.xc/1000,i_interface_1/1000,"m",label="air-ice interface")
    # plot(grid.xc/1000,i_interface_2/1000,"r",label="ocean-ice interface")
    plot(grid.xc[2:end-1]/1000,i_interface_1/1000,"m",label="air-ice interface")
    plot(grid.xc[2:end-1]/1000,i_interface_2/1000,"r",label="ocean-ice interface")
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
    # plot(grid.xc/1000,interface_1/1000,"m",label="air-ice interface")
    # plot(grid.xc/1000,interface_2/1000,"r",label="ocean-ice interface")
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

function data_to_hdf5_file(lambda::Any,hice::Any,t_halfspace::Matrix{Float64},t_rel::Matrix{Float64},t_tic::Matrix{Float64},top_dir::String)
# function data_to_hdf5_file(lambda::Any,hice::Any,t_halfspace::Vector{Any},t_rel::Vector{Any},t_tic::Vector{Any},t_rel_fitted_time::Vector{Any},t_rel_fitted_amp::Vector{Any},top_dir::String)
    lambda = vcat(map(x->x',lambda)...)
    hice = vcat(map(x->x',hice)...)
    # t_halfspace = vcat(map(x->x',t_halfspace)...)
    # t_vis = vcat(map(x->x',t_vis)...)
    # t_tic = vcat(map(x->x',t_tic)...)
    # t_vis_fitted_time = vcat(map(x->x',t_vis_fitted_time)...)
    # t_vis_fitted_amp = vcat(map(x->x',t_vis_fitted_amp)...)
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
        # g["Fitted Viscous Relaxation Time"] = t_rel_fitted_time
        # g["Viscous Relaxation Time with Fitted Amplitude(Model)"] = t_rel_fitted_amp
        # Apply an Attribute to Groups
        attrs(g)["Description"] = "This group contains only a 4 dataset"
        println("Finished Saving Data into a HDF5 File")
    end
end