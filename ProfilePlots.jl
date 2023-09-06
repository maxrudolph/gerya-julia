function get_topography_plots(grid::CartesianGrid,i_mat::Matrix{Float64},mat::Matrix{Float64},i_interface_1::Vector{Float64},interface_1::Vector{Float64},i_interface_2::Vector{Float64},interface_2::Vector{Float64},times::Vector{Any},topography::Vector{Any},time::Float64,itime::Int64,sub_dir_plots::String)
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
    figure()
    for i in 1:10:itime-1
        plot(grid.xc,topography[i])
            #,label=(L"At",@sprintf("%.3g",times[i]/3.15e7/1e3),L"kyr"))
    end
    title(L"Profile\,\,of\,\,Ice-Water\,\,Interface\,\,Topography\,\,Over\,\,Time")
    #     gca().set_xlim([0.0,1e4])
    #     gca().set_ylim([1.8e4,2.2e4])
    gca().invert_yaxis()
    # Legend is at the bottom
    # legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),fancybox="True",shadow="True", ncol=5)
    savefig(sub_dir_plots*"/Ice_Water_Inferface_Topography.pdf")            
    # close()

    # # Profile of maximum topograpgy over time 
    # figure()
    # plot(times/3.15e7/1e3,max_topo/1000)
    # title(L"Maximum\,\,Topography\,\,After\,\,%$x_time\,ka")
    # gca().set_ylabel(L"Max.\,topography\,(km)")
    # gca().set_xlabel(L"Time\,(ka)")
    # show()
end