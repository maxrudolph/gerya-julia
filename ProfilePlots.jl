using PyPlot
function profileplot(time_type::String)
    x_time = @sprintf("%.3g",time/3.15e7/1e3)
    i_it_data = get_ice_thickness(i_air_ice_interface,i_ocean_ice_interface,"initial")
    it_data = get_ice_thickness(air_ice_interface,ocean_ice_interface,"after")
    if time_type == "initial"
        # Inital Model Schematic Profile
        figure() 
        pcolor(grid.xc/1000,grid.yc/1000,i_mat)
        colorbar(cmap="viridis")
        plot(grid.xc/1000,i_air_ice_interface/1000,"m",label="air-ice interface")
        plot(grid.xc/1000,i_ocean_ice_interface/1000,"r",label="ocean-ice interface")
        title(L"Initial\,\,Model\,\,Schematic")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        legend(loc="upper right",bbox_to_anchor=[1.8,1.0])
        show()

        # Initial ice shell thickness profile along the horizontal direction
        figure() 
        plot(grid.xc/1000,i_it_data[3]/1000)
        title(L"Initial\,\,Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction")
        gca().set_ylabel(L"thickness\,(km)")
        gca().set_xlabel(L"x\,(km)")
        show()
    elseif time_type == "after"
        # Model Schematic Profile
        figure()
        pcolor(grid.xc/1000,grid.yc/1000,mat)
        colorbar(cmap="viridis")
        plot(grid.xc/1000,air_ice_interface/1000,"m",label="air-ice interface")
        plot(grid.xc/1000,ocean_ice_interface/1000,"r",label="ocean-ice interface")
        title(L"Model\,\,Schematic\,\,at\,\,%$x_time\,Myr")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        legend(loc="upper right",bbox_to_anchor=[1.8,1.0])
        show()

        # Ice shell thickness profile along the horizontal direction
        figure()
        plot(grid.xc/1000,it_data[3]/1000)
        title(L"Ice\,\,Shell\,\,Thickness\,\,Across\,\,the\,\,Horizontal\,\,Direction\,\,at\,\,%$x_time\,Myr")
        gca().set_ylabel(L"thickness\,(km)")
        gca().set_xlabel(L"x\,(km)")
        show()

        # Profile of maximum topograpgy over time 
        figure()
        plot(time_plot/seconds_in_year/1e3,max_topo/1000)
        title(L"Maximum\,\,Topography\,\,Over\,\,Time")
        gca().set_ylabel(L"Max.\,topography\,(km)")
        gca().set_xlabel(L"Time\,(ka)")
        show()

        # Temperature Profile 
        figure()
        scatter(markers.x[1,:]/1000,markers.x[2,:]/1000,c=markers.scalars[markers.scalarFields["T"],:],s=0.1)
        plot(grid.xc/1000,air_ice_interface/1000,"m",label="air-ice interface")
        plot(grid.xc/1000,ocean_ice_interface/1000,"r",label="ocean-ice interface")
        title(L"Temperature\,\,Profile\,\,at\,\,%$x_time\,Myr")
        colorbar(label=L"K",cmap="viridis")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        legend(loc="upper right",bbox_to_anchor=[1.8,1.0])
        show()

        # Viscosity Profile 
        figure()
        scatter(markers.x[1,:]/1000,markers.x[2,:]/1000,c=log10.(markers.scalars[markers.scalarFields["eta"],:]),s=0.1)
        plot(grid.xc/1000,air_ice_interface/1000,"m",label="air-ice interface")
        plot(grid.xc/1000,ocean_ice_interface/1000,"r",label="ocean-ice interface")
        title(L"Viscosity\,\,Profile\,\,at\,\,%$x_time\,Myr")
        colorbar(label=L"Pa\cdot{s}",cmap="viridis")
        gca().invert_yaxis()
        gca().set_aspect("equal")
        gca().set_ylabel(L"Height\,(km)")
        gca().set_xlabel(L"Width\,(km)")
        legend(loc="upper right",bbox_to_anchor=[1.8,1.0])
        show()
    end
end