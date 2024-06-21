########## Model Plots ##########
function get_plots(grid::CartesianGrid,S::Matrix{Float64},T::Matrix{Float64},X::Matrix{Float64},time_type::String,plot_dir::String)
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
	savefig(plot_dir*"/initial_melt_fraction.png")
        close()
        
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
	savefig(plot_dir*"/initial_temperature.png")
        close()
        
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
	savefig(plot_dir*"/initial_entropy.png")
        close()

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
	savefig(plot_dir*"/final_melt_fraction.png")
        close() 
        
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
	savefig(plot_dir*"/final_temperature.png")
        close() 
        
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
	savefig(plot_dir*"/final_entropy.png")
        close()
    end
end

function get_plots_new(grid::CartesianGrid,S::Matrix{Float64},T::Matrix{Float64},X::Matrix{Float64},time_type::String,plot_dir::String)
    if time_type == "initial"
        fig, axs = subplots(1, 3, figsize=(15, 5))  # Create a figure with 1 row and 3 columns of subplots
        fig.subplots_adjust(wspace=0.3)  # Adjust the horizontal space between subplots

        # Plot initial melt fraction
        axs[1].set_title(L"Initial\,Melt\,Fraction")
        cs = contour(axs[1], grid.xc/1e3, grid.yc/1e3, X, [0.5], colors="red")
        clabel(cs, inline=true, fontsize=8, fmt="ocean-ice interface")
        pcm1 = axs[1].pcolor(grid.xc/1e3, grid.yc/1e3, X)
        colorbar(pcm1, ax=axs[1], cmap="viridis")
        axs[1].set_ylabel(L"Height\,(km)")
        axs[1].set_xlabel(L"Width\,(km)")
        axs[1].invert_yaxis()
        axs[1].set_aspect("equal")

        # Plot initial temperature
        axs[2].set_title(L"Initial\,Temperature")
        cs = contour(axs[2], grid.xc/1e3, grid.yc/1e3, X, [0.5], colors="red")
        clabel(cs, inline=true, fontsize=8, fmt="ocean-ice interface")
        pcm2 = axs[2].pcolor(grid.xc/1e3, grid.yc/1e3, T)
        colorbar(pcm2, ax=axs[2], cmap="viridis")
        axs[2].set_ylabel(L"Height\,(km)")
        axs[2].set_xlabel(L"Width\,(km)")
        axs[2].invert_yaxis()
        axs[2].set_aspect("equal")

        # Plot initial entropy
        axs[3].set_title(L"Initial\,Entropy")
        cs = contour(axs[3], grid.xc/1e3, grid.yc/1e3, X, [0.5], colors="red")
        clabel(cs, inline=true, fontsize=8, fmt="ocean-ice interface")
        pcm3 = axs[3].pcolor(grid.xc/1e3, grid.yc/1e3, S)
        colorbar(pcm3, ax=axs[3], cmap="viridis")
        axs[3].set_ylabel(L"Height\,(km)")
        axs[3].set_xlabel(L"Width\,(km)")
        axs[3].invert_yaxis()
        axs[3].set_aspect("equal")

        tight_layout()

        savefig(plot_dir*"/initial_plots.png")  # Save the figure
        close()  # Close the figure to free up resources

    elseif time_type == "final"
        fig, axs = subplots(1, 3, figsize=(15, 5))  # Create a figure with 1 row and 3 columns of subplots
        fig.subplots_adjust(wspace=0.3)  # Adjust the horizontal space between subplots

        # Plot final melt fraction
        axs[1].set_title(L"Final\,Melt\,Fraction")
        cs = contour(axs[1], grid.xc/1e3, grid.yc/1e3, X, [0.5], colors="red")
        clabel(cs, inline=true, fontsize=8, fmt="ocean-ice interface")
        pcm1 = axs[1].pcolor(grid.xc/1e3, grid.yc/1e3, X)
        colorbar(pcm1, ax=axs[1], cmap="viridis")
        axs[1].set_ylabel(L"Height\,(km)")
        axs[1].set_xlabel(L"Width\,(km)")
        axs[1].invert_yaxis()
        axs[1].set_aspect("equal")

        # Plot final temperature
        axs[2].set_title(L"Final\,Temperature")
        cs = contour(axs[2], grid.xc/1e3, grid.yc/1e3, X, [0.5], colors="red")
        clabel(cs, inline=true, fontsize=8, fmt="ocean-ice interface")
        pcm2 = axs[2].pcolor(grid.xc/1e3, grid.yc/1e3, T)
        colorbar(pcm2, ax=axs[2], cmap="viridis")
        axs[2].set_ylabel(L"Height\,(km)")
        axs[2].set_xlabel(L"Width\,(km)")
        axs[2].invert_yaxis()
        axs[2].set_aspect("equal")

        # Plot final entropy
        axs[3].set_title(L"Final\,Entropy")
        cs = contour(axs[3], grid.xc/1e3, grid.yc/1e3, X, [0.5], colors="red")
        clabel(cs, inline=true, fontsize=8, fmt="ocean-ice interface")
        pcm3 = axs[3].pcolor(grid.xc/1e3, grid.yc/1e3, S)
        colorbar(pcm3, ax=axs[3], cmap="viridis")
        axs[3].set_ylabel(L"Height\,(km)")
        axs[3].set_xlabel(L"Width\,(km)")
        axs[3].invert_yaxis()
        axs[3].set_aspect("equal")

        tight_layout()

        savefig(plot_dir*"/final_plots.png")  # Save the figure
        close()  # Close the figure to free up resources

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