# Function to get the directory path for a given combination of parameters
function mk_main_dir(hice::Float64,lambda::Float64,amp::Float64,g::Float64)
    dir_name = joinpath(@__DIR__,"Model_Outputs","h_$hice"*"_lambda_$lambda"*"_amp_$amp"*"_g_$g")
    return dir_name
end    

# Function to get the HDF5 file path for a given combination of parameters
function get_hdf5_file_path(hice::Float64,lambda::Float64,amp::Float64,g::Float64)
    main_dir = mk_main_dir(hice,lambda,amp,g)
    file_path = joinpath(main_dir,"data.hdf5")  # Assuming the HDF5 file is named "data.hdf5"
    return file_path
end

function combine_hdf5_files(ice_shell_thickness_range::AbstractRange{Float64},wavelength_range::AbstractRange{Float64},amplitude::Float64,g::Float64,output_path::String)
    # Collect unique values
    wavelength_set = Set{Float64}()
    ice_shell_thickness_set = Set{Float64}()
    
    combined_data = Dict(
        "Wavelength" => Float64[],
        "Ice Shell Thickness" => Float64[],
        "Viscous Relaxation Time(Half-Space)" => Float64[],
        "Viscous Relaxation Time(Model)" => Float64[],
        "Fitted Viscous Relaxation Time" => Float64[],
        "Thickening Time" => Float64[],
        "Fitted Thickening Time" => Float64[]
    )

    # Loop over the range of ice shell thickness and wavelength values
    for h in ice_shell_thickness_range
        for lambda in wavelength_range
            file_path = get_hdf5_file_path(h,lambda,amplitude,g)
            if isfile(file_path)
                h5open(file_path, "r") do file
                    g = file["Model Run"]
                    append!(combined_data["Wavelength"],read(g["Wavelength"]))
                    append!(combined_data["Ice Shell Thickness"],read(g["Ice Shell Thickness"]))
                    append!(combined_data["Viscous Relaxation Time(Half-Space)"],read(g["Viscous Relaxation Time(Half-Space)"]))
                    append!(combined_data["Viscous Relaxation Time(Model)"],read(g["Viscous Relaxation Time(Model)"]))
                    append!(combined_data["Fitted Viscous Relaxation Time"],read(g["Fitted Viscous Relaxation Time"]))
                    append!(combined_data["Thickening Time"],read(g["Thickening Time"])) 
                    append!(combined_data["Fitted Thickening Time"],read(g["Fitted Thickening Time"]))
                    # Add to values to wavelength and ice shell thickness values
                    union!(wavelength_set,read(g["Wavelength"]))
                    union!(ice_shell_thickness_set,read(g["Ice Shell Thickness"]))
                end
            end
        end
    end
    
    sorted_wavelength = sort(collect(wavelength_set))
    sorted_hice = sort(collect(ice_shell_thickness_set))

    n_wavelength = length(sorted_wavelength)
    n_hice = length(sorted_hice)
    
    t_hs_matrix = zeros(n_wavelength,n_hice)
    t_rel_matrix = zeros(n_wavelength,n_hice)
    t_rel_fit_matrix = zeros(n_wavelength,n_hice)
    t_thic_matrix = zeros(n_wavelength,n_hice)
    t_thic_fit_matrix = zeros(n_wavelength,n_hice)

    for i in 1:length(combined_data["Viscous Relaxation Time(Half-Space)"])
        wavelength = combined_data["Wavelength"][i]
        hice = combined_data["Ice Shell Thickness"][i]   
        # Find the correct indices for wavelength and hice
        j_idx = findfirst(x -> x == hice, sorted_hice)
        i_idx = findfirst(x -> x == wavelength, sorted_wavelength)
        # Fill the matrices
        t_hs_matrix[i_idx,j_idx] = combined_data["Viscous Relaxation Time(Half-Space)"][i]
        t_rel_matrix[i_idx,j_idx] = combined_data["Viscous Relaxation Time(Model)"][i]
        t_rel_fit_matrix[i_idx,j_idx] = combined_data["Fitted Viscous Relaxation Time"][i]
        t_thic_matrix[i_idx,j_idx] = combined_data["Thickening Time"][i]
        t_thic_fit_matrix[i_idx,j_idx] = combined_data["Fitted Thickening Time"][i]
    end
    
    h5open(output_path,"w") do file
        g = create_group(file,"Combined Model Run")
        g["Wavelength"] = sorted_wavelength
        g["Ice Shell Thickness"] = sorted_hice
        g["Viscous Relaxation Time(Half-Space)"] = t_hs_matrix
        g["Viscous Relaxation Time(Model)"] = t_rel_matrix
        g["Fitted Viscous Relaxation Time"] = t_rel_fit_matrix
        g["Thickening Time"] = t_thic_matrix
        g["Fitted Thickening Time"] = t_thic_fit_matrix
        attrs(g)["Description"] = "This group contains combined and sorted unique datasets"
        println("Finished Saving Data into a HDF5 File")
    end
end