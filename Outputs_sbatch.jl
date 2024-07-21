function mk_main_dir(hice::Float64,lambda::Float64,amp::Float64)
    dir_name = joinpath(@__DIR__,"Model_Outputs","h_$hice"*"_lambda_$lambda"*"_amp_$amp")
    mkpath(dir_name)
    return dir_name
end

function mk_sub_dir(main_dir::String)
	  sub_dir_plots = mkdir(joinpath(main_dir,"Plots"))
    sub_dir_data = mkdir(joinpath(main_dir,"Visualization"))
    return sub_dir_plots,sub_dir_data
end
