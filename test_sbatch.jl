### Importing (using/include) packages/modules and files needed for the code to run ###
# The Base.Process module is needed to run commands
import Base.Process
include("Outputs_sbatch.jl")

function execute_cmd(cmd::Cmd)
    try
        println("Executing cmd: $cmd")
        run(cmd)
    catch e
        println("An error occurred: $e")
    end
end

function main_script()
    # Define g for icy moon
    g = 0.113
    # Define ice shell thickness range (adjust as needed) in km
    ice_start = 10.0
    ice_stop = 60.0
    nice = 15
    ice_shell_thickness_range = range(ice_start,ice_stop,nice)
    # Define wavelength range (adjust as needed) in km
    wavelength_start = 5.0
    wavelength_stop = 300.0
    nwavelength = 14
    wavelength_range = range(wavelength_start,wavelength_stop,nwavelength)
    # Define amplitude percentage range (adjust as needed)
    amplitude_start = amplitude_stop = 20.0
    amplitude_range = range(amplitude_start,amplitude_stop,1)
    # Making data table   	
    data_table_info(ice_start,ice_stop,nice,wavelength_start,wavelength_stop,nwavelength,"Model_Outputs",amplitude_start)
    # Lopping over the range of amplitude values
    for amp in amplitude_range
        # Looping over the range of ice shell thickness values
        for h in ice_shell_thickness_range
            # Looping over the range of wavelength values
            for lambda in wavelength_range
                main_dir = mk_main_dir(h,lambda,amp,g)
                # Constructing the command string
                cmd = `sbatch --account rudolphgrp --exclusive -N 1 --mem=128G --time=6-00:00:0 -p high2 --wrap "julia --threads 32 VMS_sbatch.jl $h $lambda $amp $g $main_dir"`
                # Executing the command
                execute_cmd(cmd)
            end
        end
    end
end
# Executing main function
main_script()
