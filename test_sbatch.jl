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
    # Define ice shell thickness range (adjust as needed) in km
    ice_start = 10.0
    ice_stop = 60.0
    nice = 2
    ice_shell_thickness_range = range(ice_start,ice_stop,nice)
    # Define wavelength range (adjust as needed) in km
    wavelength_start = 5.0
    wavelength_stop = 200.0 
    nwavelength = 3
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
                main_dir = mk_main_dir(h,lambda,amp)
                # Constructing the command string
                cmd = `julia VMS_sbatch.jl $h $lambda $amp $main_dir`
                # Executing the command
                execute_cmd(cmd)
            end
        end
    end
end
# Executing main function
main_script()
