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
    # Define ice shell thickness range (adjust as needed)
    ice_shell_thickness_range = range(5.0,10.0,2)
    # Define wavelength range (adjust as needed)
    wavelength_range = range(10.0,100.0,2)
    # Define amplitude percentage range (adjust as needed)
    amplitude_range = range(20.0,20.0,1)
    # Lopping over the range of amplitude values
    for amp in amplitude_range
        # Looping over the range of ice shell thickness values
        for h in ice_shell_thickness_range
            # Looping over the range of wavelength values
            for lambda in wavelength_range
                main_dir = mk_main_dir(h,lambda,amp)
                # Constructing the command string
                cmd = `julia TestOutput.jl $h $lambda $amp $main_dir`
                # Executing the command
                execute_cmd(cmd)
            end
        end
    end
end

# Executing main function
main_script()
