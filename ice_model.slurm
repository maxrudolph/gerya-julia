#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=ice_model
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=maxrudolph@ucdavis.edu
#SBATCH --partition=high2
#SBATCH --time=6-00:00:00
source /home/rudolph/miniconda3/bin/activate
which python3
which julia
julia --threads 64 -- dev_ice_model.jl 10.0 50.0 2 10.0 600.0 2

