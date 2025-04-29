#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=plume
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=ayylu@ucdavis.edu
#SBATCH --partition=high2
#SBATCH --time=6-00:00:00
#source /home/qxl51/miniconda3/bin/activate
which python3
which julia
julia --threads 64 -- PlumeModel_old.jl 150000.0 400.0 60000.3
#julia InstallDependencies.jl
