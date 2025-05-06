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
julia --threads 64 -- PlumeModel_old_10.jl 150000.0 450.0 120000.5
#julia InstallDependencies.jl
