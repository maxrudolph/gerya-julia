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
source /home/rudolph/miniconda3/bin/activate
which python3
which julia
julia --threads 64 -- EclPlume_1.jl 400.0 90000.0 600000.0
