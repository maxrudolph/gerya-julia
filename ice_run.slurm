#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=ice_model
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=caguirre@ucdavis.edu
#SBATCH --partition=high2
#SBATCH --time=6-00:00:00
source /home/caguirre/miniconda3/bin/activate
julia --threads 64 -- eat.jl 10.0 10.0 1 100.0 100.0 10.0 10.0 1
