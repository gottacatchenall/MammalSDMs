#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=RodentSDMs 
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-165 
#SBATCH --time=3:00:00         

export JULIA_DEPOT_PATH="/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees"

module load julia/1.11.3
srun --unbuffered julia 04_fit_sdm.jl