#!/bin/bash

#SBATCH -p cascade
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32768
#SBATCH --mail-user=twozniak@unimelb.edu.au
#SBATCH --mail-type=ALL
#SBATCH --job-name="s_ce"
#SBATCH --array=401-700

module load R/4.2.1
Rscript bsvars_ce.R ${SLURM_ARRAY_TASK_ID}
