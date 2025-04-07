#!/bin/bash

#SBATCH -p cascade
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32768
#SBATCH --mail-user=twozniak@unimelb.edu.au
#SBATCH --mail-type=ALL
#SBATCH --job-name="sv"
#SBATCH --array=701-1100

module load R/4.2.1
Rscript bsvar_sv.R ${SLURM_ARRAY_TASK_ID}
