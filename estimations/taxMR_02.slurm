#!/bin/bash

#SBATCH -p cascade
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32768
#SBATCH --mail-user=twozniak@unimelb.edu.au
#SBATCH --mail-type=ALL
#SBATCH --job-name="cgMR_02"

module load R/4.2.1
Rscript taxMR_02.R