#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=4G

module load R
Rscript mod_bess2.R
