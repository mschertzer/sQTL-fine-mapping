#!/bin/bash

#SBATCH --mem=240g     
#SBATCH --cpus-per-task=1
#SBATCH -p pe2,bigmem

module load R/4.1.1

Rscript ~/ipsc_sqtl/Code/susie_for_UKB_LD.R $1 $2 $3