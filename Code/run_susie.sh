#!/bin/bash

#SBATCH --mem=240g     
#SBATCH --cpus-per-task=2

module load R/3.4.4

# type the number only of the chromosome to process
Rscript /gpfs/commons/home/mschertzer/ipsc_sqtl/Code/susie_for_UKB_LD.R $1