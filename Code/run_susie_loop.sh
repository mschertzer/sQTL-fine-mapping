#!/bin/bash

module load R/3.4.4

CORES=4

# this runs the finemapping code for all chromosomes
for chr in {1..22}
do
  echo $chr
  sbatch --cpus-per-task=4 --mem=240g --wrap "Rscript /gpfs/commons/home/mschertzer/ipsc_sqtl/code/susie_for_UKB_LD.R $chr $CORES TRUE"
done