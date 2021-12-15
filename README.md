# Fine Mapping from summary statistics using SuSie

This repository contains R code that loads summary statistics from QTL mapping for a single chromosome and then finds, loads, and subsets the corresponding reference LD matrix for SuSie fine mapping

### Resource links
This is the data that I used for fine mapping:
[iPSC summary statistics](https://zenodo.org/record/4005576#.YboFQPHMK3K)

Alkes group LD matrices from UKB:
[Reference LD matrices from UKB](https://alkesgroup.broadinstitute.org/UKBB_LD)

**Notes about these specific references** 
- Some LD files are missing (e.g. chr6_29000001_32000001) and will throw an error.
- Some npz files have the extension .npz2 instead of .npz and will throw an error. 
- These errors are associated with long-range LD regions, especially on human chr6. It is recommended to exclude these regions from analysis.

Vignette for Fine mapping with summary statistics using susieR:
[susieR vignette](https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html)

### Required R packages

1. require(reticulate)
2. require(Matrix)
3. require(doMC)
4. require(foreach)
5. require(readr)
6. require(tidyverse)
7. require(susieR)
8. require(Rfast)

### Additional information
- I ran this pipeline using R version 4.1.1 on a cluster.
- This script uses a lot of memory. I requested 240g which is the max for the cluster that I use. You can use 4 cores for most chromosomes but I did have to decrease to 1 core for chr6, 8, and 16. With 4 cores, a chromosome take ~5 hours and with 1 core, it takes ~24 hours.
- In susie_for_UKB_LD.R, change directory at line 40 based on the location of your reference LD files.
- In susie_for_UKB_LD.R, change directory at line 46 based on your summary stats location and file naming system.
- Example shell scripts for submitting one chromosome at a time (run_susie.sh) or submitting all human autosomes in a loop (run_susie_loop.sh) are available for reference.

### Usage

Rscript susie_for_UKB_LD.R 'chromosome' 'cores' 'error handling'
  
Example for chr6 with 5 cores: Rscript susie_for_UKB_LD.R 6 5 TRUE
  
