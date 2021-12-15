#!/usr/bin/env Rscript

# load required packages
require(reticulate)
require(Matrix)
require(doMC)
require(foreach)
require(readr)
require(tidyverse)
require(susieR)
require(Rfast)
#require(optparse)

# command line arguments
# usage is: Rscript susie_for_UKB_LD.R <chromosome> <cores> <error handling>
# example for chr6 with 5 cores and : Rscript susie_for_UKB_LD.R 6 5 TRUE
# I recommend starting with TRUE for error handling
args = commandArgs(trailingOnly=TRUE)

# arguments <- parse_args(OptionParser(usage = "%prog [options]",
#                                      option_list=list(
#                                        make_option(c("--chrom"), default = 1, type = "double", help = "The chromosome to perform the fine mapping on. [Default %default]."),
#                                        make_option(c("--cores"), default = 1, type = "double", help = "Number of cores used. [Default %default]"),
#                                        make_option(c("--errorHandling"), default = TRUE, type = "logical", help = "This determines how the foreach loop handles errors. [Default %default], execution will be stopped, while FALSE will allow the loop to continue."))))

# function imports npz lower diagonal of reference LD matrix
read_ld = function(fn) {
  np <- import("numpy")
  npz <- np$load(fn)
  R_lower = sparseMatrix(i = npz$f[["row"]], 
                         j = npz$f[["col"]],
                         dims = npz$f[["shape"]], 
                         x = as.numeric(npz$f[["data"]]),
                         index1 = F)
  R_lower + t(R_lower) 
}

# this can be done outside of the loop
# create a tibble with all of the bin information for the reference LD matrices
basedir = "/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/"
blocks = tibble( fn = list.files(basedir, glob2rx("*.npz")), size = file.info(dir(path = basedir, full.names = TRUE, pattern = "*.npz"))$size) %>%
  mutate(stub = substr(fn, 1, nchar(fn) - 4)) %>%
  separate(stub, c("chr", "start", "end"), convert = T, remove = F)

# read in the sQTL summary statistics
stats = read_tsv(paste0("~/ipsc_sqtl/Data/full_qtl_results_", args[1], ".Pval-rescaled.txt.gz"), show_col_types = FALSE)
stats_grouped = stats %>% group_by(gene) %>% summarize() %>% ungroup()

# foreach loop will iterate over each junction from the summary statistics, find the appropriate reference LD matrix, and run SuSie
# the command line specifies how many parallel jobs you want here- because this uses a lot of memory, I found that I could have a maximum of 5 parallel jobs when I requested 240G of memory
registerDoMC(args[2])

# if command line specifies TRUE, then errors within the loop, execution will be stopped. I recommend this option when running for the first time
# chr6 throws errors because of the MHC region which has extremely long-range LD regions- for this, I would use FALSE which will still continue if there is an error but will not return the output for that iteration
if(as.logical(args[3]) == TRUE) {
  error <- "stop"
  } else{
    error <- "remove"
}

fine_map = foreach(id = 1:nrow(stats_grouped), .combine = bind_rows, .errorhandling = error) %dopar% {
  bin = as.character(stats_grouped[id, ])
  junc_data = stats %>% filter(gene == bin)
  
  chrom = junc_data$chr[1]  
  snp_range_start = min(junc_data$snp_pos)
  snp_range_end = max(junc_data$snp_pos)
  
  # find the proper npz file for the current junction-snp pairs
  theblock = blocks %>% filter(chr == paste0("chr", chrom), start < snp_range_start, snp_range_end < end ) %>% arrange(size) %>% head(1)
  matrix_size = round(theblock$size/1000000000, 2)
  
  # if statement checks that there is a corresponding reference LD matrix for the current region
  if(nrow(theblock > 0)) {
    
    # this file has the corresponding snp info for the matrix, add an index column here for merging later
    snps = read_tsv(paste0(basedir, theblock$stub, ".gz"), show_col_types = FALSE) %>% 
      mutate(ind = 1:n())
    
    # inner join, the snp info file with sQTL file- subsets to only snps that we have sQTL info for
    # important here to merge with allele info since some SNPs will have multiple alt alleles
    # there are some sQTLs that don't have a corresponding SNP id in the matrix- so these are not included- in this first example, there are 8
    joined = snps %>% inner_join(junc_data, by = c(chromosome = "chr", position = "snp_pos", allele1 = "ref", allele2 = "alt"))
    
    # this file has the lower matrix
    cat("Loading LD lower matrix for iteration", id, "of size", matrix_size, "\n")
    R = read_ld(paste0(basedir, theblock$stub, ".npz") )
    
    cat("Finished loading matrix for iteration", id, "\n")
    
    # subset the LD matrix for only snps with sQTL info
    R_sub = as.matrix(R[joined$ind,joined$ind])
    remove(R)
    
    # calculate z scores from summary stats
    z_scores = joined$beta / joined$se
    
    # run SuSie using summary statistics, can change L here
    # add pip score to output file
    cat("Running SuSie for iteration", id, "\n")
    fitted_rss <- susie_rss(z_scores, R_sub, L = 10)
    joined$pip = fitted_rss$pip
    cat("Finished loop for iteration", id, "\n")
    
    remove(R_sub)
    remove(fitted_rss)
    remove(z_scores)
    return(joined)
    
  } else{
    cat("No matrix for", bin, "\n")
  }
  
}

warnings()

output = paste0("~/ipsc_sqtl/test/ipsc_sqtl_finemapped_L10_chr", args[1], ".txt.gz")
write_tsv(fine_map, output, col_names = TRUE)
