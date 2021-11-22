#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load required packages
require(reticulate)
library(Matrix)
require(doMC)
library(foreach)
library(readr)
require(tidyverse)
library(susieR)
#library(pryr)

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
# create a dataframe with all of the bin information for the LD matrices
basedir = "/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/"
blocks = tibble( fn = list.files(basedir, glob2rx("*.gz")) ) %>% 
  mutate(stub = substr(fn, 1, nchar(fn) - 3)) %>%
  separate(stub, c("chr", "start", "end"), convert = T, remove = F)

# read in the sQTL mapped file and find the corresponding LD matrix for each junction-snp pair
registerDoMC(2)

stats = read_tsv(paste0("~/ipsc_sqtl/SplicingLevel/full_qtl_results_", args[1], ".Pval-rescaled.txt.gz"), show_col_types = FALSE)
stats_grouped = stats %>% group_by(gene) %>% summarize() %>% ungroup()

fine_map = foreach(id = 1:nrow(stats_grouped), .combine = bind_rows) %dopar% {
  #gcinfo(TRUE)
  bin = as.character(stats_grouped[id, ])
  junc_data = stats %>% filter(gene == bin)
  chrom = junc_data$chr[1]  
  snp_range_start = min(junc_data$snp_pos)
  snp_range_end = max(junc_data$snp_pos)
  
  # find the proper npz file for the current junction-snp pairs
  theblock = blocks %>% filter(chr == paste0("chr", chrom), start < snp_range_start, snp_range_end < end ) %>% head(1)
  
  # this file has the corresponding snp info for the matrix, add an index column here for merging later
  snps = read_tsv(paste0(basedir, theblock$fn), show_col_types = FALSE) %>% 
    mutate(ind = 1:n())
  #snps %>% group_by(rsid) %>% summarize(n = n()) %>% filter(n > 1)
  #snps %>% filter(rsid == "rs112108829")
  
  # inner join, the snp info file with sQTL file- subsets to only snps that we have sQTL info for
  # important here to merge with allele info since some SNPs will have multiple alt alleles
  # there are some sQTLs that don't have a corresponding SNP id in the matrix- so these are not included- in this first example, there are 8
  joined = snps %>% inner_join(junc_data, by = c(chromosome = "chr", position = "snp_pos", allele1 = "ref", allele2 = "alt"))
  
  # this file has the lower matrix
  cat("Loading LD lower matrix for iteration", id, "\n")
  R = read_ld(paste0(basedir, theblock$stub, ".npz") )
  
  # subset the LD matrix for only snps with sQTL info
  R_sub = as.matrix(R[joined$ind,joined$ind])
  remove(R)
  
  # calculate z scores
  z_scores = joined$beta / joined$se
  
  # run SuSie using summary statistics, can change L here
  cat("Running SuSie for iteration", id, "\n")
  fitted_rss <- susie_rss(z_scores, R_sub, L = 10)
  #susie_plot(fitted_rss, y="PIP")
  joined$pip = fitted_rss$pip
  cat("Finished loop for iteration", id, "\n")
  
  remove(R_sub)
  remove(fitted_rss)
  remove(z_scores)
  
  return(joined)
}

output = paste0("~/ipsc_sqtl/fine_mapped/ipsc_sqtl_finemapped_L10_chr", args[1], ".txt.gz")
write_tsv(fine_map, output, col_names = TRUE)
