
require(reticulate)
library(Matrix)

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

basedir = "/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/"
R = read_ld(paste0(basedir, "chr11_53000001_56000001.npz2") )
R[1:5,1:5]