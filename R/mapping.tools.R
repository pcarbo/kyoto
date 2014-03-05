# This file contains some functions to analyze the QTL experiment data.

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# This function creates a "scanone" object that will be used to store
# QTL mapping results based on the SNP data provided as input. (For
# further info, see the 'scanone' function in the 'qtl' package.) 
empty.scanone <- function (map) {
  d            <- map[c("chr","dist")]
  row.names(d) <- map[["snp"]]
  names(d)     <- c("chr","pos")
  class(d)     <- c("scanone","data.frame")
  return(d) 
}

# ----------------------------------------------------------------------
# Use all the markers to estimate the n x n pairwise relatedness
# matrix, which I define to be 2 times the matrix of kinship
# coefficients (for the definition of kinship coefficients, see Lange,
# 2002). To allow for uncertainty in the genotypes at the markers, I
# compute the expected value of this matrix using the genotype
# probabilities provided as output from the 'genoProb' function in
# QTLRel.
rr.matrix <- function (gp) { 

  # Get the number of samples (n) and the number of markers (p).
  d <- dim(gp$pr)
  n <- d[1]
  p <- d[3]
  
  # Get the genotype probabilities.
  pAA <- gp$pr[,1,]
  pAB <- gp$pr[,2,]
  pBB <- gp$pr[,3,]

  # Get the probability of the genotype AB averaged over all the markers.
  mAB <- matrix(rep(rowMeans(pAB),times = n),n,n)

  # Return the expected value of the pairwise relatedness matrix.
  return(2*(matrix.square(pAA)/p + matrix.square(pBB)/p) 
         + mAB + t(mAB) - matrix.square(pAB)/p)
}
