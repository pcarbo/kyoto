# Some functions for manipulating QTL experiment data.
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Convert the genotypes to a matrix of allele counts.
genotypes2counts <- function (geno)
  as.matrix(data.frame(lapply(geno,factor2integer),
                       row.names = rownames(geno),
                       check.names = FALSE))

# ----------------------------------------------------------------------
# Convert the QTL experiment data from the format used in the QTLRel
# library to the format used in qtl library. The return value is a
# 'cross' object that keeps track of all the data in a single QTL
# experiment; for more details, see the help for function read.cross
# in the qtl library. NOTE: this function currently does not properly
# handle the X chromosome.
#
# Here I assume that the alleles are A and B and the genotypes are
# labeled as AA, AB and BB. qtl requires a paternal grandmother
# ("pgm") phenotype, so a column in the table for this phenotype is
# included if it is missing, and the entries of this column are set to
# zero.
rel2qtl <- function (pheno, geno, map) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)
  
  # Convert the genotypes to a matrix of integers.
  geno <- genotypes2counts(geno)
  
  # Add the paternal grandmother ("pgm") phenotype if it does not
  # already exist.
  if (!is.list.element("pgm",pheno))
    pheno <- cbind(pheno,pgm = 0)
  
  # Initialize the cross object, setting the genotypes to an empty list.
  cross <- list(geno = list(),pheno = pheno)
  class(cross) <- c("f2","cross")

  # Set the alleles.
  attributes(cross)$alleles <- c("A","B")
  
  # Split the markers by chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers <- which(map$chr == i)
    
    # Get the genotype and map data corresponding to these markers.
    # Note that I need to coerce the genotype data to be a matrix in
    # the exceptional case when there is only a single marker
    # genotyped on the chromosome.
    m                <- map[markers,]
    d                <- list(data = as.matrix(geno[,markers]),map = m$dist)
    colnames(d$data) <- m$snp
    names(d$map)     <- m$snp

    # Store the genotypes for markers on the chromosome.
    class(d) <- "A"
    cross$geno[[i]] <- d
  }
  
  return(cross)
}

# ----------------------------------------------------------------------
# Convert the genotype probabilities from the format used by QTLRel to
# the format used by qtl, and store the genotype probabilities in the
# cross object. NOTE: I currently do not properly handle the X
# chromosome.
rel2qtl.genoprob <- function (cross, gp) {

  # Transpose the array.
  prob <- aperm(gp$pr,c(1,3,2))

  # Get the set of chromosomes.
  chromosomes <- unique(gp$chr)

  # Split the genotype probabilities by chromosome.
  for (i in chromosomes) {

    # Get the data for all markers on this chromosome.
    geno <- cross$geno[[i]]
    
    # Get the markers on the chromosome.
    markers <- which(gp$chr == i)

    # Get the genotype probabilities for these markers.
    geno$prob                 <- prob[,markers,]
    attributes(geno$prob)$map <- geno$map
    dimnames(geno$prob)       <- list(NULL,
                                      names(geno$map),
                                      c("AA","AB","BB"))

    # Store the new genotype data.
    cross$geno[[i]] <- geno
  }
  
  return(cross)
}

# ----------------------------------------------------------------------
# Get the genotype probabilities for the specified markers.
subset.genoprob <- function (gp, markers) {
  class(gp) <- c("Pr","list")
  return(within(gp,{
    pr   <- pr[,,markers]
    chr  <- chr[markers]
    dist <- dist[markers]
    snp  <- snp[markers]
  }))
}
