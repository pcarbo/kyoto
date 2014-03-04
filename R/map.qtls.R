# This script maps QTLs across the genome in a single filial
# generation of an advanced intercross line (AIL). We quantify support
# for a QTL at each of the genotyped SNPs. We use two different
# methods to assess support for an association between genotype and
# phenotype:
#
#   qtl, which does not account for varying amounts of genetic
#   sharing;
#
#   QTLRel, which attempts to account for confounding due to unequal
#   relatedness by using the marker data to estimate relatedness
#   between all pairs of mice in the sample.
#
# We also use two different approaches are to assess thresholds for
# significance:
#
#   The first approach applies a permutation-based test to estimate
#   the null distribution of LOD scores, in which unequal relatedness
#   of the animals is not taken into account;
#
#   The second approach applies the method described in Abney, Ober
#   and McPeek (American Journal of Human Genetics, 2002), in which
#   the phenotype measurements are permuted according to a specified
#   covariance matrix (which will be calculated using the markers).
#

# SCRIPT PARAMETERS
# -----------------
generation   <- "F2"       # Map QTLs in mice from this generation.
qtl.method   <- "hk"       # Which QTL mapping method to use in qtl.
map.function <- "Haldane"  # Map function to use for interval mapping.
num.perm     <- 1000       # Number of replicates for permutation test.

# Map QTLs for this phenotype
phenotypes <- "freezetocue"

# Use these covariates in the QTL mapping.
covariates <- c("sex","age","albino","agouti")

# Initialize the random number generator.
set.seed(7)

# Load function and symbol definitions from libraries and source files.
library(qtl)
capture.output(library(QTLRel))
source("misc.R")
source("read.data.R")
source("data.manip.R")

# LOAD DATA
# ---------
# Load the phenotype, phenotype and marker data for all the samples.
cat("Loading phenotype, genotype and marker data.\n");
pheno <- read.pheno("../data/pheno.csv")
map   <- read.map("../data/map.csv")
geno  <- read.geno("../data/geno.csv")

# Drop the X chromosome from the analysis.
markers <- which(map$chr != "X")
geno    <- geno[,markers]
map     <- transform(map[markers,],chr = droplevels(chr))

# Remove 3 mice from the F2 cohort because a large proportion of their
# genotypes are not available.
rows  <- which(!is.element(pheno$id,c("648A","738A","811A")))
pheno <- pheno[rows,]
geno  <- geno[rows,]

# GET F2 OR F34 CROSS
# -------------------
rows  <- which(pheno$generation == generation)
pheno <- pheno[rows,]
geno  <- geno[rows,]

stop()

# Get the markers genotyped in the F2 cross. Since I'm only retaining
# markers that are genotyped in the F34 cross, any marker genotyped in
# the F2 cross is also genotyped in the F34 cross.
F2.markers <- which(!all.missing.col(geno[F2.rows,]))

# COMPUTE GENOTYPE PROBABILITIES
# ------------------------------
# Compute the conditional genotype probabilities using QTLRel. To
# accomplish this, we need to replace the genotypes with allele counts
# (AA, AB, BB become 1, 2, 3, respectively), and we replace any 
# missing values with zeros.
cat("Calculating probabilities of missing genotypes.\n")
gp <- genoProb(zero.na(genotypes2counts(G)),map,step = Inf,
               method = "Haldane",gr = generation)

# Initialize storage for QTL mapping results (gwscan), parameter
# estimates of variance components from QTLRel analysis (vcparams),
# parameter estimates of additive QTL effects for individual markers
# (additive) and dominance QTL effects (dominance), and permutation
# tests (perms).
gwscan    <- list(F2.qtl   = NULL,
                  F2.rel   = empty.scanone(map[F2.markers,]),
                  F34      = empty.scanone(map),
                  combined = empty.scanone(map))
r         <- list(F2       = empty.scanone(map[F2.markers,]),
                  F34      = empty.scanone(map),
                  combined = empty.scanone(map))
additive  <- r
dominance <- r
pve       <- r
vcparams  <- list(F2 = NULL,F34 = NULL,combined = NULL)
perms     <- list(F2 = NULL,F34 = NULL,combined = NULL)

# ANALYSIS WITH F2 CROSS USING QTL
# --------------------------------
cat("QTL MAPPING WITH F2 CROSS USING qtl\n")
out <- map.cross.qtl(pheno,geno,phenotypes,covariates,gp,
                     F2.rows,F2.markers,qtl.method,num.perm,
                     verbose = TRUE)
gwscan$F2.qtl <- out$gwscan
perms$F2      <- out$perms

# PERMUTATION TESTS FOR F34 AND COMBINED CROSSES
# ----------------------------------------------
# Perform permutation tests to calculate thresholds for significance
# for F34 cross, ignoring relatedness between mice.
cat("CALCULATING SIGNIFICANCE THRESHOLDS FOR F34 CROSS USING qtl\n")
perms$F34 <- map.cross.qtl(pheno,geno,phenotypes,covariates,gp,
                           F34.rows,NULL,qtl.method,num.perm,
                           verbose = FALSE)$perms

# Perform a permutation test to calculate thresholds for significance
# for combined F2 and F34 cohort, ignoring relatedness between mice.
cat("CALCULATING SIGNIFICANCE THRESHOLDS FOR COMBINED COHORT USING qtl\n")
perms$combined <- map.cross.qtl(pheno,geno,phenotypes,covariates,gp,
                                qtl.method = qtl.method,num.perm = num.perm,
                                verbose = FALSE)$perms

# Repeat for each phenotype.
for (phenotype in phenotypes) {
  cat("PHENOTYPE =",toupper(phenotype),"\n")
  cols <- c(phenotype,covariates)
  
  # ANALYSIS WITH F2 CROSS USING QTLRel
  # -----------------------------------
  # Only analyze samples (i.e. rows of the genotype and phenotype
  # matrices) for which the phenotype and all the covariates are
  # observed.
  cat("(a) QTL mapping with F2 cross using QTLRel\n")
  F2.rows <- which(pheno$generation == "F2" &
                   none.missing.row(pheno[cols]))
  if (relatedness == "pedigree")
    out <- map.cross.ped(pheno,geno,map,phenotype,covariates,F2.gm,gp,
                         F2.rows,F2.markers,verbose = TRUE)
  else if (relatedness == "markers")
    out <- map.cross.rr(pheno,geno,map,phenotype,covariates,gp, 
                        F2.rows,F2.markers,verbose = TRUE)
  
  # Get the parameter estimates and QTL mapping results.
  gwscan$F2.rel[[phenotype]] <- out$gwscan$lod
  additive$F2[[phenotype]]   <- out$gwscan$additive
  dominance$F2[[phenotype]]  <- out$gwscan$dominance
  pve$F2[[phenotype]]        <- out$gwscan$pve
  vcparams$F2                <- rbind(vcparams$F2,out$vcparams)

  # ANALYSIS WITH F34 CROSS USING QTLRel
  # ------------------------------------
  # Only analyze samples (i.e. rows of the genotype and phenotype
  # matrices) for which the phenotype and all the covariates are
  # observed.
  cat("(b) QTL mapping with F34 cross using QTLRel\n")
  F34.rows <- which(pheno$generation == "F34" &
                    none.missing.row(pheno[cols]))
  if (relatedness == "pedigree")
    out <- map.cross.ped(pheno,geno,map,phenotype,covariates,F34.gm,gp,
                         rows = F34.rows,verbose = TRUE)
  else if (relatedness == "markers")
    out <- map.cross.rr(pheno,geno,map,phenotype,covariates,gp,
                        rows = F34.rows,verbose = TRUE)

  # Get the parameter estimates and QTL mapping results.
  gwscan$F34[[phenotype]]    <- out$gwscan$lod
  additive$F34[[phenotype]]  <- out$gwscan$additive
  dominance$F34[[phenotype]] <- out$gwscan$dominance
  pve$F34[[phenotype]]       <- out$gwscan$pve
  vcparams$F34               <- rbind(vcparams$F34,out$vcparams)

  # ANALYSIS WITH F2 AND F34 SAMPLES USING QTLRel
  # ---------------------------------------------
  cat("(c) QTL mapping with combined F2 and F34 cohorts using QTLRel\n")
  rows <- c(F2.rows,F34.rows)
  if (relatedness == "pedigree")
    out <- map.cross.ped(pheno,geno,map,phenotype,covariates,
                         gm,gp,rows = rows,verbose = TRUE)
  else if (relatedness == "markers")
    out <- map.cross.rr(pheno,geno,map,phenotype,covariates,gp, 
                        rows = rows,verbose = TRUE)

  # Get the parameter estimates and QTL mapping results.
  gwscan$combined[[phenotype]]    <- out$gwscan$lod
  additive$combined[[phenotype]]  <- out$gwscan$additive
  dominance$combined[[phenotype]] <- out$gwscan$dominance
  pve$combined[[phenotype]]       <- out$gwscan$pve
  vcparams$combined               <- rbind(vcparams$combined,out$vcparams)
}

# Add labels to the parameter estimates.
vcparams <- within(vcparams,{
  F2       <- data.frame(F2,check.names = FALSE,row.names = phenotypes)
  F34      <- data.frame(F34,check.names = FALSE,row.names = phenotypes)
  combined <- data.frame(combined,check.names = FALSE,row.names = phenotypes)
})
