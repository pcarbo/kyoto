# *** DESCRIPTION OF SCRIPT GOES HERE ***

# SCRIPT PARAMETERS
# -----------------
phenotype  <- "freezetocue" # Map QTLs for this phenotype.
generation <- "F2"          # Map QTLs in mice from this generation.

# Use these covariates in the QTL mapping.
covariates <- NULL # c("sex","age","albino","agouti")

# Initialize the random number generator.
set.seed(7)

# Load function and symbol definitions from libraries and source files.
library(lattice)
library(varbvs)
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
map     <- transform(map[markers,],chr = droplevels(chr))
geno    <- geno[,markers]

# GET F2 OR F34 CROSS
# -------------------
# Keep only the mice from the selected generation.
rows  <- which(pheno$generation == generation)
pheno <- pheno[rows,]
geno  <- geno[rows,]

# Get the markers genotyped in these mice.
markers <- which(!all.missing.col(geno))
map     <- map[markers,]
geno    <- geno[,markers]

# Keep only samples for which we have all observations for the
# phenotype and covariates.
cols  <- c(phenotype,covariates)
rows  <- which(none.missing.row(pheno[cols]))
pheno <- pheno[rows,]
geno  <- geno[rows,]

# COMPUTE GENOTYPE PROBABILITIES
# ------------------------------
# Compute the conditional genotype probabilities using QTLRel. To
# accomplish this, we need to replace the genotypes with allele counts
# (AA, AB, BB become 1, 2, 3, respectively), and we replace any 
# missing values with zeros.
cat("Calculating probabilities of missing genotypes.\n")
gp <- genoProb(zero.na(genotypes2counts(geno)),map,step = Inf,
               gr = as.integer(substr(generation,2,3)),
               method = "Haldane")

# varbvsoptimize <- function (X, y, sigma, sa, logodds, alpha0 = NULL,
#                             mu0 = NULL, verbose = TRUE) {

