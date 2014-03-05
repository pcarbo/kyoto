# *** DESCRIPTION OF SCRIPT GOES HERE ***

# SCRIPT PARAMETERS
# -----------------
phenotype  <- "coatcolor" # Map QTLs for this phenotype.
generation <- "F2"        # Map QTLs in mice from this generation.

# Use these covariates in the QTL mapping.
covariates <- NULL # c("sex","age","albino","agouti")

# Candidate values of the variance of the residual (sigma), the prior
# variance of the regression coefficients (sa), and the prior
# log-odds of inclusion (log10odds)
sigma     <- seq(1,10,2)
sa        <- seq(0.1,0.4,0.05)
log10odds <- seq(-2,-1,0.25)

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

# Get the mean genotypes; that is, the mean counts of the B allele.
p1   <- gp$pr[,2,]
p2   <- gp$pr[,3,]
geno <- p1 + 2*p2

# ANALYSIS USING MULTI-MARKER REGRESSION
# --------------------------------------
# Generate all combinations of the candidate values for the hyperparameters.   
grid        <- grid3d(sigma,sa,log10odds)
names(grid) <- c("sigma","sa","log10odds")

# Get the number of markers (p) and the number of combinations of the
# hyperparameters (ns).
p  <- ncol(geno)
ns <- length(grid$sigma)

# Initialize storage for the marginal log-likelihoods (lnZ),
# variational estimates of the posterior inclusion probabilities
# (alpha), and variational estimates of the posterior mean
# coefficients (mu).
lnZ   <- array(dim = dim(grid$sigma))
alpha <- matrix(NA,ncol(geno),ns)
mu    <- matrix(NA,ncol(geno),ns)

# Choose a random initialization for the variational parameters.
alpha0 <- runif(p)
alpha0 <- alpha0 / sum(alpha0)
mu0    <- rnorm(p)

# Repeat for each combination of the hyperparameters.
cat("Finding best variational approximation to marginal log-likelihood ")
cat("for",ns,"combinations of hyperparameters.\n")
for (i in 1:ns) {

  # Run the coordinate ascent algorithm.
  out <- varbvsoptimize(geno,pheno[,phenotype],sigma[i],sa[i],
                        log(10)*log10odds[i],alpha0,mu0,verbose = TRUE)
  lnZ[i]    <- out$lnZ
  alpha[,i] <- out$alpha
  mu[,i]    <- out$mu
}
