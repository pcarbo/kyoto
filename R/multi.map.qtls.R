# In this script, we explore an alternative approach to mapping QTLs
# based on simultaneously considering all SNPs as QTLs. The
# multi-marker mapping method is implemented in the varbvs package.

# SCRIPT PARAMETERS
# -----------------
phenotype  <- "freezetocue"  # Map QTLs for this phenotype.
generation <- "F34"          # Map QTLs in mice from this generation.

# Use these covariates in the QTL mapping.
covariates <- c("sex","age","albino","agouti")

# Candidate values of the variance of the residual (sigma), the prior
# variance of the regression coefficients (sa), and the prior
# log-odds of inclusion (log10odds)
sigma     <- seq(0.1,0.4,0.1)
sa        <- c(0.01,0.02,0.05,0.1,0.2)
log10odds <- seq(-3,-1,0.25)

# Initialize the random number generator.
set.seed(7)

# Load function and symbol definitions from libraries and source files.
library(lattice)
library(varbvs)
library(qtl)
capture.output(library(QTLRel))
source("misc.R")
source("read.data.R")
source("data.manip.R")
source("mapping.tools.R")

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
p1 <- gp$pr[,2,]
p2 <- gp$pr[,3,]
X  <- p1 + 2*p2

# ANALYSIS USING MULTI-MARKER REGRESSION
# --------------------------------------
# Generate all combinations of the candidate values for the hyperparameters.   
grid        <- grid3d(sigma,sa,log10odds)
names(grid) <- c("sigma","sa","log10odds")

# Get the number of markers (p) and the number of combinations of the
# hyperparameters (ns).
p  <- ncol(X)
ns <- length(grid$sigma)

# Initialize storage for the marginal log-likelihoods (lnZ),
# variational estimates of the posterior inclusion probabilities
# (alpha), and variational estimates of the posterior mean
# coefficients (mu).
lnZ   <- array(dim = dim(grid$sigma))
alpha <- matrix(NA,p,ns)
mu    <- matrix(NA,p,ns)

# Choose a random initialization for the variational parameters.
alpha0 <- runif(p)
alpha0 <- alpha0 / sum(alpha0)
mu0    <- rnorm(p)

# Compute the residuals of the phenotype given the covariates.
if (is.null(covariates)) {
  y <- pheno[,phenotype]
} else {
  f <- formula(paste(phenotype,"~",paste(covariates,collapse=" + ")))
  y <- resid(lm(f,pheno))
}

# Center the columns of X, and take into account an intercept by
# centering the outcomes Y to have mean zero.
X <- center.columns(X)
y <- y - mean(y)

# Repeat for each combination of the hyperparameters.
cat("Finding best variational approximation to marginal log-likelihood\n")
cat("for",ns,"combinations of hyperparameters.\n")
for (i in 1:ns) {
  with(grid,cat(sprintf("(%d) sigma = %0.1f, sa = %0.3f, logodds = %0.2f\n",
                        i,sigma[i],sa[i],log10odds[i])))

  # Run the coordinate ascent algorithm.
  out <- with(grid,varbvsoptimize(X,y,sigma[i],sa[i],log(10)*log10odds[i],
                                  alpha0,mu0,verbose = FALSE))
  lnZ[i]    <- out$lnZ
  alpha[,i] <- out$alpha
  mu[,i]    <- out$mu
}

# Get the posterior inclusion probabilities corresponding to the
# combination of hyperparameters with the maximum likelihood.
i   <- which.max(lnZ)
PIP <- alpha[,i]

# Create a data frame with the genome-wide scan.
gwscan        <- empty.scanone(map)
gwscan$lod    <- PIP
names(gwscan) <- c("chr","pos","PIP")

# PLOT GENOME-WIDE SCAN FROM varbvs
# ---------------------------------
trellis.device(width = 7,height = 1.75,
               title = "Multi-marker QTL mapping results")
par(ps = 9,font.lab = 1,font.main = 1,cex.main = 1,
    mai = c(0.5,0.8,0.25,0.25),tck = -0.1)

# Plot the posterior inclusion probabilities.
plot(gwscan,incl.markers = FALSE,lwd = 2,bandcol = "powderblue",
     col = "darkorange",gap = 0,xlab = "chromosome",ylab = "probability",
     main = paste0(phenotype,", ",generation," cross"))

# Add the cutoff at a posterior probability of 0.9.
abline(0.9,0,col = "black",lty = "dotted")
