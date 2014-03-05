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
phenotype  <- "albino"      # Map QTLs for this phenotype.
generation <- "F2"          # Map QTLs in mice from this generation.
num.perm   <- 100           # Number of replicates for permutation test.
threshold  <- 0.05          # Significance threshold ("alpha").
ymax       <- 9             # Height of vertical axis.

# Use these covariates in the QTL mapping.
covariates <- NULL # c("sex","age","albino","agouti")

# Initialize the random number generator.
set.seed(7)

# Load function and symbol definitions from libraries and source files.
library(lattice)
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

# ANALYSIS USING QTL
# ------------------
# Map QTLs using a simple linear regression approach that does not
# correct for possible confounding to due relatedness.
cat("Mapping QTLs for",phenotype,"in",nrow(pheno),"mice at",nrow(map),
    "candidate SNPs,\n")
if (length(covariates) > 0) {
  cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
} else {
  cat("with no covariates included.\n")
}
cat("Mapping QTLs using qtl.\n")

# Convert the experimental cross data to the format used by qtl, then
# run the genome-wide scan using the qtl function "scanone".
cross    <- rel2qtl(pheno,geno,map)
cross    <- rel2qtl.genoprob(cross,gp)
cov.data <- NULL
if (!is.null(covariates))
  cov.data <- cross$pheno[covariates]
suppressWarnings(
  gwscan.qtl <- scanone(cross,pheno.col = phenotype,addcovar = cov.data,
                        model = "normal",method = "em",use = "all.obs"))

# Estimating the null distribution of the LOD scores using qtl.
cat("Estimating null distribution of LOD scores using qtl.\n")
suppressWarnings(
  perms.qtl <- scanone(cross,pheno.col = phenotype,addcovar = cov.data,
                       model = "normal",method = "em",use = "all.obs",
                       n.perm = num.perm,verbose = FALSE))

# ESTIMATE GENETIC SHARING USING MARKERS
# --------------------------------------
# Compute the (expected) relatedness matrix using all available markers.
R <- rr.matrix(gp);
trellis.device(width = 4,height = 2,title = "Relatedness estimates")
trellis.par.set(list(fontsize = list(text = 9),
                     par.main.text = list(cex = 1)))

# Plot the distribution of the marker-based estimates of pairwise
# relatedness when the pairs correspond to the same individual.
print(histogram(diag(R),breaks = seq(1.2,1.8,0.025),col = "darkorange",
                border = "darkorange",xlab = "relatedess coef.",
                ylab = "% of entries",main = "diagonal entries",
                scales = list(tck = 0.75,x = list(at = seq(1.2,1.8,0.15)))),
      split = c(1,1,2,1),
      more = TRUE)

# Plot the distribution of the marker-based estimates of pairwise
# relatedness when the pairs correspond to different individuals.
print(histogram(offdiag(R),breaks = seq(0.6,1.4,0.025),col = "dodgerblue",
                border = "dodgerblue",xlab = "relatedness coef.",
                ylab = "% of entries",main = "off-diagonal entries",
                scales = list(tck = 0.75,x = list(at = seq(0.6,1.4,0.2)))),
      split = c(2,1,2,1))
      
# ANALYSIS USING QTLRel
# ---------------------
# Map QTLs for all markers on a single chromosome using 'scanOne' from
# the QTLRel library, in which pairwise relatedness is estimated
# using all markers except the markers on the same chromosome as the
# one being analyzed.
cat("Mapping QTLs using QTLRel.\n")
chromosomes <- levels(map$chr)

# Map QTLs separately for each chromosome.
gwscan <- list()
for (chr in chromosomes) {
  
  # Get the markers on the chromosome.
  markers <- which(map$chr == chr)
  cat("  * Mapping ",length(markers)," markers on chromosome ",chr,".\n",
      sep="")

  # Compute the (expected) relatedness matrix using all markers *except* 
  # the markers on the current chromosome.
  R <- rr.matrix(subset.genoprob(gp,which(map$chr != chr)))
  dimnames(R) <- list(pheno$id,pheno$id)

  # Use the relatedness matrix estimated from the marker data to
  # estimate the variance components.
  if (is.null(covariates)) {
    r <- estVC(pheno[,phenotype],
               v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                        MH = NULL,EE = diag(nrow(pheno))))
  } else {
    r <- estVC(pheno[,phenotype],pheno[,covariates],
               v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                        MH = NULL,EE = diag(nrow(pheno))))
  }

  # Once we have the variance components estimated, build a matrix
  # from the variance components. Note that the "AA" component here is
  # the n x n relatedness matrix estimated from the marker data (R),
  # and the "EE" component is the n x n identity matrix.
  vc <- r$par["AA"] * r$v[["AA"]] +
        r$par["EE"] * r$v[["EE"]]
  
  # Compute LOD scores for all markers on the chromosome.
  if (is.null(covariates)) {
    out <- scanOne(pheno[,phenotype],gdat = geno[,markers],
                   prdat = subset.genoprob(gp,markers),
                   vc = vc,test = "None")
  } else {
    out <- scanOne(pheno[,phenotype],pheno[,covariates],geno[,markers],
                   subset.genoprob(gp,markers),vc,test = "None")
  }
  gwscan[[chr]]     <- empty.scanone(map[markers,])
  gwscan[[chr]]$lod <- out$p/(2*log(10))
}

# Merge the QTLRel mapping results.
gwscan.rel           <- do.call(rbind,gwscan)
rownames(gwscan.rel) <- do.call(c,lapply(gwscan,rownames))

# PLOT RESULTS FROM qtl AND QTLRel
# --------------------------------
trellis.device(width = 7,height = 1.75,title = "QTL mapping results")
par(ps = 9,font.lab = 1,font.main = 1,cex.main = 1,
    mai = c(0.5,0.5,0.25,0.25),tck = -0.1)

# Plot the QTL mapping results from qtl.
plot(gwscan.qtl,incl.markers = FALSE,lwd = 4,bandcol = "powderblue",
     col = "dodgerblue",gap = 0,xlab = "",ylab = "LOD",
     main = paste0(phenotype,", ",generation," cross"))

# Plot the QTL mapping results from QTLRel.
plot(gwscan.rel,incl.markers = FALSE,lwd = 2,bandcol = "powderblue",
     col = "darkblue",gap = 0,add = TRUE)

# Add the significance threshold to the plot.
add.threshold(gwscan.qtl,perms = perms.qtl,alpha = threshold,gap = 0,
              col = "darkorange",lty = "dotted")
