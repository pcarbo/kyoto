# This file contains several functions for reading the QTL experiment
# data from files in comma-delimited ("csv") format.
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns a data frame containing the phenotype data stored in a CSV
# file. I transform pretrainfreeze, freezetocontext and freezetocue
# phenotypes to the log-odds scale using the base-10 logit function.
# Further, I transform sex to a binary trait, subtract the mean from
# age, and create traits albino and agouti from the coat colour data.
read.pheno <- function (file) {
  pheno       <- read.csv(file,comment.char="#",as.is = "id",
                          check.names = FALSE)
  cols        <- c("pretrainfreeze","freezetocontext","freezetocue")
  f           <- function (x) logit10(project.onto.interval(x,0.01,0.99))
  pheno[cols] <- lapply(pheno[cols],f)
  return(transform(pheno,
                   age    = age - mean(age),
                   sex    = as.integer(sex       == "M"),
                   albino = as.integer(coatcolor == "W"),
                   agouti = as.integer(coatcolor == "A")))
}

# ----------------------------------------------------------------------
# Returns a data frame containing the genotype data stored in a CSV
# file. I remove the first two columns of the genotype data frame
# containing the ID and generation, and set the row names for the
# genotype data to the mouse IDs.
read.geno <- function (file) {
  geno <- read.csv(file,comment.char="#",as.is = "id",check.names = FALSE)
  rownames(geno) <- geno$id
  return(geno[-(1:2)])
}

# ----------------------------------------------------------------------
# Returns a data frame containing the marker data stored in a CSV
# file. Here I convert the chromosomes to factors manually (as
# opposed to letting the read.csv function handle this) to make sure
# that the chromosomes are ordered properly. I also: adjust the map
# positions (i.e. genetic distances) slightly so that no two markers
# have the same position.
read.map <- function (file) {
  map <- read.csv(file,comment.char="#",check.names = FALSE,
                  stringsAsFactors = FALSE)
  map <- transform(map,chr = factor(chr,c(1:19,"X")))
  return(jitter.gendist(map,1e-6))
}

# ----------------------------------------------------------------------
# Adjust the map positions (i.e. genetic distances) slightly (by the
# amount "j") so that no two markers have the same position.
jitter.gendist <- function (map, j) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)

  # Repeat for each chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers  <- which(map$chr == i)
    n        <- length(markers)

    # Adjust the genetic distances slightly according to scalar j.
    map[markers,"dist"] <- map[markers,"dist"] + cumsum(rep(j,times = n))
  }

  return(map)
}
  
