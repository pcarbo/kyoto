# This file contains several functions for reading the QTL experiment
# data from files in comma-delimited ("csv") format.
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns a data frame containing the phenotype data stored in a CSV file.
read.pheno <- function (file)
  read.csv(file,comment.char="#",as.is = "id",check.names = FALSE)
