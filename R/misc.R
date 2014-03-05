# This file defines a number of functions that do not fit under any
# one category.
#
# SYMBOL DEFINITIONS
# ----------------------------------------------------------------------
# Shorthand for machine epsilon.
eps <- .Machine$double.eps

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns true if x is an element of the list.
is.list.element <- function (x, lst)
  is.element(x,names(lst))

# ----------------------------------------------------------------------
# The logit function.
logit10 <- function (x)
  log10((x + eps)/(1 - x + eps))

# ----------------------------------------------------------------------
# Projects points to the interval [a,b].
project.onto.interval <- function (x, a, b)
  pmin(b,pmax(a,x))

# ----------------------------------------------------------------------
# Convert a factor to a vector of integers.
factor2integer <- function (x)
  match(x,levels(x))


# ----------------------------------------------------------------------
# Set the missing entries (NA) of a matrix to zero.
zero.na <- function (A) {
  A[is.na(A)] <- 0
  return(A)
}

# ----------------------------------------------------------------------
# Returns a logical vector such that an entry is TRUE if all the
# entries in the corresponding list element (or column of the table)
# are missing---that is, all the entries are set to NA.
all.missing.col <- function (d)
  sapply(d,function(x) all(is.na(x)))
  
# ----------------------------------------------------------------------
# For each row of the matrix or data frame, returns true if all the
# entries in the row are provided (not missing).
none.missing.row <- function (d)
  rowSums(is.na(d)) == 0

# ----------------------------------------------------------------------
# Return a vector containing the upper triangular entries of A, or the
# off-diagonal entries of symmetric matrix A.
offdiag <- function (A)
  return(A[upper.tri(A)])

# ----------------------------------------------------------------------
# Returns the matrix product A*A'.
matrix.square <- function (A)
  return(A %*% t(A))
