# This file defines a number of functions that do not fit under any
# one category.
#
# SYMBOL DEFINITIONS
# ----------------------------------------------------------------------
# Shorthand for machine epsilon.
eps <- .Machine$double.eps

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# The logit function.
logit10 <- function (x)
  log10((x + eps)/(1 - x + eps))

# ----------------------------------------------------------------------
# Projects points to the interval [a,b].
project.onto.interval <- function (x, a, b)
  pmin(b,pmax(a,x))

