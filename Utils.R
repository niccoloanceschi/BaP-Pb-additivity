# ==============================================================================
# Utils.R — package imports and generic helper functions
# ==============================================================================

# Imports ----------------------------------------------------------------------
library(splines2)   # iSpline, bSpline
library(subplex)    # subplex() optimizer
library(ggplot2)    # plotting
library(tidyr)      # reshaping dataframes
library(tictoc)     # runtimes
library(patchwork)  # composite plots

# Link functions ---------------------------------------------------------------
invlogit <- function(z) 1 / (1 + exp(-z))
logit    <- function(p) log(p / (1 - p))  

# Softmax helpers --------------------------------------------------------------
softmax <- function(x) {
  x  <- as.numeric(x)
  x  <- x - max(x)          # numerical stability
  ex <- exp(x)
  ex / sum(ex)
}

# Row-wise softmax of [A_free | 0]: append a zero reference column, then
# normalize each row. Returns an n x (ncol(A_free) + 1) matrix.
rowSoftmax0 <- function(A_free) {
  A  <- cbind(A_free, 0)               # n x K
  mx <- apply(A, 1, max)               # length n
  A  <- A - mx
  E  <- exp(A)
  E / rowSums(E)
}

# Mixture / inversion helpers --------------------------------------------------

# Convex-mixture probabilities: for each weight w, (1 - w) * D0 + w * Dinf.
# Returns an ny x length(weights) matrix.
get_probs <- function(weights, D0, Dinf) {
  D0 %x% t(rep(1, length(weights))) + (Dinf - D0) %x% t(weights)   # ny x nx
}

# Invert a monotone I-spline weight function f(x) = basisI(x) %*% w on [xMin, xMax].
invert_f_uniroot <- function(lambda, basisI, w, xMin, xMax, tol = 1e-10) {

  my_f <- function(x, basisI, w) as.vector(basisI(x) %*% w)

  fL <- my_f(xMin, basisI, w)
  fU <- my_f(xMax, basisI, w)

  # clamp if outside range
  if (lambda == fL) return(xMin)
  if (lambda < fL)  return(NA)
  if (lambda > fU)  return(NA)

  uniroot(function(x) my_f(x, basisI, w) - lambda,
          lower = xMin, upper = xMax, tol = tol)$root
}

# Clamp values in interval
.clamp <- function(x, a, b) pmin(b, pmax(a, x))
