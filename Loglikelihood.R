# ============================================================================ #
# Loglikelihood.R — likelihoods and radial-spline machinery
# ============================================================================ #
# The objective functions that fit the convex-mixture model, plus the radial
# spline machinery they rely on. Each ordinal score distribution is a convex
# mixture (1 - lambda) F0 + lambda Finf between the unexposed and maximal-effect
# profiles, so fitting reduces to estimating the weight lambda(.):
#   loglik          single-exposure curve — lambda_1(d1), lambda_2(d2), tau2
#                   estimated jointly; self-contained (everything is an argument).
#   loglik_joint_fixed  joint-exposure likelihood for a FIXED weight surface — scores
#                       the additive nulls so they compare to the radial fit.
#   loglik_joint    joint-exposure likelihood for the radial spline surface (H1),
#                   penalized for monotonicity and across-ray roughness.
# The radial fit works in aligned radial coordinates (rho, nu); prep_joint_meas /
# ray_weights / prep_penalty_grid / eval_deriv_grid precompute the bases and the
# derivative penalty so each optimizer step stays cheap. Formulation is in the
# paper (Section 3); the joint helpers read the global model state.
# ============================================================================ #
# loglik() (single-exposure curve) is self-contained: everything is an argument.
# The joint-exposure helpers below consume MODEL STATE from two files:
#   from Data_setup.R (fixed scaffolding):
#     dInf_BaP, d1_MAX, d2_MAX                    dose ranges and region bounds
#     base_Delta, base_Phi, deriv_Delta, deriv_Phi  radial spline bases
#     basis_BaP, R_MAX                            reference (BaP) basis / scale
#     K_phi, K_delta                              spline dimensions
#   from Run_fit.R (fitted single-exposures):
#     g, g_inv, fit_tau2                          dose alignment and its inverse
#     fit_w_BaP                                   reference (BaP) weight vector
# Source Utils.R and Data_setup.R before this file; the joint helpers additionally
# require Run_fit.R (fitted state) before they are called.
# ============================================================================ #

# Single-exposure log-likelihood (simultaneous fit) ----------------------------
#
# Negative log-likelihood for the two single-exposure weight functions, fit jointly.
# Both single-exposure data share the extremal profiles F0/Finf and the Pb->BaP dose
# alignment g(.), so BaP weights, Pb weights, and the Pb potency scalar tau2 are
# estimated in one optimization. Minimizes over `params` (via subplex).
#
# Args:
#   params  : initial values for free parameters, length (nI1-1) + (nI2-1) + 1, 
#             laid out as [ BaP spline logits | Pb spline logits | logit(tau2) ]. 
#             The final spline weight of each set is fixed to 0
#             (softmax reference), hence nI-1 free logits per chemical
#   y1, y2  : count matrices (ny x nx), columns = doses in x1/x2, rows = scores.
#   x1, x2  : unique dose values for BaP / Pb (including 0 and the max dose).
#   nI1,nI2 : number of I-spline basis functions for BaP / Pb.
#   basisI1, basisI2 : I-spline basis constructors, basisI(x) -> (length(x) x nI) matrix.
#   dInf1   : max BaP dose D1^{(max)}, used to place aligned Pb doses on the BaP axis.
#   F0, Finf: CDFs (length ny) at zero and maximal exposure; the mixture endpoints.
#   eps     : floor on probabilities before the log, for numerical stability.
#
# Returns:
#   Scalar negative log-likelihood (sum over both chemicals)
#
# Details:
#   w1/w2 are softmax-normalized spline weights; tau2 = invlogit(last param).
#   Weights h = basisI %*% w give the mixing weight lambda(d) per dose; Pb doses
#   enter via the aligned argument dInf1 * tau2 * g(d2). get_probs() turns each
#   weight into the convex mixture (1-h) F0 + h Finf as pmf increments (D0/Dinf),
#   and the log-probabilities are weighted by the observed counts.
#
loglik <- function(params,
                   y1, x1, nI1, basisI1, dInf1,
                   y2, x2, nI2, basisI2,
                   F0, Finf, eps = 1e-10) {

  # separate parameters
  params1 <- params[1:(nI1 - 1)]                  # nI1 - 1
  params2 <- params[(nI1 - 1) + 1:(nI2 - 1)]      # nI2 - 1

  # reconstruct normalized weights
  w1 <- softmax(c(params1, 0))                    # nI1
  w2 <- softmax(c(params2, 0))                    # nI2

  # reconstruct normalized shift
  tau2 <- invlogit(params[nI1 + nI2 - 1])

  # pmf increments at the extremal profiles
  D0   <- F0   - c(0, F0[-length(F0)])            # ny
  Dinf <- Finf - c(0, Finf[-length(Finf)])        # ny

  # mixture weights
  h_d1 <- as.vector(basisI1(x1) %*% w1)                    # nx1
  g_d2 <- as.vector(basisI2(x2) %*% w2)                    # nx2
  h_d2 <- as.vector(basisI1(dInf1 * tau2 * g_d2) %*% w1)   # nx2

  # probabilities
  prob1 <- get_probs(h_d1, D0, Dinf)              # ny1 x nx1
  prob2 <- get_probs(h_d2, D0, Dinf)              # ny2 x nx2

  # log probabilities
  logP1 <- log(pmax(prob1, eps))                  # ny1 x nx1
  logP2 <- log(pmax(prob2, eps))                  # ny2 x nx2

  logL <- sum(y1 * logP1) + sum(y2 * logP2)

  -logL
}

# Penalty on radial-spline derivatives -----------------------------------------
#
# Three-stage pipeline for the monotonicity penalty on lambda_12: make_grid()
# lays out the evaluation points, prep_penalty_grid() precomputes the spline
# bases there once, and eval_deriv_grid() evaluates the partial derivatives of
# lambda_12 at those points for a given parameter matrix (called every optimizer
# step). Splitting the fixed precomputation from the per-step evaluation is what
# keeps the penalized fit computationally cheaper.

## make_grid --------------------------------------------------------------------
## 
## Build the set of (d1, d2) dose points at which the derivative penalty is
## enforced. 
## 
## Args:
##   d1_MAX, d2_MAX : upper dose bounds for BaP (d1) and Pb (d2).
##   fit_LTR   : if TRUE, keep only the lower-triangular region actually
##                    probed (d2 <= d2_MAX - (d2_MAX/d1_MAX) * d1).
##   n_out          : number of grid steps per axis before filtering.
## 
## Returns:
##   Data frame with columns Var1 (d1) and Var2 (d2), the (0,0) point removed,
##   optionally restricted to the triangle.
## 
## State:
##   Reads dInf_BaP, fit_tau2 (BaP-equivalent scaling) and g_inv (Pb dose that
##   maps to a given aligned BaP dose).
##
make_grid <- function(d1_MAX, d2_MAX, fit_LTR = TRUE, n_out = 35) {

  u_grid <- seq(0, d1_MAX, length.out = n_out)
  u_grid <- (u_grid / d1_MAX)^4 * (d1_MAX)
  v_grid <- sapply(u_grid / (dInf_BaP * fit_tau2), g_inv)

  pen_grid <- expand.grid(u_grid, v_grid)
  pen_grid <- subset(pen_grid, Var2 + Var1 > 0)

  if (fit_LTR) {
    pen_grid <- subset(pen_grid, Var2 <= d2_MAX - (d2_MAX / d1_MAX) * Var1)
  }

  pen_grid
}

## prep_penalty_grid -----------------------------------------------------------
##
## Precompute everything about the penalty grid that does NOT depend on the
## parameters: aligned coordinates (rho, nu) and the four spline design matrices
## (Delta/Phi bases and their derivatives). Done once, then reused by
## eval_deriv_grid at every optimizer step. Points with non-finite nu or 
## rho <= 0 are dropped up front.
##
## Args:
##   pen_grid              : output of make_grid (columns Var1 = d1, Var2 = d2).
##   base_Delta, deriv_Delta : B-spline basis over nu and its derivative.
##   base_Phi, deriv_Phi   : I-spline basis over rho and its derivative.
##   g                     : Pb->BaP alignment function.
##   dInf_BaP, fit_tau2    : BaP-equivalent scaling for the aligned Pb coordinate.
##
## Returns:
##   List (the "pen_cache") with, over the retained points:
##     d1, d2       : original doses.
##     x, y         : aligned coordinates, x = d1, y = dInf_BaP*fit_tau2*g(d2).
##     rho, nu      : total aligned dose x+y, and mixing proportion y/rho.
##     DeltaMat, dDeltaMat  : Delta basis over nu and its derivative.
##     PhiMat, dPhiMat      : Phi basis over rho and its derivative.
##
## State:
##   Depends only on its arguments (dInf_BaP, fit_tau2, g are passed in, not 
##   read from the global model state).
##
prep_penalty_grid <- function(pen_grid, base_Delta, deriv_Delta,
                              base_Phi, deriv_Phi, g, dInf_BaP, fit_tau2) {

  x   <- pen_grid$Var1
  y   <- dInf_BaP * fit_tau2 * g(pen_grid$Var2)   # vectorized
  rho <- x + y
  nu  <- y / rho

  ok  <- is.finite(nu) & rho > 0
  x <- x[ok]; y <- y[ok]; rho <- rho[ok]; nu <- nu[ok]

  list(
    d1 = pen_grid$Var1[ok], d2 = pen_grid$Var2[ok],
    x = x, y = y, rho = rho, nu = nu,
    DeltaMat  = base_Delta(nu),
    dDeltaMat = deriv_Delta(nu),
    PhiMat    = base_Phi(rho),
    dPhiMat   = deriv_Phi(rho)
  )
}

## eval_deriv_grid -------------------------------------------------------------
##
## Evaluate the partial derivatives of lambda_12 with respect to d1 and d2 at
## every cached grid point, for a given parameter matrix. 
##
## Args:
##   par_mat   : (K_phi-1) x K_delta matrix of radial-spline parameters.
##   pen_cache : output of prep_penalty_grid.
##
## Returns:
##   n x 2 matrix; column 1 = d lambda_12 / d d1, column 2 = d lambda_12 / d d2,
##   one row per cached grid point.
##
## Details:
##   bMat = rowSoftmax0(Delta %*% par) are the per-ray mixing weights;
##   dbMat is their derivative in nu.
##
eval_deriv_grid <- function(par_mat, pen_cache) {

  logits_free <- pen_cache$DeltaMat %*% t(par_mat)

  bMat <- rowSoftmax0(logits_free)

  wd_free <- pen_cache$dDeltaMat %*% t(par_mat)
  wdMat   <- cbind(wd_free, 0)

  dphi_rho <- rowSums(pen_cache$dPhiMat * bMat)

  dbMat <- bMat * (wdMat - rowSums(bMat * wdMat))

  dphi_nu <- rowSums(pen_cache$PhiMat * dbMat)

  cbind(dphi_rho * pen_cache$rho^2 + dphi_nu * pen_cache$x,
        dphi_rho * pen_cache$rho^2 - dphi_nu * pen_cache$y)
}

# Joint-exposure log-likelihood — fixed weight function ------------------------
#
# Negative log-likelihood of the joint-exposure counts under a FIXED weight
# surface l_func (no free parameters).
#
# Args:
#   l_func  : weight function lambda_12(d1, d2) -> [0,1], evaluated per dose pair.
#   x_all   : dose combinations, nx x 2 matrix (column 1 = BaP d1, column 2 = Pb d2).
#   y_all   : count matrix ny x nx (rows = ordinal scores, columns = dose pairs).
#   D0, Dinf: pmf increments (length ny) of the extremal profiles F0 / Finf.
#   eps     : floor on probabilities before the log, for numerical stability.
#
# Returns:
#   Scalar negative log-likelihood.
#
# Details:
#   h_d1d2 is the per-dose mixing weight from l_func; get_probs turns it into the
#   convex mixture (1-h) D0 + h Dinf, and the log-probabilities are weighted by
#   the observed counts.
#
# State:
#   None directly — but l_func (the null models) reads the global model state.
#
loglik_joint_fixed <- function(l_func, x_all, y_all, D0, Dinf, eps = 1e-10) {

  # mixture weights
  h_d1d2 <- apply(x_all, 1, function(vec) l_func(vec[1], vec[2]))   # nx

  # probabilities
  probs <- get_probs(h_d1d2, D0, Dinf)            # ny x nx

  # log probabilities
  logPs <- log(pmax(probs, eps))                  # ny x nx

  logL <- sum(y_all * logPs)

  -logL
}

# Joint-exposure log-likelihood — radial-splines weight function ---------------
#
# Fit-time pipeline for the radial spline surface lambda_12: prep_joint_meas()
# caches the parameter-independent bases at the observed doses once,
# ray_weights() maps a parameter matrix to the per-dose weights, and
# loglik_joint() wraps them into the penalized objective handed to the optimizer.

## prep_joint_meas --------------------------------------------------------------
##
## Precompute everything about the observed joint-exposure doses that does NOT
## depend on the parameters: aligned coordinates (rho, nu) and the Delta/Phi
## design matrices. Built once, then reused by ray_weights at every optimizer
## step. The counts z are carried alongside so the objective has its target.
##
## Args:
##   d1d2_combos    : dose combinations, nx x 2 (column 1 = BaP, column 2 = Pb).
##   y_all_adj    : count matrix ny x nx aligned to d1d2_combos columns.
##   fit_LTR : if TRUE, keep only doses with d1 < d1_MAX and d2 < d2_MAX
##                  (the probed region), subsetting both doses and counts.
##
## Returns:
##   List ("data" / "joint_meas") with:
##     z        : count matrix ny x nx (possibly subset).
##     x, y     : aligned coordinates, x = d1, y = dInf_BaP*fit_tau2*g(d2).
##     rho, nu  : total aligned dose x+y and mixing proportion y/rho
##                (nu set to 0 where rho = 0, i.e. the origin).
##     DeltaMat : Delta basis over nu   (depends only on nu_i).
##     PhiMat   : Phi basis over rho    (depends only on rho_i).
##
## State:
##   Reads d1_MAX, d2_MAX (triangle bounds), and dInf_BaP, fit_tau2, g,
##   base_Delta, base_Phi (alignment and bases).
##
prep_joint_meas <- function(d1d2_combos, y_all_adj, fit_LTR = TRUE) {

  if (fit_LTR) {
    idx_keep  <- (d1d2_combos[, 1] < d1_MAX) & (d1d2_combos[, 2] < d2_MAX)
    d1d2_combos <- d1d2_combos[idx_keep, ]
    y_all_adj <- y_all_adj[, idx_keep]
  }

  joint_meas <- list()

  joint_meas$z <- y_all_adj                       # ny x nx  (d1d2_combos: nx x 2)

  joint_meas$x <- d1d2_combos[, 1]
  joint_meas$y <- dInf_BaP * fit_tau2 * g(d1d2_combos[, 2])

  joint_meas$rho <- joint_meas$x + joint_meas$y
  joint_meas$nu  <- joint_meas$y / joint_meas$rho
  joint_meas$nu[(!is.finite(joint_meas$nu)) | is.na(joint_meas$nu)] <- 0

  joint_meas$DeltaMat <- base_Delta(joint_meas$nu)   # depends only on nu_i
  joint_meas$PhiMat   <- base_Phi(joint_meas$rho)    # depends only on rho_i

  joint_meas
}

## ray_weights -----------------------------------------------------------------
##
## Map a radial-spline parameter matrix to the joint weight lambda_12 at each
## cached dose. Along each ray, softmax mixing weights over the Phi basis give the
## effective dose index phi(rho, nu); this is rescaled to the BaP axis and pushed
## through the reference weight function l1, so lambda_12 = l1(R_MAX * phi).
##
## Args:
##   mat_W : (K_phi-1) x K_delta parameter matrix.
##   data  : output of prep_joint_meas (uses DeltaMat and PhiMat).
##
## Returns:
##   Numeric vector of weights (length nx), one per cached dose combination.
##
## State:
##   Reads basis_BaP, fit_w_BaP (reference BaP weight function) and R_MAX.
##
ray_weights <- function(mat_W, data) {
  exp_W    <- data$DeltaMat %*% t(mat_W)
  logits_W <- rowSoftmax0(exp_W)
  phi      <- rowSums(data$PhiMat * logits_W)

  as.vector(basis_BaP(R_MAX * phi) %*% fit_w_BaP)    # nx
}

## loglik_joint -----------------------------------------------------------------
##
## Penalized negative log-likelihood minimized to fit the radial spline surface
## (H1). Combines the joint-exposure fit with two optional regularizers: a
## monotonicity penalty on negative (d1, d2) derivatives, and a roughness penalty
## on the across-ray mixing weights. This is the objective passed to subplex.
##
## Args:
##   params       : flat parameter vector, reshaped to (K_phi-1) x K_delta.
##   data         : output of prep_joint_meas (weights target z + cached bases).
##   pen_cache    : output of prep_penalty_grid (only used when gamma_pen > 0).
##   gamma_smooth : roughness penalty strength on second differences of the
##                  mixing logits across rays (0 disables it).
##   gamma_pen    : monotonicity penalty strength on negative partial derivatives
##                  of lambda_12 over the penalty grid (0 disables it).
##   D0, Dinf     : pmf increments (length ny) of F0 / Finf.
##   eps          : floor on probabilities before the log.
##
## Returns:
##   Scalar objective, -logL + pen + smooth.
##
## Details:
##   ray_weights gives the per-dose weights; get_probs forms the convex mixture
##   and the counts z weight the log-probabilities. The penalty squares the worst
##   (most negative) derivative in each direction; the smoothness term penalizes
##   curvature of the mixing weights across rays.
##
## State:
##   Reads K_phi, K_delta (parameter reshape). Downstream ray_weights and
##   eval_deriv_grid read further model state (see their headers).
##
loglik_joint <- function(params, data, pen_cache,
                         gamma_smooth, gamma_pen,
                         D0, Dinf, eps = 1e-10) {

  mat_W <- matrix(params, nrow = K_phi - 1, ncol = K_delta)

  # mixture weights
  h_d1d2 <- ray_weights(mat_W, data)              # nx

  # probabilities
  probs <- get_probs(h_d1d2, D0, Dinf)            # ny x nx

  # log probabilities
  logPs <- log(pmax(probs, eps))                  # ny x nx

  logL <- sum(data$z * logPs)

  # penalty on negative partial derivatives
  pen <- 0
  if (gamma_pen > 0) {
    g_dxy <- eval_deriv_grid(mat_W, pen_cache)    # grid evaluation
    pen <- gamma_pen * (max(0, -min(g_dxy[, 1]))^2 + max(0, -min(g_dxy[, 2]))^2)
  }

  # smoothness across rays
  smooth <- 0
  if (gamma_smooth > 0) {
    smooth <- gamma_smooth *
      sum(apply(mat_W[, -1, drop = FALSE], 1, function(b) sum(diff(b, differences = 2)^2)))
  }

  -logL + pen + smooth
}
