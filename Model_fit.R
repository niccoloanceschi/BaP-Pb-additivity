# ============================================================================ #
# Model_fit.R — margins fit, weight functions, and joint radial-spline fit
# ============================================================================ #
# Fits the two stages of the model on the scaffolding built in Data_setup.R,
# and leaves behind the FITTED model state. Run top to bottom, it proceeds as:
#   1. fit the single-drug margins simultaneously (loglik) — recovering the BaP
#      and Pb spline weights and the potency scalar tau2.
#   2. define the fitted weight functions l1, g, l2 and their inverses.
#   3. build the parameter-independent caches (measurement + penalty grid) that
#      depend on the fitted alignment.
#   4. fit the radial spline surface (loglik_joint) under the monotonicity
#      penalty, giving par_fitJ.
# ============================================================================ #
# Requires: Utils.R, Loglikelihood.R, and Data_setup.R sourced (the latter
# provides the bases, constants, PMFs, and data subsets read below).
# Creates the FITTED model state — fit_w_BaP, fit_w_Pb, fit_tau2, l1/g/l2 and
# inverses, data_list/_no0/_core, pen_cache, par_fitJ — in the global environment,
# consumed by Null_models.R, the joint helpers in Loglikelihood.R, and Plots.R.
# ============================================================================ #

# Fit single exposures data (simultaneously for BaP and Pb) --------------------

## optimizer setup
init     <- c(rep(0, nI_BaP + nI_Pb - 2), logit(0.5)) # initial free parameters: spline weights + logit(tau2)
parscale <- rep(0.5, nI_BaP + nI_Pb - 1) # scales of the free parameters in the optimizer

cat('Fit Single-Exposure Data: ')
tic()
fit_margins <- subplex(init, loglik,
  y1 = yBaP_input, x1 = xBaP, nI1 = nI_BaP, basisI1 = basis_BaP, dInf1 = dInf_BaP,
  y2 = yPb_input,  x2 = xPb,  nI2 = nI_Pb,  basisI2 = basis_Pb,
  F0 = F0, Finf = Finf,
  hessian = FALSE, control = list(parscale = parscale, maxit = 100000, reltol = 1e-4))
toc()

## extract fitted parameters (same layout as init)
fit_par_BaP <- fit_margins$par[1:(nI_BaP - 1)]
fit_par_Pb  <- fit_margins$par[(nI_BaP - 1) + 1:(nI_Pb - 1)]

## reconstruct normalized weights and shift from the free parameters
fit_w_BaP <- softmax(c(fit_par_BaP, 0))
fit_w_Pb  <- softmax(c(fit_par_Pb,  0))
fit_tau2  <- invlogit(fit_margins$par[nI_BaP + nI_Pb - 1])

# Weight functions and inverses ------------------------------------------------

# The fitted margins as closures over the spline weights. l1 is the reference
# (BaP) weight; g maps a Pb dose to a [0,1] fraction of its aligned range; l2
# reuses l1 on the aligned Pb dose, so both drugs share one dose-response shape.
l1 <- function(d1) as.vector(basis_BaP(d1) %*% fit_w_BaP)   # lambda_1(d1) in [0,1]
g  <- function(d2) as.vector(basis_Pb(d2)  %*% fit_w_Pb)    # g(d2) in [0,1]
l2 <- function(d2) l1(dInf_BaP * fit_tau2 * g(d2))          # lambda_2(d2) = lambda_1(D1^max * tau2 * g(d2))

# Inverses by root-finding on the monotone bases (scalar in, scalar out; NA if
# the target is outside the basis range).
l1_inv <- function(y) invert_f_uniroot(y, basisI = basis_BaP, w = fit_w_BaP, xMin = 0, xMax = dInf_BaP)
g_inv  <- function(u) invert_f_uniroot(u, basisI = basis_Pb,  w = fit_w_Pb,  xMin = 0, xMax = dInf_Pb)
l2_inv <- function(y) {
  u <- l1_inv(y) / (dInf_BaP * fit_tau2)
  if (!is.finite(u) || u < 0 || u > 1) return(NA_real_)
  g_inv(u)
}

# Fit joint exposures ----------------------------------------------------------

## optimizer setup
p        <- (K_phi - 1) * K_delta # number of free parameters
par0     <- rep(0, p)             # initial free parameters: (K_phi-1) x K_delta matrix
parscale <- rep(0.5, p)           # scales of the free parameters in the optimizer
fit_LTR  <- TRUE                  # restrict the joint fit to the probed lower-triangular region

## penalty parameters setup
##   NB:  the penalty parameters are scaled by the total count so it stays 
##        commensurate with the log-likelihood term.
ga_sm  <- 0 # across-ray smoothness (no penalty used here)
ga_pen <- sum(y_all_adj) * 10^(-2) # monotonicity 

## penalty grid cache
##  make_grid lays out the points where the monotonicity penalty is enforced;
##  prep_penalty_grid precomputes their bases + derivatives once (reused every
##  optimizer step). Independent of the measured doses above.
pen_grid_v0 <- make_grid(d1_MAX, d2_MAX, fit_LTR = fit_LTR)
pen_cache   <- prep_penalty_grid(pen_grid_v0, base_Delta, deriv_Delta,
                                 base_Phi, deriv_Phi, g, dInf_BaP, fit_tau2)

## measurement caches (basis precomputed on the observed dose grid)
##  prep_joint_meas caches the parameter-independent bases per dose set, so the
##  objective can be re-evaluated cheaply. data_list is the full fit target; the
##  _no0 / _core caches feed the log-likelihood comparison across data subsets.
data_list <- prep_joint_meas(d1d2_meas, y_all_adj, fit_LTR = fit_LTR)
data_no0  <- prep_joint_meas(x_tri_no0,  y_tri_no0)
data_core <- prep_joint_meas(x_core_tri, y_core_tri)

cat('Fit Joint-Exposure Data: ')
tic()
fitJ <- subplex(par0, loglik_joint,
  data = data_list, pen_cache = pen_cache,
  D0 = D0, Dinf = Dinf,
  gamma_smooth = ga_sm, gamma_pen = ga_pen,
  hessian = FALSE, control = list(parscale = parscale, maxit = 100000, reltol = 1e-2))
toc()

## reshape to the (K_phi-1) x K_delta weight matrix
par_fitJ       <- matrix(fitJ$par, nrow = K_phi - 1, ncol = K_delta)

## Check the most negative partial derivative over the penalty grid
min_deriv_fitJ <- min(eval_deriv_grid(par_fitJ, pen_cache))
