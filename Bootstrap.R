# ============================================================================ #
# Bootstrap_test.R — parametric bootstrap for the additivity test (REFERENCE)
# ============================================================================ #
# This script documents the parametric-bootstrap calibration of the likelihood-
# ratio test of additivity (H0: additive null vs H1: radial fit). 
# The bootstrap requires refitting the radial model on each of many replicates, 
# which is computationally intensive. The analysis in the manuscript drew 10,000 
# replicates per null model and alternative hypothesis, and was run in parallel
# across multiple nodes on a cluster (later aggregated).
#
# This script presents that procedure
# as a single linear pass for readability; it is provided for transparency and is
# not sourced by the analysis pipeline.
# The analysis is mainly divided in three stages:
#   1. simulate  — per (null model, seed) job, bootstrap replicates and refit
#                  the radial model, saving per-replicate log-likelihoods.
#   2. aggregate — collect the per-job outputs into per-hypothesis tables.
#   3. summarize — form the LR test statistics, p-values, and power.
# ============================================================================ #
# Paths are repo-relative via source_dir Requires the same model state as the main
# pipeline (Data_setup.R, Null_models.R, Loglikelihood.R, Model_fit.R).
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# 0. SETUP — baseline quantities and the observed test statistic -------------
# ---------------------------------------------------------------------------- #
# Prepares everything the bootstrap needs from the already-fitted model state:
# the extremal-profile PMFs, the probed triangular region and its joint-only
# core, the reference CDFs under each hypothesis (the data-generating profiles
# for the replicates), and the optimizer/penalty settings. Also computes the
# observed LR statistic T_obs from the actual fit, against which the bootstrap
# null distributions are later compared.

## hyper-parameters 
which_HP <- "H0_bliss" # must be one of "H1", "H0_bliss", "H0_mean", and "H0_hand"
mc_seed  <- 123
nRep     <- 10000
ga_pen0  <- -2

if (!which_HP %in% c("H1", "H0_bliss", "H0_mean", "H0_hand"))
  stop("which_HP must be one of 'H1','H0_bliss','H0_mean','H0_hand'")

set.seed(mc_seed)

## baseline quantities -------------------------------------------------------- #
# pmf increments at the extremal profiles
D0   <- F0   - c(0, F0[-length(F0)])       # ny
Dinf <- Finf - c(0, Finf[-length(Finf)])   # ny

## triangular (probed) region and its joint-only core
idx_tri   <- (d1d2_meas[, 1] < d1_MAX) & (d1d2_meas[, 2] < d2_MAX)
x_all_tri <- d1d2_meas[idx_tri, ]
y_all_tri <- y_all_adj[, idx_tri]
F_all_tri <- F_all_adj[, idx_tri]

idx_core   <- apply(x_all_tri, 1, min) > 0
x_core_tri <- x_all_tri[idx_core, ]
y_core_tri <- y_all_tri[, idx_core]
F_core_tri <- F_all_tri[, idx_core]

data_core <- prep_joint_meas(x_core_tri, y_core_tri)

## reference CDFs under each hypothesis (data-generating profile for the bootstrap)
h_d1d2_mean  <- apply(x_core_tri, 1, function(v) l_mean(v[1], v[2]))
h_d1d2_bliss <- apply(x_core_tri, 1, function(v) l_bliss(v[1], v[2]))
h_d1d2_hand  <- apply(x_core_tri, 1, function(v) l_hand(v[1], v[2]))
h_d1d2_fit   <- ray_weights(par_fitJ, data_core)

F_mean  <- get_probs(h_d1d2_mean,  F0, Finf)
F_bliss <- get_probs(h_d1d2_bliss, F0, Finf)
F_hand  <- get_probs(h_d1d2_hand,  F0, Finf)
F_fit   <- get_probs(h_d1d2_fit,   F0, Finf)

## optimizer / penalty setup
nCat     <- nrow(y_core_tri)
ga_pen   <- sum(y_all_adj) * (10^ga_pen0)
p        <- (K_phi - 1) * K_delta
par0     <- rep(0, p)
parscale <- rep(0.5, p)

logL_all <- logL_core <- matrix(NA, ncol = 5, nrow = nRep)
colnames(logL_all) <- colnames(logL_core) <- c("hand", "mean", "bliss", "fit", "fit_pen")
conv_check <- rep(NA, nRep)

## observed test statistic from the actual fit (already in the model state)
## `logL_vals` is produced when the log-likelihood comparison is computed
## (see plot_loglik / Model_fit.R); it holds logL per (data subset, method).
lv <- logL_vals
idx_fit   <- lv$data == "core" & lv$meth == "fit"
idx_bliss <- lv$data == "core" & lv$meth == "bliss"
idx_mean  <- lv$data == "core" & lv$meth == "mean"
idx_hand  <- lv$data == "core" & lv$meth == "hand"

T_Fit_bliss <- -2 * (lv$logL[idx_fit] - lv$logL[idx_bliss])
T_Fit_mean  <- -2 * (lv$logL[idx_fit] - lv$logL[idx_mean])
T_Fit_hand  <- -2 * (lv$logL[idx_fit] - lv$logL[idx_hand])

# ---------------------------------------------------------------------------- #
# 1. SIMULATE — bootstrap replicates and refits ------------------------------
# ---------------------------------------------------------------------------- #
# Looping over the four hypotheses, storing the core-dose
# log-likelihoods per hypothesis in `ris_core` (no files written).

ris_core <- list()                              # one logL_core matrix per hypothesis

for (which_HP in c("H1", "H0_bliss", "H0_mean", "H0_hand")) {
  
  F_ref <- switch(which_HP,
                  "H0_mean"  = F_mean,
                  "H0_bliss" = F_bliss,
                  "H0_hand"  = F_hand,
                  "H1"       = F_fit)
  
  logL_all <- logL_core <- matrix(NA, ncol = 5, nrow = nRep,
                                  dimnames = list(NULL, c("hand","mean","bliss","fit","fit_pen")))
  conv_check <- rep(NA, nRep)
  
  ## bootstrap loop ---------------------------------------------------------- #
  for (rr in seq_len(nRep)) {
    
    ## draw a bootstrap replicate from the reference CDFs on the core doses
    y_core_sim <- sapply(seq_len(ncol(y_core_tri)), function(ii) {
      raw <- sample.int(nCat, replace = TRUE, size = sum(y_core_tri[, ii]),
                        prob = c(F_core_tri[1, ii], diff(F_ref[, ii])))
      table(factor(raw - 1, levels = 0:(nCat - 1)))
    })
    colnames(y_core_sim) <- colnames(y_core_tri)
    
    ## embed the simulated core back into the full triangular counts
    y_all_sim <- y_all_tri
    y_all_sim[, idx_core] <- y_core_sim
    
    data_sim  <- prep_joint_meas(x_all_tri,  y_all_sim)
    data_core <- prep_joint_meas(x_core_tri, y_core_sim)
    pen_cache <- prep_penalty_grid(pen_grid_v0, base_Delta, deriv_Delta,
                                   base_Phi, deriv_Phi, g, dInf_BaP, fit_tau2)
    
    ## refit the radial model on the replicate
    fitJ_sim <- subplex(par_fitJ, loglik_joint,
                        data = data_sim, pen_cache = pen_cache,
                        D0 = D0, Dinf = Dinf, gamma_smooth = 0, gamma_pen = ga_pen,
                        hessian = FALSE, control = list(parscale = parscale, maxit = 100000, reltol = 1e-2))
    par_simJ <- matrix(fitJ_sim$par, nrow = K_phi - 1, ncol = K_delta)
    
    ## log-likelihoods of each model on the replicate (full triangle and core)
    logL_all[rr, ] <- c(
      loglik_joint_fixed(l_hand,  x_all_tri, y_all_sim, D0, Dinf),
      loglik_joint_fixed(l_mean,  x_all_tri, y_all_sim, D0, Dinf),
      loglik_joint_fixed(l_bliss, x_all_tri, y_all_sim, D0, Dinf),
      loglik_joint(par_simJ, data_sim, pen_cache, 0, ga_pen, D0, Dinf),
      loglik_joint(par_simJ, data_sim, pen_cache, 0, 0,      D0, Dinf))
    logL_core[rr, ] <- c(
      loglik_joint_fixed(l_hand,  x_core_tri, y_core_sim, D0, Dinf),
      loglik_joint_fixed(l_mean,  x_core_tri, y_core_sim, D0, Dinf),
      loglik_joint_fixed(l_bliss, x_core_tri, y_core_sim, D0, Dinf),
      loglik_joint(par_simJ, data_core, pen_cache, 0, ga_pen, D0, Dinf),
      loglik_joint(par_simJ, data_core, pen_cache, 0, 0,      D0, Dinf))
    
    conv_check[rr] <- fitJ_sim$convergence == 0
  }
  
  ris_core[[which_HP]] <- logL_core
}

# ---------------------------------------------------------------------------- #
# 2. AGGREGATE — collect the per-hypothesis replicate log-likelihoods --------
# ---------------------------------------------------------------------------- #
# In the reported analysis each hypothesis/seed was a separate cluster job whose
# output was saved and later row-bound. Presented here as an in-memory list, one
# entry per hypothesis, holding the core-dose log-likelihoods used below.
# (Stage 1 above computes `logL_core` for a single `which_HP`; here we assume it
#  has been evaluated for each, e.g. by looping stage 1 over the four values.)

ris_H0_bliss <- as.data.frame(ris_core[["H0_bliss"]])
ris_H0_mean  <- as.data.frame(ris_core[["H0_mean"]])
ris_H0_hand  <- as.data.frame(ris_core[["H0_hand"]])
ris_H1       <- as.data.frame(ris_core[["H1"]])

# ---------------------------------------------------------------------------- #
# 3. SUMMARIZE — LR test statistics, p-values, and power ---------------------
# ---------------------------------------------------------------------------- #
# T = -2 (logL_fit - logL_null) on the core doses. Null distribution from the H0
# replicates; p-value = P(T_H0 >= T_obs); power = P(T_H1 >= critical value).
# NB: in ris_core the column `fit_pen` is the UNPENALIZED refit (gamma_pen = 0),
# i.e. the fitted statistic entering the LR test.

## null distributions of the test statistic
T_H0_bliss <- -2 * (ris_H0_bliss$fit_pen - ris_H0_bliss$bliss)
T_H0_mean  <- -2 * (ris_H0_mean$fit_pen  - ris_H0_mean$mean)
T_H0_hand  <- -2 * (ris_H0_hand$fit_pen  - ris_H0_hand$hand)

## critical values at level 0.05
cA_H0_bliss <- quantile(T_H0_bliss, 1 - 0.05, na.rm = TRUE)
cA_H0_mean  <- quantile(T_H0_mean,  1 - 0.05, na.rm = TRUE)
cA_H0_hand  <- quantile(T_H0_hand,  1 - 0.05, na.rm = TRUE)

## p-values (observed statistic vs its null distribution)
pV_H1_bliss <- mean(T_H0_bliss >= T_Fit_bliss, na.rm = TRUE)
pV_H1_mean  <- mean(T_H0_mean  >= T_Fit_mean,  na.rm = TRUE)
pV_H1_hand  <- mean(T_H0_hand  >= T_Fit_hand,  na.rm = TRUE)

## power (alternative statistics vs the null critical value)
T_H1_bliss <- -2 * (ris_H1$fit_pen - ris_H1$bliss)
T_H1_mean  <- -2 * (ris_H1$fit_pen - ris_H1$mean)
T_H1_hand  <- -2 * (ris_H1$fit_pen - ris_H1$hand)

power_H1_bliss <- mean(T_H1_bliss >= cA_H0_bliss, na.rm = TRUE)
power_H1_mean  <- mean(T_H1_mean  >= cA_H0_mean,  na.rm = TRUE)
power_H1_hand  <- mean(T_H1_hand  >= cA_H0_hand,  na.rm = TRUE)

## summary table (p-value and power per null model)
ris_summary <- data.frame(
  row.names = c("H0_bliss", "H0_mean", "H0_hand"),
  p_value   = c(pV_H1_bliss, pV_H1_mean, pV_H1_hand),
  power     = c(power_H1_bliss, power_H1_mean, power_H1_hand)
)
ris_summary