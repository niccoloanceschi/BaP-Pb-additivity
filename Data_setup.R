# ============================================================================ #
# Data_setup.R — data import, spline bases, and data subsets
# ============================================================================ #
# Loads the preprocessed assay data and builds the fixed scaffolding the fit
# runs on, none of it parameter-dependent:
#   - single-exposure I-spline bases (basis_BaP, basis_Pb, deriv_Pb) and their sizes
#   - radial bases for the joint surface: I-spline in rho (base_Phi, deriv_Phi)
#     and B-spline in the mixing proportion nu (base_Delta, deriv_Delta)
#   - dose constants (dInf_BaP, dInf_Pb, d1_MAX, d2_MAX, R_MAX) and spline
#     dimensions (K_phi, K_delta)
#   - extremal-profile PMFs (D0, Dinf)
#   - nested joint-exposure data subsets (tri / no0 / core)
# Built once and reused unchanged; the fitted profiles and the basis caches that
# depend on them are created later in Run_fit.R.
# ============================================================================ #
# Requires: source_dir, data_file (set in Analysis.Rmd); Utils.R sourced.
# Creates the fixed MODEL STATE (bases, constants, PMFs, data subsets) that lives
# in the global environment and is consumed by Run_fit.R, Null_models.R, the
# joint helpers in Loglikelihood.R, and Plots.R.
# ============================================================================ #

# Import data ------------------------------------------------------------------

# Load the preprocessing output and expose its elements as top-level objects.
# Expects the list:
#   yBaP_input, yPb_input  single-exposure count matrices (ny x nx) for BaP / Pb
#   xBaP, xPb              their unique dose values (incl. 0 and max)
#   F0, Finf               extremal CDFs (ny): unexposed and maximal-effect profiles
#   F_all_adj, y_all_adj   batch-adjusted joint-exposure CDFs / counts (ny x nx)
#   d1d2_meas              measured joint dose pairs (nx x 2: BaP, Pb)
ris_fit <- readRDS(file.path(source_dir, data_file))
list2env(ris_fit, envir = .GlobalEnv)
rm(ris_fit)

# PMF at the two extremal profiles: p_l = F_l - F_{l-1} (mass per score)
D0   <- F0   - c(0, F0[-length(F0)])       # ny
Dinf <- Finf - c(0, Finf[-length(Finf)])   # ny

# Spline bases — single exposure -----------------------------------------------

# Max doses D1^max / D2^max, rounded up so the spline domain covers all data.
dInf_BaP <- ceiling(max(xBaP))
dInf_Pb  <- ceiling(max(xPb))

# Monotone I-spline bases for the single-exposure weights lambda_1, lambda_2.
# degree 3, interior knots hand-placed in each dose range, intercept = FALSE so
# the basis (and thus the weight) starts at 0; softmax weights then give a
# monotone non-decreasing curve on [0, 1]. deriv_Pb is the same basis' first
# derivative g'(.), needed by the Hand null.
basis_BaP <- function(x) iSpline(x, degree = 3, knots = c(1, 2),  Boundary.knots = c(0, dInf_BaP), intercept = FALSE)
basis_Pb  <- function(x) iSpline(x, degree = 3, knots = c(5, 10), Boundary.knots = c(0, dInf_Pb),  intercept = FALSE)
deriv_Pb  <- function(x) iSpline(x, degree = 3, knots = c(5, 10), Boundary.knots = c(0, dInf_Pb),  intercept = FALSE, derivs = 1)

# number of basis functions per chemical (free spline weights = nI - 1 after softmax)
nI_BaP <- ncol(basis_BaP(0))
nI_Pb  <- ncol(basis_Pb(0))

# Spline bases — radial (joint exposure) ---------------------------------------

d1_MAX <- 3.5        # BaP dose bound for the fitted surface
d2_MAX <- 120        # Pb dose bound
R_MAX  <- dInf_BaP   # max aligned total dose rho (reference = BaP scale)

# Radial component: I-spline in the total aligned dose rho (monotone along a ray),
# with deriv_Phi its derivative for the monotonicity penalty.
base_Phi  <- function(x) iSpline(x, degree = 3, knots = c(1, 2), Boundary.knots = c(0, dInf_BaP), intercept = FALSE)
deriv_Phi <- function(x) iSpline(x, degree = 3, knots = c(1, 2), Boundary.knots = c(0, dInf_BaP), intercept = FALSE, derivs = 1)
K_phi <- ncol(base_Phi(0))

# Angular component: B-spline in the mixing proportion nu on [0, 1], with a
# leading constant column (cbind 1) for the intercept; its derivative starts with
# 0 accordingly. Governs how the ray weights vary across rays.
base_Delta  <- function(x) cbind(1, bSpline(x, degree = 3, Boundary.knots = c(0, 1), intercept = TRUE))
deriv_Delta <- function(x) cbind(0, bSpline(x, degree = 3, Boundary.knots = c(0, 1), intercept = TRUE, derivs = 1))
K_delta <- ncol(base_Delta(0))

# Joint-exposure data subsets --------------------------------------------------
# Subset the measured joint doses into nested views (tri / no0 / core), each as a
# matched triple x_* (doses) / y_* (counts) / F_* (empirical CDFs):
#   `tri`  : probed lower-rectangular region (d1 < d1_MAX, d2 < d2_MAX)
#   `no0`  : tri with the (0,0) control dropped   (row max > 0)
#   `core` : tri joint exposures only, both doses > 0   (row min > 0)

idx_tri   <- (d1d2_meas[, 1] < d1_MAX) & (d1d2_meas[, 2] < d2_MAX)
x_all_tri <- d1d2_meas[idx_tri, ]
y_all_tri <- y_all_adj[, idx_tri]
F_all_tri <- F_all_adj[, idx_tri]

idx_no0   <- apply(x_all_tri, 1, max) > 0
x_tri_no0 <- x_all_tri[idx_no0, ]
y_tri_no0 <- y_all_tri[, idx_no0]
F_tri_no0 <- F_all_tri[, idx_no0]

idx_core   <- apply(x_all_tri, 1, min) > 0
x_core_tri <- x_all_tri[idx_core, ]
y_core_tri <- y_all_tri[, idx_core]
F_core_tri <- F_all_tri[, idx_core]









