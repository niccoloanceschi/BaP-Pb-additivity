
# ============================================================================ #
# Additive null models — lambda_12^(o)(d1, d2)
# ============================================================================ #
# The three baselines for drug additive activion, against which synergy and 
# antagonism are tested. Each combines # the single-drug weights l1, l2 under 
# a different notion of additivity:
#   Bliss  — independent action (survival factorizes).
#   Loewe+ — dose equivalence via the Explicit Mean Equation (EME).
#   Hand   — additive instantaneous effect rates along each ray.
# Derivations are in the paper (Hand: Appendix C); the code here just evaluates
# them given the fitted margins. All read the global model state.
# ============================================================================ #
# Consumes MODEL STATE from two files:
#   from Data_setup.R (fixed scaffolding):
#     dInf_BaP, dInf_Pb       maximum doses (D1^max, D2^max)
#     deriv_Pb                Pb I-spline derivative basis (Hand only)
#   from Run_fit.R (fitted margins):
#     l1, l2, g, g_inv        single-drug weight functions and inverses
#     fit_tau2                Pb potency-compression scalar
#     fit_w_Pb                Pb spline weights (Hand only)
# Source Utils.R and Data_setup.R before this file; the models additionally
# require Run_fit.R (fitted state) before they are called.
# ============================================================================ #

#  Bliss -----------------------------------------------------------------------
#
# Bliss independence: viewing l_s as a failure probability, the joint survival
# factorizes, 1 - l12 = (1 - l1)(1 - l2), i.e. l12 = l1 + l2 - l1*l2.
# Vectorized: d1, d2 may be equal-length vectors (l1/l2 act elementwise).
#
# Args:
#   d1, d2 : BaP and Pb doses (scalar or equal-length vectors).
# Returns:
#   Weight lambda_12 in [0,1], same length as the inputs.
# State:
#   Reads l1, l2.
#
l_bliss <- function(d1, d2){
  l1(d1) + l2(d2) - l1(d1) * l2(d2)
}

# Loewe+ — Explicit Mean Equation (EME) ----------------------------------------
#
# EME null (Lederer et al. 2018):
#   l_mean(d1,d2) = 1/2 [ l1{ d1 + l1^{-1}(l2(d2)) } + l2{ d2 + l2^{-1}(l1(d1)) } ]
# Each term augments one drug by the equal-effect dose of the other, using the
# dose alignment to evaluate the inverses in closed form:
#   l1^{-1}(l2(d2)) = dInf_BaP * tau2 * g(d2)
#   l2^{-1}(l1(d1)) = g^{-1}( d1 / (dInf_BaP * tau2) )
# Augmented doses are clamped to each drug's range. Returns NA when d1 exceeds
# the aligned BaP maximum (dInf_BaP*tau2), where the Pb-equivalent of l1(d1) is
# undefined.
# SCALAR ONLY: the range guards use scalar if(); call one (d1, d2) pair at a time
# (e.g. via apply over rows), unlike l_bliss / l_hand.
#
# Args:
#   d1, d2 : single BaP and Pb dose.
# Returns:
#   Weight lambda_12 in [0,1], or NA_real_ if out of the defined domain.
# State:
#   Reads dInf_BaP, dInf_Pb, fit_tau2, g, g_inv, l1, l2, and .clamp.
#
l_mean <- function(d1, d2) {
  # term 1: l1(d1 + l1^{-1}(l2(d2))) = l1(d1 + D1^max * tau2 * g(d2))
  t1_arg <- .clamp(d1 + dInf_BaP * fit_tau2 * g(d2), 0, dInf_BaP)
  term1  <- l1(t1_arg)

  # term 2: l2(d2 + l2^{-1}(l1(d1))) = l2(d2 + g^{-1}( d1 / (D1^max * tau2) ))
  u_needed <- d1 / (dInf_BaP * fit_tau2)
  if (!is.finite(u_needed) || u_needed < 0 || u_needed > 1) return(NA_real_)
  shift_d2 <- g_inv(u_needed)
  if (is.na(shift_d2)) return(NA_real_)
  t2_arg <- .clamp(d2 + shift_d2, 0, dInf_Pb)
  term2  <- l2(t2_arg)

  0.5 * (term1 + term2)
}

# Hand model -------------------------------------------------------------------
#
# The Hand null adds instantaneous effect rates along each fixed-ratio ray,
# working in the UNALIGNED space (nu = d1/(d1+d2), rho = d1+d2). hand_tables()
# precomputes, per ray, the mapping between effect level eta and total dose rho;
# l_hand() then reads lambda_12 off that mapping by interpolation.

## hand_tables ------------------------------------------------------------------
##
## For a set of ray proportions nu_u, tabulate rho = psi^{-1}(eta) via the
## simplified Hand integral of Appendix C:
##   psi^{-1}(l2(s)) = \int_0^s  A g'(sigma) / ( nu + (1-nu) A g'(sigma) ) d sigma,
## with A = dInf_BaP * tau2. The integrand is evaluated on a grid of the Pb-
## equivalent coordinate s and integrated by the trapezoidal rule; the effect
## grid eta = l2(s) is shared across rays (0 -> eta2_max).
##
## Args:
##   nu_u    : unique ray proportions nu = d1/(d1+d2) (length K).
##   n_grid  : number of integration nodes over s in [eta_thr, dInf_Pb].
##   eta_thr : small lower bound on s, avoiding the degenerate s = 0 endpoint.
## Returns:
##   List with
##     ETA : effect grid (length n_grid), shared across rays.
##     RHO : n_grid x K matrix of cumulative total dose; RHO[, k] is increasing in
##           ETA, and RHO[n_grid, k] = rho_max(nu_k), the largest reachable dose
##           on ray k (beyond which Hand is undefined, eta > eta2_max).
## State:
##   Reads dInf_Pb, dInf_BaP, fit_tau2, deriv_Pb, fit_w_Pb, l2.
## 
hand_tables <- function(nu_u, n_grid = 200L, eta_thr = 1e-10) {
  s   <- seq(eta_thr, dInf_Pb, length.out = n_grid)
  G2  <- dInf_BaP * fit_tau2 * as.vector(deriv_Pb(s) %*% fit_w_Pb)     # A * g'(s) >= 0
  eta <- l2(s)                                               # shared; 0 -> eta2_max
  Fm  <- outer(G2, nu_u, function(g, v) g / (v + (1 - v) * g))         # n_grid x K
  dR  <- 0.5 * diff(s) * (Fm[-1, , drop = FALSE] + Fm[-n_grid, , drop = FALSE])
  rho <- rbind(0, apply(dR, 2, cumsum))                               # n_grid x K
  list(ETA = eta, RHO = rho)                                 # RHO[n_grid, ] = rho_max(nu)
}

## l_hand ----------------------------------------------------------------------
##
## Evaluate the Hand null lambda_12 at dose pairs (d1, d2). Pure-drug rays reduce
## to the marginals l1 / l2; mixed rays are handled by building the eta<->rho
## tables once per unique nu (via hand_tables) and linearly interpolating eta at
## the queried total dose rho = d1 + d2. Points beyond rho_max(nu) return NA
## (Hand undefined there — NOT clamped to 1).
## Vectorized: d1, d2 may be equal-length vectors.
##
## Args:
##   d1, d2  : BaP and Pb doses (scalar or equal-length vectors).
##   n_grid  : integration resolution, forwarded to hand_tables.
##   eta_thr : ray-purity tolerance (nu within eta_thr of 0 or 1 is treated as a
##             marginal) and integration lower bound.
## Returns:
##   Weight lambda_12 in [0,1] (NA beyond the Hand domain), same length as inputs.
##
## State:
##   Reads l1, l2, dInf_BaP, dInf_Pb; calls hand_tables.
##
l_hand <- function(d1, d2, n_grid = 200L, eta_thr = 1e-10) {

  rho <- d1 + d2
  nu  <- ifelse(rho > 0, d1 / rho, 1)
  out <- numeric(length(rho))

  isB <- nu >= 1 - eta_thr                    # pure BaP -> marginal (defined up to 1)
  isP <- nu <= eta_thr                        # pure Pb  -> marginal (up to eta2_max)

  if (any(isB)) out[isB] <- l1(pmin(d1[isB], dInf_BaP))
  if (any(isP)) out[isP] <- l2(pmin(d2[isP], dInf_Pb))

  mix <- which(!isB & !isP)
  if (length(mix)) {
    nu_u <- unique(nu[mix])
    col  <- match(nu[mix], nu_u)
    tb   <- hand_tables(nu_u, n_grid, eta_thr)
    ETA  <- tb$ETA; RHO <- tb$RHO
    L    <- nrow(RHO); m <- length(mix)
    rr   <- rho[mix]

    sub   <- RHO[, col, drop = FALSE]              # L x m
    j     <- colSums(sub <= rep(rr, each = L))     # 0..L
    undef <- j >= L                                # rho > rho_max(nu): Hand undefined
    j0    <- pmin(pmax(j, 1L), L - 1L)
    ii <- seq_len(m)
    rl <- sub[cbind(j0, ii)]
    rh <- sub[cbind(j0 + 1L, ii)]
    el <- ETA[j0]
    eh <- ETA[j0 + 1L]
    t  <- ifelse(rh > rl, (rr - rl) / (rh - rl), 0)   # guard flat g' (g'=0) ties
    eta <- el + t * (eh - el)
    eta[undef] <- NA_real_                          # NOT clamped to 1
    out[mix] <- eta
  }
  out
}
