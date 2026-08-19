# ============================================================================ #
# Plots.R — figures and log-likelihood comparison
# ============================================================================ #
# Builds the four result figures of the paper from the FITTED model state. Each
# is a self-contained function with a common (show, save, filename, dir)
# signature: show prints to the active device, save writes a PDF to `dir`, and
# the plot object is returned invisibly. Sourcing this file only DEFINES the
# functions; nothing is drawn until they are called.
#   plot_weight_functions  single-exposure weights lambda_1, lambda_2 and the
#                          Pb->BaP dose alignment (base R, 1x3).
#   plot_weight_surfaces   joint weight surfaces (2x2): the three additive nulls
#                          (Loewe+, Hand, Bliss) and the radial fit.
#   plot_cdf_scatter       fitted vs empirical CDFs, all data and joint only.
#   plot_loglik            log-likelihood across data subsets (all / no-(0,0) /
#                          joint) and models (H0 nulls, H1 fit +/- penalty).
# Helpers: contour_plot_01 (single surface panel), l_fit_rays_fast (pointwise
# radial weight), CDFs_scatter_plots (single scatter panel), and the lab_fun /
# lab_fun_logL plotmath label maps. ggplot panels are composed with patchwork.
# ============================================================================ #
# Requires: Utils.R, Loglikelihood.R, Data_setup.R, Null_models.R, and
# Model_fit.R sourced (the last provides the fitted state: fit_w_BaP/fit_w_Pb,
# fit_tau2, l1/g/l2, par_fitJ, data_list/_no0/_core, pen_cache, ga_pen; plus the
# null models l_mean/l_bliss/l_hand and the data subsets x_*_tri / y_*_tri).
# Writes PDFs to `fig_dir` (set in main.R) when save = TRUE.
# ============================================================================ #

# Single exposure — weight functions -------------------------------------------

## `plot_weight_surfaces`
## 
## Compose the four joint weight-surface panels into a 2x2 figure: the three
## additive-null surfaces (Loewe+, Hand, Bliss) and the fitted radial surface,
## each from contour_plot_01, sharing one lambda_12 colorbar (guides = "collect").
## 
## 
## Args:
##   show, save : print to the device / write a PDF.
##   filename, dir, width, height : output file and page size.
## Returns:
##   The patchwork object (invisibly).
## State:
##   Reads the null models l_mean / l_hand / l_bliss and par_fitJ (via fit_fun);
##   contour_plot_01 additionally reads d1_MAX / d2_MAX.
##
plot_weight_functions <- function(show = TRUE, save = FALSE,
                                  filename = "1D_weight_functions", dir = fig_dir,
                                  width = 12, height = 4) {
  draw <- function() {
    op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))
    on.exit(par(op))
    
    ## BaP weight function
    x.grid.BaP <- seq(0, dInf_BaP, by = 0.01)
    matplot(x = x.grid.BaP, y = basis_BaP(x.grid.BaP),
            type = 'l', lty = 2, col = 1 + 1:nI_BaP,
            ylim = c(0, 1), xlab = expression("BaP dose [" * mu * "M]"),
            ylab = expression(lambda[1](d[1])), main = 'BaP Weight Function')
    lines(x.grid.BaP, y = basis_BaP(x.grid.BaP) %*% fit_w_BaP, lwd = 2)
    points(xBaP, y = basis_BaP(xBaP) %*% fit_w_BaP, lwd = 2, pch = 1)
    
    ## Pb weight function
    x.grid.Pb <- seq(0, dInf_Pb, by = 0.01)
    matplot(x = x.grid.Pb, y = basis_BaP(dInf_BaP * fit_tau2 * basis_Pb(x.grid.Pb) %*% fit_w_Pb),
            type = 'l', lty = 2, col = 1 + 1:nI_Pb,
            ylim = c(0, 1), xlab = expression("Pb dose [" * mu * "M]"),
            ylab = expression(lambda[2](d[2])), main = 'Pb Weight Function')
    lines(x.grid.Pb, y = basis_BaP(dInf_BaP * fit_tau2 * basis_Pb(x.grid.Pb) %*% fit_w_Pb) %*% fit_w_BaP, lwd = 2)
    points(xPb, y = basis_BaP(dInf_BaP * fit_tau2 * basis_Pb(xPb) %*% fit_w_Pb) %*% fit_w_BaP, lwd = 2, pch = 1)
    
    ## BaP-equivalent dose (alignment)
    x.grid.Pb <- seq(0, dInf_Pb, by = 0.01)
    matplot(x = x.grid.Pb, y = dInf_BaP * fit_tau2 * basis_Pb(x.grid.Pb),
            type = 'l', lty = 2, col = 1 + 1:nI_Pb,
            ylim = c(0, dInf_BaP), xlab = expression("Pb dose [" * mu * "M]"), 
            ylab = expression(D[1]^{(max)} * tau[2] * g(d[2])), main = 'BaP-equivalent Dose')
    lines(x.grid.Pb, y = dInf_BaP * fit_tau2 * basis_Pb(x.grid.Pb) %*% fit_w_Pb, lwd = 2)
    points(xPb, y = dInf_BaP * fit_tau2 * basis_Pb(xPb) %*% fit_w_Pb, lwd = 2, pch = 1)
  }
  
  if (show) draw()
  if (save) { pdf(file.path(dir, paste0(filename, ".pdf")), width = width, height = height); draw(); dev.off() }
  invisible(NULL)
}

# Joint exposure — weight-surface contour plots --------------------------------

## `contour_plot_01`
##
## Evaluate a weight function on a dense (d1, d2) grid and render it as a filled
## raster with white iso-effect contours at 0.1 spacing. Optionally trims to the
## probed lower triangle. Panel title and colorbar name are passed in so the same
## helper serves both the null surfaces and the fit.
##
## Args:
##   l_func    : weight function lambda_12(d1, d2), evaluated per grid point.
##   leg_title : colorbar name (an expression, e.g. lambda[12]).
##   trim      : if TRUE, keep only d2 <= d2_MAX - (d2_MAX/d1_MAX) d1.
##   title     : panel heading (an expression).
## Returns:
##   A ggplot object.
## State:
##   Reads d1_MAX, d2_MAX.
##
contour_plot_01 <- function(l_func, leg_title = 'l12', trim = FALSE, title='Weight Surface') {
  
  d1_grid <- seq(0, d1_MAX, by = 0.05)
  d2_grid <- seq(0, d2_MAX, by = 1)
  grid_l  <- expand.grid(d1 = d1_grid, d2 = d2_grid)
  
  grid_l <- subset(grid_l, d1 + d2 > 0)
  
  grid_l$z <- mapply(l_func, grid_l$d1, grid_l$d2)
  grid_l   <- subset(grid_l, is.finite(z))
  
  if (trim) {
    grid_l <- subset(grid_l, d2 <= d2_MAX - (d2_MAX / d1_MAX) * d1)
  }
  
  ggplot(grid_l, aes(d1, d2, fill = z)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = z), color = "white",
                 breaks = seq(0, 1, by = 0.1), linewidth = 0.5) +
    scale_fill_viridis_c(name = leg_title, limits = c(0, 1)) +
    labs(x = expression("BaP dose [" * mu * "M]"),
         y = expression("Pb dose [" * mu * "M]"), 
         title = title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

##  `l_fit_rays_fast`
## 
##  Pointwise fitted radial weight lambda_12(d1, d2) for a given parameter matrix.
##  Scalar counterpart of ray_weights: aligns the dose pair to (rho, nu), forms the
##  softmax mixing weights over the Phi basis, and pushes the effective dose index
##  through the reference weight l1. Used by plot_weight_surfaces to render the fit.
## 
##  Args:
##    d1, d2  : single BaP and Pb dose (scalars).
##    par_mat : (K_phi-1) x K_delta radial-spline parameter matrix.
##  Returns:
##    Weight in [0,1]; 0 at the origin, NA where rho <= 0 or nu is non-finite.
##  State:
##    Reads dInf_BaP, fit_tau2, g, base_Delta, base_Phi, R_MAX, l1.
##
l_fit_rays_fast <- function(d1, d2, par_mat) {
  
  if (d1 + d2 == 0) return(0)
  
  x   <- d1
  y   <- dInf_BaP * fit_tau2 * g(d2)
  rho <- x + y
  nu  <- y / rho
  
  # handle rho == 0 safely
  if (!is.finite(nu) || rho <= 0) return(NA_real_)
  
  Delta <- base_Delta(nu)                  # 1 x K_delta
  Phi   <- base_Phi(rho)                   # 1 x K_phi
  
  w_0     <- as.vector(Delta %*% t(par_mat))     # length K_phi - 1
  logit_0 <- softmax(c(w_0, 0))
  phi     <- as.vector(Phi %*% logit_0)
  
  l1(R_MAX * phi)
}

##  `plot_weight_surfaces`
## 
##  Compose the four joint weight-surface panels into a 2x2 figure: the three
##  additive-null surfaces (Loewe+, Hand, Bliss) and the fitted radial surface,
##  each from contour_plot_01, sharing one lambda_12 colorbar (guides = "collect").
## 
##  Args:
##    show, save : print to the device / write a PDF.
##    filename, dir, width, height : output file and page size.
##  Returns:
##    The patchwork object (invisibly).
##  State:
##    Reads the null models l_mean / l_hand / l_bliss and par_fitJ (via fit_fun);
##    contour_plot_01 additionally reads d1_MAX / d2_MAX.
##  
plot_weight_surfaces <- function(show = TRUE, save = FALSE,
                                 filename = "2D_weight_surfaces", dir = fig_dir,
                                 width = 9, height = 8) {
  leg_ttl <- expression(lambda[12])
  fit_fun <- function(u1, u2) l_fit_rays_fast(u1, u2, par_mat = par_fitJ)
  
  p <- contour_plot_01(l_mean,  leg_title = leg_ttl, trim = TRUE,
                       title=expression(lambda[12]^(o) ~ bold("Loewe"^"+"))) +
    contour_plot_01(l_hand,  leg_title = leg_ttl, trim = TRUE,
                    title=expression(lambda[12]^(o) ~ bold("Bliss"))) +
    contour_plot_01(l_bliss, leg_title = leg_ttl, trim = TRUE,
                    title=expression(lambda[12]^(o) ~ bold("Hand"))) +
    contour_plot_01(fit_fun, leg_title = leg_ttl, trim = TRUE,
                    title=expression(lambda[12] ~ bold("Radial Fit"))) +
    plot_layout(ncol = 2, guides = "collect") 
  
  if (show) print(p)
  if (save) ggsave(file.path(dir, paste0(filename, ".pdf")), p, width = width, height = height)
  invisible(p)
}

# Model fit comparison — log-likelihood ----------------------------------------

##  `lab_fun`
## 
##  Plotmath legend labeller for the CDF scatter: maps the model keys
##  (Mean / Hand / Bliss / RayFit) to expressions like H[0] "Loewe"^"+". Passed as
##  `labels` to the color/fill/shape scales in CDFs_scatter_plots.
## 
##  Args:
##    x : character (or factor) vector of model keys.
##  Returns:
##    An expression vector, one entry per element of x.
##
lab_fun <- function(x) {
  map <- list(
    Mean   = expression(H[0] * " Loewe"^"+"),
    Hand   = expression(H[0] * " Hand"),
    Bliss  = expression(H[0] * " Bliss"),
    RayFit = expression(H[1] * " Radial Fit")
  )
  do.call(expression, lapply(as.character(x), function(k) map[[k]][[1]]))
}

##  `CDFs_scatter_plots`
## 
##  One fitted-vs-empirical CDF scatter panel: points colored/shaped by model, the
##  y = x diagonal for reference, a "CDFs" annotation, and plotmath legend labels.
##  Called twice by plot_cdf_scatter (all data, joint only) and combined there.
## 
##  Args:
##    df    : long data frame with columns Actual, Predicted, Model
##            (Model a factor whose level order sets the legend order).
##    title : panel heading.
##  Returns:
##    A ggplot object.
##  State:
##    None (uses lab_fun; all data passed in).
## 
CDFs_scatter_plots <- function(df, title=NULL) {
  cols_meth <- c("Hand" = "#57B9B2", "Mean" = "#DC267F",
                 "Bliss" = "#FFB000", "RayFit" = "#648FFF")
  
  ggplot(df, aes(x = Actual, y = Predicted,
                 color = Model, fill = Model, shape = Model)) +
    geom_point(alpha = 0.6, size = 1.5, stroke = 0.4) +
    theme_minimal() +
    annotate("text", x = 0.05, y = 0.975, label = "CDFs",
             hjust = 0, vjust = 1, size = 4) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2, color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_manual(values = cols_meth, labels = lab_fun, name = expression(lambda[12])) +
    scale_fill_manual (values = cols_meth, labels = lab_fun, name = expression(lambda[12])) +
    scale_shape_manual(values = c("Hand" = 25, "Mean" = 16,"Bliss" = 15, "RayFit" = 17),
                       labels = lab_fun, name = expression(lambda[12])) + 
    labs(x = "Empirical", y = "Fitted", title = title) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 3)),
           fill  = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3)))
}

##  `plot_cdf_scatter`
## 
##  Assemble the two-panel CDF figure. Computes model weights at the measured joint
##  doses (nulls via l_mean/l_bliss/l_hand, fit via ray_weights), turns them into
##  fitted CDFs (get_probs), reshapes to long form, and pairs the all-data and
##  joint-only panels with a single shared bottom legend.
## 
##  Args:
##    show, save : print / write PDF.
##    filename, dir, width, height : output file and page size.
##  Returns:
##    The patchwork object (invisibly).
##  State:
##    Reads x_all_tri / F_all_tri, x_core_tri? (F_core_tri, idx_core), F0, Finf,
##    par_fitJ, data_list, and the null models l_mean / l_bliss / l_hand.
## 
plot_cdf_scatter <- function(show = TRUE, save = FALSE,
                             filename = "CDF_scatter", dir = fig_dir,
                             width = 10, height = 5) {
  h_mean  <- apply(x_all_tri, 1, function(v) l_mean(v[1], v[2]))
  h_bliss <- apply(x_all_tri, 1, function(v) l_bliss(v[1], v[2]))
  h_hand  <- apply(x_all_tri, 1, function(v) l_hand(v[1], v[2]))
  h_fit   <- ray_weights(par_fitJ, data_list)
  
  F_mean  <- get_probs(h_mean,  F0, Finf)
  F_bliss <- get_probs(h_bliss, F0, Finf)
  F_hand  <- get_probs(h_hand,  F0, Finf)
  F_fit   <- get_probs(h_fit,   F0, Finf)
  
  lvl <- c("Mean", "Hand", "Bliss", "RayFit")   # desired legend order
  
  df_all <- data.frame(Actual = as.vector(F_all_tri),
                       Hand = as.vector(F_hand), Mean = as.vector(F_mean),
                       Bliss = as.vector(F_bliss), RayFit = as.vector(F_fit)) %>%
    pivot_longer(cols = -Actual, names_to = "Model", values_to = "Predicted")
  df_all$Model <- factor(df_all$Model, levels = lvl)
  
  df_core <- data.frame(Actual = as.vector(F_core_tri),
                        Hand = as.vector(F_hand[, idx_core]), Mean = as.vector(F_mean[, idx_core]),
                        Bliss = as.vector(F_bliss[, idx_core]), RayFit = as.vector(F_fit[, idx_core])) %>%
    pivot_longer(cols = -Actual, names_to = "Model", values_to = "Predicted") 
  df_core$Model <- factor(df_core$Model, levels = lvl)
  
  p <- CDFs_scatter_plots(df_all, title='All data') + CDFs_scatter_plots(df_core, title='Joint exposures only') +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if (show) print(p)
  if (save) ggsave(file.path(dir, paste0(filename, ".pdf")), p, width = width, height = height)
  invisible(p)
}

# Model fit comparison — CDFs scatter plot -------------------------------------

##  `lab_fun_logL`
## 
##  Plotmath legend labeller for the log-likelihood figure: same idea as lab_fun
##  but keyed to the five method levels (mean / hand / bliss / fit / fit_pen),
##  adding the H[1] "Penalized" entry.
## 
##  Args:
##    x : character (or factor) vector of method keys.
##  Returns:
##    An expression vector, one entry per element of x.
## 
lab_fun_logL <- function(x) {
  map <- list(
    mean    = expression(H[0] * " Loewe"^"+"),
    hand    = expression(H[0] * " Hand"),
    bliss   = expression(H[0] * " Bliss"),
    fit     = expression(H[1] * " Radial Fit"),
    fit_pen = expression(H[1] * " Penalized")
  )
  do.call(expression, lapply(as.character(x), function(k) map[[k]][[1]]))
}

##  `plot_loglik`
## 
##  Fit-quality figure: log-likelihood of each model on each data subset. Scores
##  the radial fit (penalized and unpenalized, via loglik_joint) and the three
##  nulls (via loglik_joint_fixed) on the all / no-(0,0) / joint subsets, then
##  plots them faceted by subset with free y-scales.
## 
##  Args:
##    show, save : print / write PDF.
##    filename, dir, width, height : output file and page size.
##  Returns:
##    The ggplot object (invisibly).
##  State:
##    Reads par_fitJ, data_list / data_no0 / data_core, pen_cache, ga_pen, D0, Dinf,
##    the null models l_mean / l_bliss / l_hand, and the data subsets x_*/y_*.
## 
plot_loglik <- function(show = TRUE, save = FALSE,
                        filename = "loglikelihood_comparison", dir = fig_dir,
                        width = 8, height = 3.5) {
  logL_vals <- data.frame(
    data = rep(c("all", "no00", "core"), each = 5),
    meth = rep(c("fit_pen", "fit", "mean", "bliss", "hand"), 3),
    logL = c(
      loglik_joint(par_fitJ, data_list, pen_cache, 0, ga_pen, D0, Dinf),
      loglik_joint(par_fitJ, data_list, pen_cache, 0, 0,      D0, Dinf),
      loglik_joint_fixed(l_mean,  x_all_tri, y_all_tri, D0, Dinf),
      loglik_joint_fixed(l_bliss, x_all_tri, y_all_tri, D0, Dinf),
      loglik_joint_fixed(l_hand,  x_all_tri, y_all_tri, D0, Dinf),
      loglik_joint(par_fitJ, data_no0, pen_cache, 0, ga_pen, D0, Dinf),
      loglik_joint(par_fitJ, data_no0, pen_cache, 0, 0,      D0, Dinf),
      loglik_joint_fixed(l_mean,  x_tri_no0, y_tri_no0, D0, Dinf),
      loglik_joint_fixed(l_bliss, x_tri_no0, y_tri_no0, D0, Dinf),
      loglik_joint_fixed(l_hand,  x_tri_no0, y_tri_no0, D0, Dinf),
      loglik_joint(par_fitJ, data_core, pen_cache, 0, ga_pen, D0, Dinf),
      loglik_joint(par_fitJ, data_core, pen_cache, 0, 0,      D0, Dinf),
      loglik_joint_fixed(l_mean,  x_core_tri, y_core_tri, D0, Dinf),
      loglik_joint_fixed(l_bliss, x_core_tri, y_core_tri, D0, Dinf),
      loglik_joint_fixed(l_hand,  x_core_tri, y_core_tri, D0, Dinf)
    )
  )
  logL_vals$meth <- factor(logL_vals$meth, levels = c("mean", "hand", "bliss", "fit", "fit_pen"))
  logL_vals$data <- factor(logL_vals$data, levels = c("all", "no00", "core"))
  
  labs_meth <- c(fit_pen = "H1 - Fit+Penalty", fit = "H1 - Fit",
                 mean = "H0 - Mean", bliss = "H0 - Bliss", hand = "H0 - Hand")
  labs_data <- c(all = "All Data", no00 = "All Data - No (0,0)", core = "Only Joint Exposure")
  cols_meth <- c(hand = "#57B9B2", mean = "#DC267F", bliss = "#FFB000",
                 fit = "#648FFF", fit_pen = "#785EF0")
  pch_meth  <- c(hand = 25, mean = 16, bliss = 15, fit = 17, fit_pen = 18)
  
  p <- ggplot(logL_vals, aes(meth, logL, color = meth, shape = meth, fill = meth)) +
    geom_point(size = 3) + theme_bw() +
    facet_wrap(~data, scales = "free_y", nrow = 1, labeller = as_labeller(labs_data)) +
    scale_color_manual(values = cols_meth, labels = lab_fun_logL, name = expression(lambda[12])) +
    scale_fill_manual(values = cols_meth, labels = lab_fun_logL, name = expression(lambda[12])) +
    scale_shape_manual(values = pch_meth, labels = lab_fun_logL, name = expression(lambda[12])) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2, "pt"), legend.key.spacing.x = unit(10, "pt")) +
    labs(x = NULL, y = "Log-Likelihood", color = "Model", fill = "Model", shape = "Model") +
    guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  
  if (show) print(p)
  if (save) ggsave(file.path(dir, paste0(filename, ".pdf")), p, width = width, height = height)
  invisible(p)
}