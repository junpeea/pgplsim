############################################################
# Compare estimated link functions from Lu GPLSIM and Yu GPLSIM
# on a common single-index axis
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(splines)
})

# ----------------------------------------------------------
# 1. Choose which outcome/model family to compare
#    Use ONE of the following blocks
# ----------------------------------------------------------

## ---- For Binomial outcome ----
fit_glm_cmp    <- fit_logit_glm
fit_gplsim_cmp <- fit_bino_gplsim
family_name    <- "binomial"
title_text     <- "Estimated link functions along Lu's single index (Binomial)"

# ----------------------------------------------------------
# 2. Build common single-index axis using Lu's alpha_hat
# ----------------------------------------------------------
Lu_alpha_hat <- as.numeric(fit_gplsim_cmp$alpha_hat)
Lu_alpha_hat <- Lu_alpha_hat / sqrt(sum(Lu_alpha_hat^2))
Lu_idx       <-  as.numeric(fit_gplsim_cmp$u_hat)
Yu_alpha_hat <- as.numeric(fit_yu$theta)
Yu_idx       <- as.numeric(as.matrix(Z2) %*% Yu_alpha_hat)
idx_hat = c(Lu_idx,Yu_idx)

grid <- seq(
  quantile(idx_hat, 0.01, na.rm = TRUE),
  quantile(idx_hat, 0.99, na.rm = TRUE),
  length.out = 200
)

Z_Lu_grid <- outer(grid, Lu_alpha_hat)
Z_Yu_grid <- outer(grid, Yu_alpha_hat)

colnames(Z_Lu_grid) <- c("BMI_s", "BPSysAve_s", "DirectChol_s")
colnames(Z_Yu_grid) <- c("BMI_s", "BPSysAve_s", "DirectChol_s")

# Reference profile for factor covariates
ref_gender <- levels(dat1$Gender)[1]
ref_race   <- levels(dat1$Race_new)[1]

newdat_Lu <- data.frame(
  Age_s         = rep(0,nrow(Z_Lu_grid)),
  BMI_s         = Z_Lu_grid[, "BMI_s"],
  BPSysAve_s    = Z_Lu_grid[, "BPSysAve_s"],
  DirectChol_s  = Z_Lu_grid[, "DirectChol_s"],
  Gender        = factor(ref_gender, levels = levels(dat1$Gender)),
  Race_new      = factor(ref_race,   levels = levels(dat1$Race_new))
)

X_ref <- matrix(
  0,
  nrow = length(grid),
  ncol = ncol(X2)
)
colnames(X_ref) <- colnames(X2)

newdat_Yu <- data.frame(
  a = grid,
  x = I(as.matrix(Z_Yu_grid)),
  z = I(as.matrix(X_ref)),
  Age_s         = rep(0,nrow(Z_Yu_grid)),
  BMI_s         = Z_Yu_grid[, "BMI_s"],
  BPSysAve_s    = Z_Yu_grid[, "BPSysAve_s"],
  DirectChol_s  = Z_Yu_grid[, "DirectChol_s"],
  Gender        = factor(ref_gender, levels = levels(dat1$Gender)),
  Race_new      = factor(ref_race,   levels = levels(dat1$Race_new))
)

# ----------------------------------------------------------
# 3. GLM curve on link scale with 95% CI
# ----------------------------------------------------------
pr_Yu_glm <- predict(fit_glm_cmp, newdata = newdat_Yu, type = "link", se.fit = TRUE)
pr_Lu_glm <- predict(fit_glm_cmp, newdata = newdat_Lu, type = "link", se.fit = TRUE)

df_Yu_glm <- data.frame(
  s = grid,
  fit = as.numeric(pr_Yu_glm$fit),
  lo  = as.numeric(pr_Yu_glm$fit - 1.96 * pr_Yu_glm$se.fit),
  hi  = as.numeric(pr_Yu_glm$fit + 1.96 * pr_Yu_glm$se.fit),
  method = "GLM"
)
df_Lu_glm <- data.frame(
  s = grid,
  fit = as.numeric(pr_Lu_glm$fit),
  lo  = as.numeric(pr_Lu_glm$fit - 1.96 * pr_Lu_glm$se.fit),
  hi  = as.numeric(pr_Lu_glm$fit + 1.96 * pr_Lu_glm$se.fit),
  method = "GLM"
)

# ----------------------------------------------------------
# 5. Yu GPLSIM curve with pointwise CI
# ----------------------------------------------------------
pred_yu_link <- predict(
  fit_yu,
  newdata = newdat_Yu,
  type = "link",
  se.fit = TRUE
)

df_yu_gplsim = data.frame(
  s = grid,
  eta   = as.numeric(pred_yu_link$fit),
  se    = as.numeric(pred_yu_link$se.fit)
) %>%
  mutate(
    fit   = eta,
    lo = eta - 1.96 * se,
    hi = eta + 1.96 * se,
    method = "Yu GPLSIM"
  ) %>% dplyr::select(s,fit,lo,hi,method)

# ----------------------------------------------------------
# 6. Lu GPLSIM curve with pointwise CI
# ----------------------------------------------------------
# Assumes the fitted GPLSIM object contains:
#   u_hat      : observed index values
#   phi_hat    : estimated link function at u_hat
#   BB_hat     : spline basis matrix evaluated at u_hat
#   gamma_hat  : spline coefficients
#   vcov_model / vcov_robust / vcov_levadj / vcov_finite
#
# For your current implementation, the last length(gamma_hat) rows/cols
# of the vcov matrix correspond to the spline coefficients.

extract_gplsim_curve <- function(fit_obj, s_grid = NULL,
                                 vcov_type = c("model", "robust", "levadj", "finite"),
                                 smooth_curve = TRUE) {
  vcov_type <- match.arg(vcov_type)
  
  needed <- c("u_hat", "phi_hat", "BB_hat", "gamma_hat")
  stopifnot(all(needed %in% names(fit_obj)))
  
  V <- switch(
    vcov_type,
    model  = fit_obj$vcov_model,
    robust = fit_obj$vcov_robust,
    levadj = fit_obj$vcov_levadj,
    finite = fit_obj$vcov_finite
  )
  
  # spline-coefficient covariance block
  k <- length(fit_obj$gamma_hat)
  Vg <- V[(nrow(V) - k + 1):nrow(V), (ncol(V) - k + 1):ncol(V), drop = FALSE]
  
  u   <- as.numeric(fit_obj$u_hat)
  phi <- as.numeric(fit_obj$phi_hat)
  B   <- as.matrix(fit_obj$BB_hat)
  
  ok <- is.finite(u) & is.finite(phi) & apply(B, 1, function(x) all(is.finite(x)))
  u   <- u[ok]
  phi <- phi[ok]
  B   <- B[ok, , drop = FALSE]
  
  # pointwise delta-method SE for phi(u)
  var_phi <- rowSums((B %*% Vg) * B)
  se_phi  <- sqrt(pmax(var_phi, 0))
  
  df_obs <- data.frame(
    s   = u,
    fit = phi,
    lo  = phi - 1.96 * se_phi,
    hi  = phi + 1.96 * se_phi
  )
  
  # sort by index
  df_obs <- df_obs[order(df_obs$s), ]
  
  # collapse duplicate s values if needed
  df_obs <- df_obs |>
    dplyr::group_by(s) |>
    dplyr::summarise(
      fit = mean(fit),
      lo  = mean(lo),
      hi  = mean(hi),
      .groups = "drop"
    )
  
  # return observed-scale curve if no grid requested
  if (is.null(s_grid)) {
    return(df_obs)
  }
  
  if (!smooth_curve) {
    # simple interpolation without spline smoothing
    fit_interp <- approx(df_obs$s, df_obs$fit, xout = s_grid, rule = 2)$y
    lo_interp  <- approx(df_obs$s, df_obs$lo,  xout = s_grid, rule = 2)$y
    hi_interp  <- approx(df_obs$s, df_obs$hi,  xout = s_grid, rule = 2)$y
    
    return(data.frame(
      s   = s_grid,
      fit = fit_interp,
      lo  = lo_interp,
      hi  = hi_interp
    ))
  }
  
  # spline smoothing for plotting
  sm_fit <- smooth.spline(df_obs$s, df_obs$fit)
  sm_lo  <- smooth.spline(df_obs$s, df_obs$lo)
  sm_hi  <- smooth.spline(df_obs$s, df_obs$hi)
  
  data.frame(
    s   = s_grid,
    fit = predict(sm_fit, s_grid)$y,
    lo  = predict(sm_lo,  s_grid)$y,
    hi  = predict(sm_hi,  s_grid)$y
  )
}

df_lu_gplsim <- extract_gplsim_curve(
  fit_obj    = fit_gplsim_cmp,
  s_grid     = grid,
  vcov_type  = "model",   # change to "robust", "levadj", or "finite" if desired
  smooth_curve = TRUE
) |>
  dplyr::mutate(method = "Lu GPLSIM")

# ----------------------------------------------------------
# 6. Combine and plot
# ----------------------------------------------------------

library(dplyr)
library(ggplot2)

# Mean-correct each method-specific curve
df_Yu_glm <- df_Yu_glm %>%
  mutate(
    mean_fit = mean(fit, na.rm = TRUE),
    fit = fit - mean_fit,
    lo  = ifelse(is.finite(lo), lo - mean_fit, lo),
    hi  = ifelse(is.finite(hi), hi - mean_fit, hi)
  )
df_Lu_glm <- df_Lu_glm %>%
  mutate(
    mean_fit = mean(fit, na.rm = TRUE),
    fit = fit - mean_fit,
    lo  = ifelse(is.finite(lo), lo - mean_fit, lo),
    hi  = ifelse(is.finite(hi), hi - mean_fit, hi)
  )

# Create panel-specific data
plot_left <- bind_rows(
  df_Yu_glm %>% mutate(panel = "GLM vs Yu GPLSIM"),
  df_yu_gplsim %>% mutate(panel = "GLM vs Yu GPLSIM")
)

plot_right <- bind_rows(
  df_Lu_glm %>% mutate(panel = "GLM vs Lu GPLSIM"),
  df_lu_gplsim %>% mutate(panel = "GLM vs Lu GPLSIM")
)

plot_df <- bind_rows(plot_left, plot_right)

# common y-limits across both panels
ymin_all <- min(plot_df$fit, plot_df$lo, na.rm = TRUE)
ymax_all <- max(plot_df$fit, plot_df$hi, na.rm = TRUE)

p_link_compare <- ggplot() +
  
  # ---------------- Left panel: Yu GPLSIM CI ----------------
geom_line(
  data = subset(plot_df, panel == "GLM vs Yu GPLSIM" & method == "Yu GPLSIM" & is.finite(lo)),
  aes(x = s, y = lo, color = method),
  linewidth = 0.8,
  linetype = "dashed"
) +
  geom_line(
    data = subset(plot_df, panel == "GLM vs Yu GPLSIM" & method == "Yu GPLSIM" & is.finite(hi)),
    aes(x = s, y = hi, color = method),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # ---------------- Right panel: Lu GPLSIM CI ----------------
geom_line(
  data = subset(plot_df, panel == "GLM vs Lu GPLSIM" & method == "Lu GPLSIM" & is.finite(lo)),
  aes(x = s, y = lo, color = method),
  linewidth = 0.8,
  linetype = "dashed"
) +
  geom_line(
    data = subset(plot_df, panel == "GLM vs Lu GPLSIM" & method == "Lu GPLSIM" & is.finite(hi)),
    aes(x = s, y = hi, color = method),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # ---------------- Main fitted curves ----------------
geom_line(
  data = subset(plot_df, method == "GLM"),
  aes(x = s, y = fit, color = method),
  linewidth = 1,
  linetype = "solid"
) +
  geom_line(
    data = subset(plot_df, method == "Yu GPLSIM"),
    aes(x = s, y = fit, color = method),
    linewidth = 1,
    linetype = "solid"
  ) +
  geom_line(
    data = subset(plot_df, method == "Lu GPLSIM"),
    aes(x = s, y = fit, color = method),
    linewidth = 1.1,
    linetype = "solid"
  ) +
  
  facet_wrap(~ panel, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c(
    "GLM" = "black",
    "Yu GPLSIM" = "red",
    "Lu GPLSIM" = "blue"
  )) +
  coord_cartesian(ylim = c(ymin_all, ymax_all)) +
  labs(
    title = "",
    x = "single index",
    y = "estimated link function"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_link_compare)


