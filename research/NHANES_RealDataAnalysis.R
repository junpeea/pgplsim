############################################################
# Real Data Analysis for Section 5.2
# Dataset assumed already created as:
#
#
# Outcomes:
# Diabetes        : binary outcome
#
# Analysis plan:
#   - Construct X (linear part) and Z (index-related variables)
#   - Fit standard GLMs as baseline models
#   - Fit GAMs as flexible benchmarks for nonlinear effects
#   - Summarize coefficients, model fit, and visualize nonlinear patterns
#
# Notes:
#   - This script does NOT implement your custom GPLSIM estimator directly.
#   - Instead, it prepares the data and provides strong baseline/flexible
#     analyses that can accompany or benchmark your spline-based method.
############################################################

suppressPackageStartupMessages({
  libs <- c("NHANES","dplyr", "ggplot2", "mgcv", "broom", "patchwork","mgcv")
  to_install <- setdiff(libs, rownames(installed.packages()))
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(libs, require, character.only = TRUE))
})

options(stringsAsFactors = FALSE)

############################################################
# 1. Data preparation
############################################################
# Data
data(NHANES)
df <- NHANES
dat <- df %>%
  dplyr::select(
    Diabetes,        # outcome (Yes/No)
    Age,             # numeric
    BMI,             # numeric
    BPSysAve,        # systolic BP
    # TotChol,         # Total cholesterol
    DirectChol,      # HDL cholesterol
    Gender,          # factor
    Race1            # factor
  ) %>%
  na.omit()

# Defensive copy
dat0 <- dat

# Recode factors cleanly
dat1 <- dat0 %>%
  mutate(
    Gender   = factor(Gender),
    Race1    = factor(Race1),
    Diabetes = factor(Diabetes, levels = c("No", "Yes")),
    Race_new = as.character(Race1),
    Race_new = ifelse(Race_new %in% c("Mexican", "Hispanic"),
                      "Hispanic",
                      Race_new),
    Race_new = factor(
      Race_new,
      levels = c("White", "Black", "Hispanic", "Other")
    )
  ) %>%
  filter(
    !is.na(Diabetes),
    !is.na(Age),
    !is.na(BMI),
    !is.na(BPSysAve),
    !is.na(DirectChol),
    !is.na(Gender),
    !is.na(Race_new)
  )
table(dat1$Race_new)
levels(dat1$Race_new)

# Binary outcome
dat1 <- dat1 %>%
  mutate(
    Diabetes_bin = ifelse(Diabetes == "Yes", 1L, 0L)
  )

# Standardized continuous predictors for stable estimation
dat1 <- dat1 %>%
  mutate(
    Age_s     = as.numeric(scale(Age)),
    BMI_s     = as.numeric(scale(BMI)),
    BPSysAve_s= as.numeric(scale(BPSysAve)),
    DirectChol_s = as.numeric(scale(DirectChol))
  )

# Index-style matrix Z and linear-design matrix X
Z <- as.matrix(dat1[, c("BMI_s", "BPSysAve_s", "DirectChol_s")])
X <- model.matrix(~ Age_s + Gender + Race_new, data = dat1)

cat("Sample size:", nrow(dat1), "\n")
cat("Diabetes prevalence:", mean(dat1$Diabetes_bin), "\n\n")
summary(dat1)

############################################################
# 2. Descriptive summaries
############################################################

cat("===== Descriptive Summary =====\n")
print(summary(dat1[, c("Diabetes_bin", "Age", "BMI", "BPSysAve", "DirectChol")]))
cat("\nCounts by Gender:\n")
print(table(dat1$Gender))
cat("\nCounts by Race_new:\n")
print(table(dat1$Race_new))
cat("\n")

############################################################
# 3. Baseline Logistic GLM for Diabetes
############################################################

fit_logit_glm <- glm(
  Diabetes_bin ~ Age_s + BMI_s + BPSysAve_s + DirectChol_s + Gender + Race_new,
  family = binomial(link = "logit"),
  data = dat1
)

cat("===== Logistic GLM Summary =====\n")
print(summary(fit_logit_glm))
cat("\n")

############################################################
# 4. Flexible Logistic GAM for nonlinear effects
############################################################

fit_logit_gam <- mgcv::gam(
  Diabetes_bin ~
    s(BMI_s, k = 10) +
    s(BPSysAve_s, k = 10) +
    s(DirectChol_s, k = 10) +
    Age_s + Gender + Race_new,
  family = binomial(link = "logit"),
  method = "REML",
  data = dat1
)

cat("===== Logistic GAM Summary =====\n")
print(summary(fit_logit_gam))
coef(fit_logit_gam)
cat("\n")

############################################################
# 5. Proposed GPLSIM for Logistic
############################################################
source("Updates2026/sim_lib_260305.R")
source("Updates2026/Binomial_Lu_GPLSIM.R")
X2 <- X
Z2 <- Z

is_binary <- function(v) {
  ux <- unique(na.omit(v))
  length(ux) <= 2 && all(sort(ux) %in% c(0, 1))
}

for (j in seq_len(ncol(X2))) {
  if (!is_binary(X2[, j])) {
    X2[, j] <- scale(X2[, j])
  }
}

for (j in seq_len(ncol(Z2))) {
  if (!is_binary(Z2[, j])) {
    Z2[, j] <- scale(Z2[, j])
  }
}

fit_bino_gplsim <- Binomial_Lu_GPLSIM(
  Y = dat1$Diabetes_bin,
  X = X2,
  Z = Z2,
  M = 80,
  err = 1e-6,
  lambda_init = 1,
  make_plot = TRUE,
  verbose = TRUE
)
fit_bino_gplsim$coef_table
fit_bino_gplsim$alpha_hat
fit_bino_gplsim$beta_hat 


############################################################
# 6. Publication-style coefficient tables
############################################################

# Logistic GLM: exponentiated coefficients = odds ratios
logit_tab <- broom::tidy(fit_logit_glm, conf.int = TRUE) %>%
  mutate(
    estimate_or = exp(estimate),
    conf.low_or = exp(conf.low),
    conf.high_or = exp(conf.high)
  )

cat("===== Logistic Regression Table (Odds Ratios) =====\n")
print(logit_tab)
cat("\n")

# source("Updates2026/NHANES_Figure3.R")
source("Updates2026/NHANES_Phi_Figure.R")

# Main reports
print(summary(fit_logit_glm))
print(summary(fit_logit_gam))
fit_bino_gplsim$coef_table
colnames(X2)
