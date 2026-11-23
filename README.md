# pgplsim

<!-- badges: start -->
[![R-CMD-check](https://github.com/junpeea/pgplsim/actions/workflows/rhub.yaml/badge.svg)](https://github.com/junpeea/pgplsim/actions/workflows/rhub.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

## Overview

**pgplsim** is an R package for fitting **Generalized Partially Linear Single-Index Models (GPLSIMs)** using penalized B-splines.

The package provides a unified framework for modeling nonlinear covariate effects through a single-index structure while retaining interpretable linear effects. Automatic smoothing parameter estimation is performed using the Fellner–Schall update, and the package includes model-based, sandwich, and bootstrap variance estimation.

The current release supports:

- Binary outcomes via logistic regression
- Count outcomes via Poisson regression
- Penalized B-spline estimation
- Automatic smoothing parameter selection
- Formula-based model specification
- Prediction for new observations
- Bootstrap inference
- Sandwich covariance estimation
- Publication-quality plotting

---

## Model

`pgplsim` fits models of the form

\[
g\{E(Y \mid X,Z)\}
=
X^\top\beta
+
\phi(Z^\top\alpha),
\]

where

- \(X\) contains covariates entering the linear component,
- \(Z\) contains variables entering the nonlinear single index,
- \(\beta\) contains linear regression coefficients,
- \(\alpha\) is the normalized index coefficient vector, and
- \(\phi(\cdot)\) is an unknown smooth function estimated using penalized splines.

For identifiability, the package normalizes \(\alpha\) so that

\[
\|\alpha\| = 1
\]

and constrains its first component to be positive.

---

## Installation

### Development version

```r
install.packages("remotes")

remotes::install_github("junpeea/pgplsim")
```

or

```r
devtools::install_github("junpeea/pgplsim")
```

---

## Quick example

```r
library(pgplsim)

set.seed(123)

n <- 500

dat <- data.frame(
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.5),
  z1 = rnorm(n),
  z2 = rnorm(n),
  z3 = rnorm(n)
)

alpha <- c(1, 1, 1)
alpha <- alpha / sqrt(sum(alpha^2))

X <- as.matrix(dat[, c("x1", "x2")])
Z <- as.matrix(dat[, c("z1", "z2", "z3")])

eta <- X %*% c(1, -0.5) +
  sin(as.vector(Z %*% alpha))

prob <- plogis(eta)

dat$y <- rbinom(
  n,
  size = 1,
  prob = prob
)

fit <- pgplsim(
  y ~ 0 + x1 + x2,
  index = ~ z1 + z2 + z3,
  data = dat,
  family = binomial()
)

summary(fit)

plot(fit)

pred <- predict(
  fit,
  newdata = dat,
  type = "response"
)

head(pred)
```

---

## Main function

The primary public fitting interface is:

```r
pgplsim(
  formula,
  index,
  data,
  family = binomial(),
  offset = NULL,
  M = NULL,
  lambda = 1,
  control = pgplsim_control()
)
```

For example:

```r
fit <- pgplsim(
  y ~ x1 + x2,
  index = ~ z1 + z2 + z3,
  data = dat,
  family = binomial()
)
```

The arguments have the following roles:

- `formula` specifies the response and linear component.
- `index` specifies variables entering the single-index component.
- `data` supplies the model variables.
- `family` specifies the response distribution.
- `offset` supplies an additive link-scale offset where supported.
- `M` controls the spline basis dimension.
- `lambda` supplies the initial smoothing parameter.
- `control` specifies numerical fitting options.

Computational controls can be specified with:

```r
pgplsim_control(
  maxit = 100,
  tol = 1e-6,
  trace = FALSE
)
```

---

## Supported families

Currently supported:

- `gaussian()` (identity link)
- `binomial()` (logit link)
- `poisson()` (log link)

For Poisson models, `offset` follows the usual GLM convention and is supplied on the **link scale**. For example, with a positive exposure variable:

```r
fit <- pgplsim(
  y ~ x1 + x2,
  index = ~ z1 + z2 + z3,
  data = dat,
  family = poisson(),
  offset = log(dat$exposure)
)
```

Potential future extensions include:

- Survival GPLSIMs

---

## Main features

- Generalized partially linear single-index models
- Penalized B-spline estimation
- Gaussian identity-link GPLSIM
- Logistic GPLSIM
- Poisson GPLSIM
- Formula-based model specification
- Automatic smoothing parameter estimation
- Fisher scoring optimization
- Model-based covariance estimation
- Sandwich covariance estimation
- Bootstrap standard errors
- Prediction using `newdata`
- Publication-quality graphics

---

## Prediction

Prediction follows standard R modeling conventions.

```r
predict(
  fit,
  newdata = dat,
  type = "response"
)
```

Available prediction types include:

```r
type = "response"
type = "link"
type = "index"
type = "smooth"
type = "linear"
```

For example:

```r
u_hat <- predict(
  fit,
  newdata = dat,
  type = "index"
)

smooth_hat <- predict(
  fit,
  newdata = dat,
  type = "smooth"
)
```

The lower-level `newX` and `newZ` arguments are retained for backward compatibility, but `newdata` is recommended for new code.

---

## Bootstrap inference

Bootstrap inference is available through:

```r
boot <- bootstrap_pgplsim(
  object = fit,
  B = 500,
  method = "stratified",
  seed = 2026
)
```

Bootstrap covariance estimates can also be obtained through `vcov()`.

---

## Package vignettes

The package includes three tutorials.

### Introduction

```r
vignette(
  "introduction",
  package = "pgplsim"
)
```

Introduces GPLSIM estimation, prediction, plotting, Poisson offsets, and inference.

### Simulation study

```r
vignette(
  "simulation-study",
  package = "pgplsim"
)
```

Illustrates a Monte Carlo simulation study and compares estimated nonlinear functions with the true data-generating mechanism.

### NHANES diabetes example

```r
vignette(
  "NHANES-diabetes",
  package = "pgplsim"
)
```

Demonstrates a real-data application using NHANES diabetes data.

---

## Documentation

```r
help(package = "pgplsim")

?pgplsim
?pgplsim_control
?predict.pgplsim
?plot.pgplsim
?summary.pgplsim
?vcov.pgplsim
?bootstrap_pgplsim
```

---

## Citation

If you use **pgplsim** in published work, please cite:

> Jun, Y.-B. (2026). **pgplsim: Penalized Spline Estimation for Generalized Partially Linear Single-Index Models**. R package.

A formal software paper is currently under preparation.

---

## Development roadmap

Potential future developments include:

- Cross-validation for smoothing selection
- Confidence bands for the smooth function
- Multiple-index GPLSIMs
- Parallel bootstrap
- Case weights
- Survey-weighted GPLSIMs
- Shiny interface
- CRAN release

---

## Contributing

Bug reports, feature requests, and pull requests are welcome.

Please open an issue at:

https://github.com/junpeea/pgplsim/issues

---

## License

MIT © 2026 Yoon-Bae Jun


## Gaussian PGPLSIM (identity link)

Gaussian support uses the same normalized single-index model and
penalized B-spline method as the generalized families. The smoothing
multiplier stored in `fit$lambda` is the coefficient multiplying the
penalty (the manuscript's squared lambda). The adjusted residual
variance is RSS divided by the effective residual degrees of freedom.

```r
set.seed(2026)
n <- 160
d <- data.frame(x1 = rnorm(n), z1 = runif(n, -1, 1),
                z2 = runif(n, -1, 1), z3 = runif(n, -1, 1))
u <- as.matrix(d[c("z1", "z2", "z3")]) %*% rep(1 / sqrt(3), 3)
d$y <- 0.7 * d$x1 + sin(2 * u) + rnorm(n, sd = 0.5)
fit <- pgplsim(y ~ x1, index = ~ z1 + z2 + z3,
               data = d, family = gaussian(), M = 5)
summary(fit)
predict(fit, type = "response") # identical to link scale
# Four analytic covariance choices from Sections 4–5 of Lu (2026):
V <- vcov(fit, type = "leverage")
boot <- bootstrap_pgplsim(fit, B = 20, method = "case", seed = 42)
vcov(fit, type = "bootstrap", bootstrap = boot)
```

For Gaussian models, offsets and outcome-stratified bootstrapping are
not supported. The Gaussian summary reports asymptotic normal-reference linear-coefficient
inference; effective residual degrees of freedom are used in estimating
residual dispersion, not to select a t reference distribution. Index-parameter
standard errors are supplied separately. Analytic and bootstrap
uncertainty estimates should be assessed by simulation.
