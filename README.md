# PSPLINE

<!-- badges: start -->
[![R-CMD-check](https://github.com/junpeea/PSPLINE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/junpeea/PSPLINE/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

## Overview

**PSPLINE** is an R package for fitting **Generalized Partially Linear Single-Index Models (GPLSIMs)** using penalized B-splines.

The package provides a unified framework for modeling nonlinear covariate effects through a single-index structure while retaining interpretable linear effects. Automatic smoothing parameter estimation is performed using the Fellner–Schall update, and the package includes model-based, sandwich, and bootstrap variance estimation.

The current release supports

- Binary outcomes via logistic regression
- Count outcomes via Poisson regression
- Penalized B-spline estimation
- Automatic smoothing parameter selection
- Prediction for new observations
- Bootstrap inference
- Sandwich covariance estimation
- Publication-quality plotting

---

## Model

PSPLINE fits models of the form

\[
g\{E(Y|X,Z)\}
=
X^\top\beta
+
\phi(Z^\top\alpha),
\]

where

- \(X\) contains linear covariates,
- \(Z\) contains variables entering the nonlinear single index,
- \(\alpha\) is the index coefficient,
- \(\phi(\cdot)\) is an unknown smooth function estimated by penalized splines.

---

## Installation

### Development version

```r
install.packages("remotes")

remotes::install_github("junpeea/PSPLINE")
```

or

```r
devtools::install_github("junpeea/PSPLINE")
```

---

## Quick example

```r
library(PSPLINE)

set.seed(123)

n <- 500

X <- cbind(
  x1 = rnorm(n),
  x2 = rbinom(n,1,0.5)
)

Z <- cbind(
  z1 = rnorm(n),
  z2 = rnorm(n),
  z3 = rnorm(n)
)

alpha <- c(1,1,1)
alpha <- alpha/sqrt(sum(alpha^2))

eta <- X %*% c(1,-0.5) +
       sin(as.vector(Z %*% alpha))

prob <- plogis(eta)

y <- rbinom(n,1,prob)

fit <- LuGPLSIM(
  y = y,
  X = X,
  Z = Z,
  family = binomial()
)

summary(fit)

plot(fit)

pred <- predict(
  fit,
  newX = X,
  newZ = Z,
  type = "response"
)
```

---

## Main function

```r
LuGPLSIM(
    y,
    X,
    Z,
    family = binomial(),
    offset = NULL,
    M = NULL,
    lambda = 1,
    maxit = 100,
    tol = 1e-6
)
```

---

## Supported families

Currently supported

- `binomial()`
- `poisson()`

Future releases will include

- Gaussian
- Negative Binomial
- Zero-inflated Poisson
- Cox-type survival GPLSIMs

---

## Main features

- Penalized B-spline estimation
- Generalized partially linear single-index models
- Logistic GPLSIM
- Poisson GPLSIM
- Automatic smoothing parameter estimation
- Fisher scoring optimization
- Model-based covariance estimation
- Sandwich covariance estimation
- Bootstrap standard errors
- Prediction for new observations
- Publication-quality graphics

---

## Package vignettes

The package includes three tutorials.

### Introduction

```r
vignette("01-introduction", package="PSPLINE")
```

Introduces GPLSIM estimation, prediction, plotting, and inference.

---

### Simulation study

```r
vignette("02-simulation-study", package="PSPLINE")
```

Illustrates simulation studies and compares estimated nonlinear functions with the true data-generating mechanism.

---

### NHANES diabetes example

```r
vignette("03-NHANES-diabetes", package="PSPLINE")
```

Demonstrates a real-data application using NHANES diabetes data.

---

## Documentation

```r
help(package="PSPLINE")

?LuGPLSIM

?predict.LuGPLSIM

?plot.LuGPLSIM

?bootstrap_LuGPLSIM
```

---

## Citation

If you use **PSPLINE** in published work, please cite

> Jun, Y.-B. (2026). **PSPLINE: Penalized Spline Estimation for Generalized Partially Linear Single-Index Models**. R package.

A formal software paper is currently under preparation.

---

## Development roadmap

Upcoming features include

- Cross-validation for smoothing selection
- Confidence bands for the smooth function
- Multiple-index GPLSIM
- Parallel bootstrap
- Case weights
- Survey-weighted GPLSIM
- Negative Binomial GPLSIM
- Zero-inflated GPLSIM
- Shiny interface
- CRAN release

---

## Contributing

Bug reports, feature requests, and pull requests are welcome.

Please open an issue at

https://github.com/junpeea/PSPLINE/issues

---

## License

MIT © 2026 Yoon-Bae Jun
