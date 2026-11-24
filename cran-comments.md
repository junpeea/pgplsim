
# CRAN submission comments — pgplsim 0.1.0

## First submission

This is the first CRAN submission of pgplsim.

The package implements penalized generalized partially linear
single-index models for Gaussian, binomial, and Poisson outcomes.

## Test environments

The package was checked using R-hub on:

- Linux, R-devel (GCC)
- Windows, R-devel
- macOS Intel, R-devel
- macOS ARM64, R-devel
- Ubuntu Linux, R-release
- Ubuntu Linux, R-devel (Clang)
- Ubuntu Linux, R-devel (Clang ASAN)

All seven checks completed with:

- 0 ERRORs
- 0 WARNINGs
- 0 NOTEs

The test suite reported 403 passing expectations, with one
intentional CRAN-specific test skip.

## Local Windows check

Windows 11 x64, R 4.4.0:

- 0 ERRORs
- 0 WARNINGs
- 1 NOTE

The NOTE was:

checking for future file timestamps ... unable to verify current time

This NOTE appears to be related to local environment time
verification. It did not occur in the seven R-hub checks.

## Package functionality

The package supports:

- Gaussian models with an identity link
- Binomial models with a logit link
- Poisson models with a log link

It provides model fitting, prediction, analytical covariance
estimation, and bootstrap inference.

## Additional information

- Package version: 0.1.0
- License: MIT
- No reverse dependencies are known.
