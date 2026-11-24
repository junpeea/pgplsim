# CRAN submission comments — pgplsim 0.1.0

## Current local check (November 23, 2026)

- Windows 11 x64, R 4.4.0 (2024-04-24 ucrt).
- `devtools::test()`: 308 passed, 0 failed, 0 warnings, 0 skipped (user-provided Step 1.6 output).
- `devtools::check()`: 0 errors, 0 warnings, 1 NOTE.
- NOTE: `checking for future file timestamps ... unable to verify current time`.
- All package tests, documentation, examples, and vignette rebuilding passed in that local check.

The timestamp NOTE arises during R's environment time verification; no source-code error was reported. We will recheck with CRAN incoming settings and on additional platforms prior to submission.


## Test environments

The package was checked using R-hub on the following configurations:

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

The test suite reported 403 passing expectations,
with one intentional CRAN-specific test skip.

Local Windows checks using R 4.4.0 also passed
with 0 ERRORs and 0 WARNINGs.

The local source-package check reported two NOTEs:
- A GitHub Actions workflow URL that has since
  been corrected.
- An inability to verify the current time.

Neither issue was reproduced in the seven R-hub
configurations.


## New functionality

The first CRAN release candidate now includes Gaussian identity-link PGPLSIM alongside binomial and Poisson models. Gaussian methods include Fisher-scoring estimation, adaptive smoothing, dispersion diagnostics, analytical covariance options, prediction, and case-bootstrap inference.

## Prior R-hub checks — historical only

The following earlier R-hub results predate the Gaussian integration and the second-difference penalty correction, so they **must not** be described as current cross-platform validation:

| Platform | Errors | Warnings | NOTEs |
|---|---:|---:|---:|
| Ubuntu Linux, R 4.6.1 | 0 | 0 | 0 |
| Windows, R-devel | 0 | 0 | 0 |
| Ubuntu Linux, R 4.6.2 beta | 0 | 0 | 1 |

The earlier Ubuntu NOTE reported CPU time 3.6 times elapsed time in testthat. These results must be refreshed using the final Gaussian-enabled source.

## Pending final verification

- Run current three-family compatibility tests and full suite after documentation updates.
- Repeat local check with strict CRAN incoming settings, where feasible.
- R-hub checks on current Linux, Windows, and macOS targets.
- Review build output, manual, license, README, and DESCRIPTION before submission.

## Simulation evidence and scope

A Gaussian pilot simulation (n = 400, 30 replications) converged in 30/30 runs; the full 1,000-replication study was intentionally skipped. This is preliminary evidence, not definitive coverage calibration.

## First submission

This is intended as the first CRAN submission of pgplsim 0.1.0. No reverse dependencies are known.
