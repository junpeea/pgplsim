# CRAN submission comments

## First submission

This is the first CRAN submission of `pgplsim`, version 0.1.0.

## R CMD check results

The package was checked on the following platforms:

| Platform | Errors | Warnings | NOTES |
|---|---:|---:|---:|
| Windows 11, R 4.4.0 (local) | 0 | 0 | 1 |
| Ubuntu Linux, R 4.6.1 (R-hub) | 0 | 0 | 0 |
| Windows, R-devel (R-hub) | 0 | 0 | 0 |
| Ubuntu Linux, R 4.6.2 beta (R-hub) | 0 | 0 | 1 |

The local Windows NOTE was:

```
checking for future file timestamps ... NOTE
unable to verify current time
```

This appears to be related to the local environment's inability to verify the current time.

The Ubuntu next NOTE was:

```
checking tests ... NOTE
Running R code in 'testthat.R' had CPU time
3.6 times elapsed time
```

All tests completed successfully. The difference between CPU time and elapsed time may be related to multithreaded numerical computations.

## Additional validation

- All 188 tests passed in the local test suite.
- Package examples and documentation checks passed.
- All three vignettes built successfully.
- Independent R-hub checks completed successfully on all three selected platforms.
- No reverse dependencies are known for this first submission.
