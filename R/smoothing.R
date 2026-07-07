############################################################
# Fellner-Schall smoothing update
############################################################


update_lambda_FS <- function(
    gamma,
    penalty,
    lambda
){


  numerator <-
    sum(diag(
      lambda * penalty
    ))


  denominator <-
    as.numeric(
      t(gamma) %*%
        penalty %*%
        gamma
    )


  lambda_new <-
    numerator /
    (denominator + 1e-8)


  lambda_new
}
