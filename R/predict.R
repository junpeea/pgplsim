############################################################
# Prediction method
############################################################


#' @export

predict.LuGPLSIM <- function(
    object,
    newX,
    newZ,
    type=c(
      "response",
      "link",
      "index"
    ),
    ...
){


  type <- match.arg(type)


  u <-
    as.numeric(
      newZ %*%
        object$alpha_hat
    )


  B <-
    predict_basis(
      u,
      object$basis$B
    )


  eta <-
    newX %*%
    object$beta_hat +
    B %*%
    object$gamma_hat



  if(type=="index")
    return(u)


  if(type=="link")
    return(eta)



  inverse_link(
    eta,
    object$family
  )

}
