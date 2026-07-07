############################################################
# Variance estimation
############################################################


vcov.LuGPLSIM <- function(
    object,
    type =
      c(
        "model",
        "sandwich"
      ),
    ...
){

  type <- match.arg(type)


  if(type=="model"){

    return(
      object$vcov_model
    )

  }


  if(type=="sandwich"){

    return(
      object$vcov_robust
    )

  }

}
