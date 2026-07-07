############################################################
# Plot estimated link function
############################################################


#' @export

plot.LuGPLSIM <- function(
    x,
    ...
){


  ord <- order(
    x$u_hat
  )


  plot(
    x$u_hat[ord],
    x$basis$B[ord,] %*%
      x$gamma_hat,
    type="l",
    xlab=
      expression(
        Z^T*alpha
      ),
    ylab=
      expression(
        hat(phi)
      ),
    ...
  )

}
