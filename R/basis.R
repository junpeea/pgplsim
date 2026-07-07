############################################################
# B-spline basis construction
############################################################


#' Construct B-spline basis and penalty matrix
#'
#' @param u Numeric vector of single-index values.
#' @param M Number of B-spline basis functions.
#' @param degree Degree of the B-spline basis. Default is 3.
#'
#' @return A list containing the B-spline basis matrix and penalty matrix.
#'
#' @importFrom splines bs
#' @export
compute_basis <- function(
    u,
    M = 80,
    degree = 3
){

  B <- splines::bs(
    u,
    df = M,
    degree = degree,
    intercept = TRUE
  )


  D <- difference_matrix(
    ncol(B),
    order = 2
  )


  P <- t(D) %*% D


  list(
    B = as.matrix(B),
    penalty = P
  )
}



# prediction basis
predict_basis <- function(
    newu,
    old_basis
){

  splines::bs(
    newu,
    df=ncol(old_basis),
    intercept=TRUE
  )
}
