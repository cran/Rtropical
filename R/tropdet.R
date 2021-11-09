#' Tropical Determinant of a Matrix
#'
#' Compute the tropical determinant for a given matrix. This is equivalent to
#' solving an assignment problem.
#'
#' @importFrom lpSolve lp.assign
#' @param x a square matrix
#' @return The determinant of the given matrix,
#' @examples
#' R <- matrix(sample(1:9, 9), nrow = 3)
#' tropdet(R)
#' @export
#' @export tropdet
tropdet <- function(x) {
  if (nrow(x) <= 5) {
    switch(ncol(x),
      x,
      trop_det2d(x),
      trop_det3d(x),
      trop_det4d(x),
      trop_det5d(x)
    )
  } else {
    lp.assign(x, "max")$objval
  }
}
