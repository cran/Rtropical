#' @name dets
#' @rdname dets
#'
#' @title Tropical Determinants for Small Matrices
#'
#' @param X a square matrix
#' @return Tropical determinant of a given matrix
#' @keywords internal
trop_det2d <- function(X) {
  max(X[1, 1] + X[2, 2], X[2, 1] + X[1, 2])
}
#' @rdname dets
trop_det3d <- function(X) {
  a1 <- max(X[2, 2] + X[3, 3], X[2, 3] + X[3, 2])
  a2 <- max(X[2, 1] + X[3, 3], X[2, 3] + X[3, 1])
  a3 <- max(X[2, 1] + X[3, 2], X[2, 2] + X[3, 1])
  return(max(X[1, ] + c(a1, a2, a3)))
}
#' @rdname dets
trop_det4d <- function(X) {
  XX <- X[-1, ]
  a1 <- trop_det3d(XX[, -1])
  a2 <- trop_det3d(XX[, -2])
  a3 <- trop_det3d(XX[, -3])
  a4 <- trop_det3d(XX[, -4])
  return(max(X[1, ] + c(a1, a2, a3, a4)))
}
#' @rdname dets
trop_det5d <- function(X) {
  XX <- X[-1, ]
  a1 <- trop_det4d(XX[, -1])
  a2 <- trop_det4d(XX[, -2])
  a3 <- trop_det4d(XX[, -3])
  a4 <- trop_det4d(XX[, -4])
  a5 <- trop_det4d(XX[, -5])
  return(max(X[1, ] + c(a1, a2, a3, a4, a5)))
}
