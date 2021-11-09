#' Projection on Tropical Polytope
#'
#' Project a point onto a given tropical polytope.
#'
#' @param x a data vector, of length e.
#' @param tconv a data matrix, of size e x s, with each column a vertex of the tropical polytope.
#' e is the dimension of the tropical space and s is the number of vertices of the polytope
#'
#' @return A projected vector on the given tropical polytope.
#'
#'
#' @examples
#'
#' # Generate a tropical polytope consisting of three trees each with 5 leaves
#' library(ape)
#' pltp <- sapply(1:3, function(i) {
#'   as.vector(rcoal(5))
#' })
#' # Generate an observation and vectorize it
#' tree <- rcoal(5)
#' tree_vec <- as.vector(tree)
#' tropproj.poly(tree_vec, pltp)
#' @export
#' @export tropproj.poly
#'
tropproj.poly <- function(x, tconv) {
  if (is.null(dim(tconv))) {
    lambda <- min(x - tconv)
    pi_D <- c(t(lambda + t(tconv)))
  } else {
    lambda <- colMins(x - tconv, value = TRUE)
    pi_D <- rowMaxs(eachrow(tconv, lambda, "+"), value = TRUE)
  }
  return(pi_D)
}
