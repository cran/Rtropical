#' Extract Optimal Tropical Hyperplane from a tropsvm object
#'
#' Obtain the optimal tropical hyperplane in the form of vectors from a tropsvm object.
#'
#' @param object a fitted \code{"tropsvm"} object.
#' @param \dots Not used. Other arguments.
#'
#' @return An output of the apex of the fitted optimal tropical hyperplane.
#' @method coef tropsvm
#' @export
#' @export coef.tropsvm
coef.tropsvm <- function(object, ...) {
  cat("The apex of seprarating hyperplane is: \n", -object$apex, "\n")
}
