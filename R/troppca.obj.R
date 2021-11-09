#' @name troppca.obj
#'
#' @title Compute Tropical PCA Objective
#' @keywords internal
#' @importFrom parallel parLapply
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast rowMins
#'
#' @param pc a matrix of principal components
#' @param x_list a list of vectors
#' @param cl cluster for parallel computing
#'
#' @return The numeric value of the objective function of tropical principle component analysis.
#' This is the sum of all tropical distance from each point to its projection on the tropical polytope.
#' @rdname troppca.obj
troppca.obj <- function(pc, x_list, cl) {
  proj <- parLapply(cl, x_list, tropproj.poly, tconv = pc)
  temp <- do.call("rbind", x_list) - do.call("rbind", proj)
  sum(rowMaxs(temp, value = T) - rowMins(temp, value = T))
}
#' @rdname troppca.obj
troppca.obj2 <- function(pc, x_list, cl) {
  pc <- linsp_to_poly(pc)
  proj <- parLapply(cl, x_list, tropproj.poly, tconv = t(pc))
  temp <- do.call("rbind", x_list) - do.call("rbind", proj)
  sum(rowMaxs(temp, value = T) - rowMins(temp, value = T))
}
