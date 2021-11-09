#' Dimension Reduction by Isometry
#' @keywords internal
#' Map points in row space of a tropical polytope to its column space.
#'
#' @param D a data matrix of size s x e whose rows are vertices of a tropical polytope.
#' @param P a vector of a point in row spaces of \code{D}.
#'
#' @return A vector of a point corresponding to \code{P} in column span of \code{D}.
#
polytope_iso <- function(D, P) {
  e <- length(P)
  s <- dim(D)[[1]]
  Q <- mat.or.vec(1, s)
  for (i in seq(s)) {
    maxvalue <- D[i, 1] - P[[1]]
    for (j in seq(e)) {
      maxvalue <- max(maxvalue, D[i, j] - P[[j]])
    }
    Q[[i]] <- maxvalue
  }
  return(Q)
}
