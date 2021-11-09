#' Convert Tropical Linear Space to Convex Hull
#' @keywords internal
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom Rfast rowSort
#' @param V a data matrix, of dimension s x e, with each row a basis of tropical linear space.
#' e is the dimension of the tropical space and s is the dimension of the linear space.
#' @return the return is a e choose s-1 x e matrix.
#' @export linsp_to_poly
linsp_to_poly <- function(V) {
  pcs <- nrow(V)
  e <- ncol(V)
  all_dets <- array(-1e10, dim = rep(e, pcs))
  all_combns <- comboGeneral(1:e, pcs)
  all_combns_list <- lapply(1:nrow(all_combns), function(i) {
    all_combns[i, ]
  })
  all_dets[all_combns] <- unlist(lapply(all_combns_list, function(i) {
    tropdet(V[, i])
  }))
  all_combns2 <- comboGeneral(1:e, pcs - 1)
  all_ind <- cbind(rep(1:e, each = nrow(all_combns2)), all_combns2[rep(seq(nrow(all_combns2)), e), ])
  matrix(all_dets[rowSort(all_ind)], ncol = e)
}
