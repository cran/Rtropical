#' Projection on Tropical Linear Space
#'
#' Compute projection of data points on a given tropical linear space.
#'
#' @importFrom Rfast rowSort
#'
#' @param x a data matrix, of size n x e, with each row an observation.
#' @param V a data matrix, of dimension s x e, with each row a basis of tropical linear space.
#' e is the dimension of the tropical space and s is the dimension of the linear space.
#'
#' @return A matrix of projections of all data points.
#' @examples
#' library(Rfast)
#' n <- 100
#' e <- 10
#' sig2 <- 1
#' s <- 3
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(sig2, e)))
#' V <- matrix(runif(s * e, -10, 10), nrow = s, ncol = e)
#' x_proj <- tropproj.linsp(x, V)
#' head(x_proj)
#' @export
#' @export tropproj.linsp
tropproj.linsp <- function(x, V) {
  if (is.vector(x)) x <- t(as.matrix(x))
  n <- nrow(x)
  pcs <- nrow(V)
  e <- ncol(x)
  all_dets <- array(NA, dim = (e - pcs +1): e)
  all_combns <- comboGeneral(1:e, pcs)
  all_combns_list <- lapply(1:nrow(all_combns), function(i) {
    all_combns[i, ]
  })
  all_dets[all_combns] <- unlist(lapply(all_combns_list, function(i) {
    # i = 1
    tropdet(as.matrix(V[, i]))
  }))
  if (pcs == 1){
    data_proj <- matrix(rep(V, n), nrow = n, byrow = T)
  } else{
    data_proj <- matrix(0, nrow = n, ncol = e)
    for (i in 1:e) {
      # i = 1
      all_tau <- comboGeneral(c(1:e)[-i], pcs - 1)
      all_tau <- lapply(1:nrow(all_tau), function(i) {
        all_tau[i, ]
      })
      temp2 <- lapply(all_tau, function(tau) {
        # tau = all_tau[[1]]
        all_j <- c(1:e)[-tau]
        temp_block <- eachrow(matrix(x[, all_j], ncol = e - pcs + 1),
                              all_dets[rowSort(cbind(all_j, matrix(tau, ncol = length(tau), nrow = (e - pcs + 1), byrow = T)))], "-")
        return(rowMins(temp_block, T) + tropdet(V[, c(tau, i)]))
      })
      temp2 <- do.call(cbind, temp2)
      data_proj[, i] <- rowMaxs(temp2, T)
    }
  }
  rownames(data_proj) <- rownames(x)
  colnames(data_proj) <- colnames(x)
  data_proj
}
# tropproj.linsp <- function(x, V) {
#   if (is.vector(x)) x <- t(as.matrix(x))
#   n <- nrow(x)
#   pcs <- nrow(V)
#   e <- ncol(x)
#   all_dets <- array(NA, dim = c(e - 2, e - 1, e))
#   all_combns <- comboGeneral(1:e, pcs)
#   all_combns_list <- lapply(1:nrow(all_combns), function(i) {
#     all_combns[i, ]
#   })
#   all_dets[all_combns] <- unlist(lapply(all_combns_list, function(i) {
#     tropdet(V[, i])
#   }))
#
#   data_proj <- matrix(0, nrow = n, ncol = e)
#   for (i in 1:e) {
#     # i = 1
#     all_tau <- comboGeneral(c(1:e)[-i], pcs - 1)
#     all_tau <- lapply(1:nrow(all_tau), function(i) {
#       all_tau[i, ]
#     })
#     temp2 <- lapply(all_tau, function(tau) {
#       all_j <- c(1:e)[-tau]
#       temp_block <- eachrow(matrix(x[, all_j], ncol = e - pcs + 1), all_dets[rowSort(cbind(all_j, matrix(tau, ncol = 2, nrow = (e - pcs + 1), byrow = T)))], "-")
#       return(rowMins(temp_block, T) + tropdet(V[, c(tau, i)]))
#     })
#     temp2 <- do.call(cbind, temp2)
#     data_proj[, i] <- rowMaxs(temp2, T)
#   }
#   rownames(data_proj) <- rownames(x)
#   colnames(data_proj) <- colnames(x)
#   data_proj
# }
