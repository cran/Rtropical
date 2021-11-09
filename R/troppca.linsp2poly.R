#' Tropical Principal Component Analysis by Polytope Converted from Linear Space
#'
#' Approximate the principal component as a tropical polytope converted from tropical linear space for a given data matrix
#' via MCMC and return the results as an object of class \code{troppca}.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom stats runif
#'
#' @param x a data matrix, of size n x e, with each row an observation vector.
#' e is the dimension of the tropical space#'
#' @param pcs a numeric value indicating the order of principal component. (default: 2)
#' @param nsample a numeric value indicating the number of samples of MCMC. (default: 1000)
#' @param ncores a numeric value indicating the number of threads utilized for multi-cored CPUs. (default: 2)
#'
#' @return A list of S3 class \code{"troppca"}, including:
#' \item{pc}{The principal component as a tropical linear space}
#' \item{obj}{The tropical PCA objective, the sum of tropical distance from each point to the projection.}
#' \item{projection}{The projections of all data points.}
#' \item{type}{The geometry of principal component.}
#'
#'
#' @examples
#' \donttest{
#' library(Rfast)
#' n <- 50
#' e <- 50
#' s <- 5
#' x <- rbind(
#'   rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
#'   rmvnorm(n, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
#' )
#' troppca_fit <- troppca.linsp2poly(x)
#' }
#'
#' @export
#' @export troppca.linsp2poly
# troppcalinsp2poly = function(x, pcs = 2, iteration = list(), ncores = 2){
#   con <- list(
#     exhaust = FALSE,
#     niter = 100
#   )
#   con[names(iteration)] <- iteration
#
#   x_list <- lapply(seq_len(nrow(x)), function(i) x[i, ])
#   exhaust <- con$exhaust
#   niter <- con$niter
#   pcs <- pcs + 1
#   all_choices <- comboGeneral(nrow(x), pcs)
#   if (exhaust){
#     all_choices <- lapply(1: nrow(all_choices), function(i) all_choices[i, ])
#   } else{
#     all_choices <- lapply(sample(1: nrow(all_choices), niter, replace = F), function(i) all_choices[i, ])
#   }
#   cl <- makeCluster(ncores)
#   all_objs <- unlist(parLapply(cl, all_choices, function(ind){
#     V <- x[ind, ]
#     proj = do.call("rbind", lapply(x_list, tropproj.poly, t(linsp_to_poly(V))))
#     temp <- x - proj
#     sum(rowMaxs(temp, T) - rowMins(temp, T))
#   }))
#   stopCluster(cl)
#   best_choice <- all_choices[[which.min(all_objs)]]
#   pc <- x[best_choice, ]
#   rownames(pc) <- paste("pc", 1: pcs, sep = "")
#   proj_points <- do.call("rbind", lapply(x_list, tropproj.poly, t(linsp_to_poly(pc))))
#   troppca.out <- list("pc" = pc,
#                      "obj" = min(all_objs),
#                      "projection" = proj_points,
#                      "type" = "linear space")
#   class(troppca.out) <- "troppca"
#   troppca.out
# }
troppca.linsp2poly <- function(x, pcs = 2, nsample = 1000, ncores = 2) {
  pcs <- pcs + 1
  n <- nrow(x)
  cl <- makeCluster(ncores)
  x_list <- lapply(seq_len(n), function(i) x[i, ])
  troppca_objs <- vector(mode = "numeric", nsample)
  samples <- matrix(NA, nrow = nsample, ncol = pcs)
  samples[1, ] <- sample(1:n, pcs)
  troppca_objs[1] <- troppca.obj2(x[samples[1, ], ], x_list, cl)

  t <- 1
  while (t < nsample) {
    # Find a new proposal by changing a randomly selected vertex of the current polytope
    current_choice <- samples[t, ]
    current_obj <- troppca_objs[t]

    change_ind <- sample(pcs, 1)
    out_change <- sample(c(1:n)[-current_choice], 1)
    new_choice <- c(current_choice[-change_ind], out_change)
    new_obj <- troppca.obj2(x[new_choice, ], x_list, cl)

    # Compute the probability we accept the new PCA base
    p <- min(1, current_obj / new_obj)

    if (sample(c(0, 1), 1, prob = c(1 - p, p)) == 1) {
      samples[(t + 1), ] <- new_choice
      troppca_objs[(t + 1)] <- new_obj
      t <- t + 1
    }
  }
  min_index <- which(troppca_objs == min(troppca_objs))[1]
  best_obj <- troppca_objs[min_index]
  pc <- x[samples[min_index, ], ]
  pc2 <- linsp_to_poly(pc)
  proj_points <- do.call("rbind", parLapply(cl, x_list, tropproj.poly, tconv = t(pc2)))
  stopCluster(cl)
  rownames(pc) <- paste("pc", 1:pcs, sep = "")
  troppca.out <- list(
    "pc" = pc,
    "obj" = troppca_objs[min_index],
    "projection" = proj_points,
    "type" = "linear space"
  )
  class(troppca.out) <- "troppca"
  troppca.out
}
