#' Tropical Principal Component Analysis by Tropical Polytope
#'
#' Approximates the principal component as a tropical polytope for a given data matrix
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
#' troppca_fit <- troppca.poly(x)
#' plot(troppca_fit)
#' }
#'
#' @export
#' @export troppca.poly

troppca.poly <- function(x, pcs = 2, nsample = 1000, ncores = 2) {
  pcs <- pcs + 1
  n <- nrow(x)
  cl <- makeCluster(ncores)
  x_list <- lapply(seq_len(n), function(i) x[i, ])
  troppca_objs <- vector(mode = "numeric", nsample)
  samples <- matrix(NA, nrow = nsample, ncol = pcs)
  samples[1, ] <- sample(1:n, pcs)
  troppca_objs[1] <- troppca.obj(t(x[samples[1, ], ]), x_list, cl)

  t <- 1
  while (t < nsample) {
    # Find a new proposal by changing a randomly selected vertex of the current polytope
    current_choice <- samples[t, ]
    current_obj <- troppca_objs[t]

    change_ind <- sample(pcs, 1)
    out_change <- sample(c(1:n)[-current_choice], 1)
    new_choice <- c(current_choice[-change_ind], out_change)
    new_obj <- troppca.obj(t(x[new_choice, ]), x_list, cl)

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
  proj_points <- do.call("rbind", parLapply(cl, x_list, tropproj.poly, tconv = t(pc)))
  stopCluster(cl)
  rownames(pc) <- paste("pc", 1:pcs, sep = "")
  rownames(proj_points) <- rownames(x)
  colnames(proj_points) <- colnames(x)
  troppca.out <- list(
    "pc" = pc,
    "obj" = troppca_objs[min_index],
    "projection" = proj_points,
    "type" = "polytope"
  )
  class(troppca.out) <- "troppca"
  troppca.out
}
