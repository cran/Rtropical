#' Tropical Principal Component Analysis by Tropical Linear Space
#'
#' Approximate the principal component as a tropical linear space
#' for a given data matrix and returns the results as an object of class \code{troppca}.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#'
#' @param x a data matrix, of size n x e, with each row an observation vector.
#' e is the dimension of the tropical space
#' @param pcs a numeric value indicating the order of principal component. (default: 2)
#' @param iteration a list with arguments controlling the iteration of the algorithm.
#' \describe{
#' \item{exhaust}{a logical variable indicating if to iterate over all possible combinations of the linear space
#' based on the given data matrix \code{x}. If FALSE, please input a number of iteration for \code{niter}.
#' If TRUE, please enter 0 for \code{niter} and this function will iterate over all possible combinations of linear space.
#' This could be time consuming when \code{x} is large. (default: FALSE)}
#' \item{niter}{a numeric variable indicating the number of iterations. (default: 100)}
#' }
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
#' n <- 100
#' e <- 10
#' sig2 <- 1
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(sig2, e)))
#' troppca_fit <- troppca.linsp(x)
#' }
#'
#' @export
#' @export troppca.linsp
troppca.linsp <- function(x, pcs = 2, iteration = list(), ncores = 2) {
  con <- list(
    exhaust = FALSE,
    niter = 100
  )
  con[names(iteration)] <- iteration

  exhaust <- con$exhaust
  niter <- con$niter
  pcs <- pcs + 1
  if (exhaust) {
    warning("Iterating all possible PC choices enabled, this could take long for high order PCs.")
    all_choices <- comboGeneral(nrow(x), pcs)
    all_choices <- lapply(1:nrow(all_choices), function(i) all_choices[i, ])
  } else {

    all_choices <- lapply(1: niter, function(i){sample(1: nrow(x), pcs, replace = F)})
  }
  cl <- makeCluster(ncores)
  all_objs <- unlist(parLapply(cl, all_choices, function(ind) {
    # ind = 1
    V <- matrix(x[ind, ], nrow = length(ind))
    proj <- tropproj.linsp(x, V)
    temp <- x - proj
    sum(rowMaxs(temp, T) - rowMins(temp, T))
  }))
  stopCluster(cl)
  best_choice <- all_choices[[which.min(all_objs)]]
  pc <- matrix(x[best_choice, ], nrow = pcs)
  rownames(pc) <- paste("pc", 1:pcs, sep = "")
  proj_points <- tropproj.linsp(x, pc)
  rownames(proj_points) <- rownames(x)
  colnames(proj_points) <- colnames(x)
  troppca.out <- list(
    "pc" = pc,
    "obj" = min(all_objs),
    "projection" = proj_points,
    "type" = "linear space"
  )
  class(troppca.out) <- "troppca"
  troppca.out
}
