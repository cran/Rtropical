#' Vectorize a Set of Phylognetic Trees
#'
#' Unifies tip labels of all phylogenetic trees in \code{multiPhylo}
#' object the same as the first tree and returns the
#' cophenetic distance of their corresponding chronogram.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom ape root
#' @importFrom ape chronos
#' @param x an object of class \code{multiPhylo}
#' @param tipOrder a numeric vector of order of leaf names to which all trees
#' in the \code{multiPhylo} object will unified. If not specified on purpose,
#' the tip order of the first tree will be used.
#' @param parallel a logical value indicating if parallel computing should be used. (default: FALSE)
#' @param ncores a numeric value indicating the number of threads
#' utilized for multi-cored CPUs. (default: 2)
#' @param \dots Not used. Other arguments to as.vector
#'
#' @return A data matrix with each row a vector representation of a chronogram. Each element of the vector is the distance between two leaves.
#'
#' @examples
#' data(apicomplexa)
#' data <- as.matrix(apicomplexa[1: 10]) # matrixize first ten trees
#'
#' @method as.matrix multiPhylo
#' @export
#' @export as.matrix.multiPhylo
as.matrix.multiPhylo<- function(x, tipOrder = x[[1]]$tip.label, parallel = FALSE, ncores = 2, ...) {
  trees_root <- root(x, outgroup = tipOrder[1], resolve.root = TRUE)
  if (parallel){
    cl <- makeCluster(ncores)
    chronotrees <- parLapply(cl, trees_root, chronos)
    distVec_all <- parLapply(cl, chronotrees, as.vector.phylo)
    stopCluster(cl)
  } else{
    chronotrees <- lapply(trees_root, chronos)
    distVec_all <- lapply(chronotrees, as.vector.phylo)
  }

  do.call("rbind", distVec_all)
}
