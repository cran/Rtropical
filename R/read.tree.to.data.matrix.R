#' Read Newick-formatted trees in two categories into a data matrix
#'
#' @importFrom ape read.tree
#' @param data.file1 A file containing trees in Newick form in a category.
#' @param data.file2 A file containing trees in Newick form in an assumed different category.
#'
#' @return \code{read.tree.to.data.matrix} has the same return as \code{read.nexus.to.data.matrix}.
#' @export
#'
read.tree.to.data.matrix <- function(data.file1, data.file2) {
  G1 <- read.tree(data.file1)
  G2 <- read.tree(data.file2)

  n <- length(G1[[1]]$tip.label)
  to <- G1[[1]]$tip.label
  N1 <- length(G1)
  N2 <- length(G2)

  distVec_all1 <- as.matrix(G1, to)
  distVec_all2 <- as.matrix(G2, to)
  rownames(distVec_all1) <- NULL
  rownames(distVec_all2) <- NULL

  class <- as.factor(c(rep(1, N1), rep(2, N2)))
  D <- data.frame(class, rbind(distVec_all1, distVec_all2))
  return(D)
}
