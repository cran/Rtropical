#' Read NEXUS-formatted trees from two categories into a data matrix
#'
#' @importFrom ape read.nexus
#' @param data.file1 A data set with trees from one category.
#' @param data.file2 A data set with trees from the other category.
#'
#' @return A data matrix with the first x rows corresponding the x trees in the first file and the last y rows are the trees from the second file.
#' @export
#'
read.nexus.to.data.matrix <- function(data.file1, data.file2) {
  G1 <- read.nexus(data.file1)
  G2 <- read.nexus(data.file2)

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
