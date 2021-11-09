#' Vectorize a Phylogenetic Tree
#'
#' Computes the cophenetic distance and outputs them in a vector
#' of a phylogenetic tree in \code{phylo} object
#'
#' @importFrom ape cophenetic.phylo
#' @param x A object of class phylo
#' @param mode The same as \code{base::as.vector}. But only numeric
#' output in vector form is accepted for other functions in \code{Rtropical}
#'
#' @return A vector with its elements the distance between two leaves of the tree.
#'
#' @examples
#' library(ape)
#' tree <- rcoal(5)
#' tree_vec <- as.vector(tree)
#'
#' @method as.vector phylo
#' @export
#' @export as.vector.phylo
#'
as.vector.phylo <- function(x, mode = "any") {
  cophe <- cophenetic.phylo(x)
  labs <- outer(rownames(cophe), colnames(cophe), paste, sep="--")
  phyvec <- cophe[lower.tri(cophe)]
  names(phyvec) <- labs[lower.tri(cophe)]
  as.vector(phyvec, mode)
}
