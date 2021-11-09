#' Efficient Finder of Best Assignment
#' @keywords internal
#' @importFrom Rfast colsums
#'
#' @param P a data matrix of the first category.
#' @param Q a data matrix of the second category.
#' @param t a numeric value controlling the number of selected best assignments.
#'
#' @return A matrix with each row a unique assignment, starting from the best assignment
#' and ending with the worst.
#'
assignment_finder <- function(P, Q, t = 25) {
  e <- ncol(P)
  n1 <- nrow(P)
  n2 <- ncol(Q)
  P_colsum <- colsums(P)
  Q_colsum <- colsums(Q)

  ip_stat <- 2 * P_colsum - Q_colsum
  iq_stat <- 2 * Q_colsum - P_colsum
  jp_stat <- -Q_colsum
  jq_stat <- -P_colsum

  ip_rank <- order(ip_stat, decreasing = TRUE)
  iq_rank <- order(iq_stat, decreasing = TRUE)
  jp_rank <- order(jp_stat, decreasing = TRUE)
  jq_rank <- order(jq_stat, decreasing = TRUE)
  t <- ifelse(e > t, t, e)
  all_rank <- rbind(ip_rank, iq_rank, jp_rank, jq_rank)
  temp <- as.matrix(expand.grid(
    ip = ip_rank[1:t], iq = iq_rank[1:t],
    jp = jp_rank[1:t], jq = jq_rank[1:t]
  ))
  temp <- temp[temp[, 1] - temp[, 2] != 0, ]
  temp <- temp[temp[, 1] - temp[, 3] != 0, ]
  temp <- temp[temp[, 2] - temp[, 4] != 0, ]
  temp[order(ip_stat[temp[, 1]] + iq_stat[temp[, 2]] + jp_stat[temp[, 3]] + jq_stat[temp[, 4]], decreasing = T), ]
}
