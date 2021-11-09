#' Tropical Fermat-Weber Point
#'
#' Compute the tropical Fermat-Weber (FW) point for a given data matrix.
#' The FW point minimizes the summed tropical distance to the trees described
#' in the data matrix.
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom lpSolveAPI make.lp
#' @importFrom lpSolveAPI set.constr.type
#' @importFrom lpSolveAPI set.rhs
#' @importFrom lpSolveAPI solve.lpExtPtr
#' @importFrom lpSolveAPI get.variables
#' @importFrom lpSolveAPI set.column
#' @importFrom lpSolveAPI get.objective
#' @importFrom lpSolveAPI set.objfn
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#'
#' @return A list containing:
#' \item{fw}{The fermat-weber point.}
#' \item{distsum}{The sum of distance from each observation to the fermat-weber point.}
#'
#'
#'
#'
#' @examples
#' x <- matrix(rnorm(100), ncol = 10)
#' tropFW(x)
#' @export
#' @export tropFW
tropFW <- function(x) {
  nn <- nrow(x)
  e <- ncol(x)
  jk <- comboGeneral(1:e, 2)
  combn_size <- nrow(jk)
  obj <- c(rep(1, nn), rep(0, e))
  conY <- matrix(0, nrow = nn * combn_size, ncol = e)
  all_v <- x[, jk[, 2]] - x[, jk[, 1]]
  conY[cbind(1:(nn * combn_size), rep(jk[, 1], each = nrow(x)))] <- -1
  conY[cbind(1:(nn * combn_size), rep(jk[, 2], each = nrow(x)))] <- 1
  conD <- matrix(0, nrow = 2 * nn * combn_size, ncol = nn)
  conD[cbind(1:(2 * nn * combn_size), 1:nn)] <- -1
  con <- cbind(conD, rbind(conY, -conY))
  rhs <- c(matrix(all_v), -matrix(all_v))

  lprec <- make.lp(nrow(con), nn + e)
  for (i in 1:ncol(con)) {
    set.column(lprec, i, -con[, i])
  }
  set.constr.type(lprec, rep(">=", (2 * nn * combn_size)))
  set.rhs(lprec, -rhs)
  set.objfn(lprec, c(rep(1, nn), rep(0, e)))
  solve.lpExtPtr(lprec)
  sols <- get.variables(lprec)
  list("fw" = sols[-c(1:nn)], "distsum" = sum(sols[1:nn]))
}
