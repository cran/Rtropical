#' Tropical Support Vector Machines
#'
#' Fit a discriminative two-class classifier via linear programming defined by the tropical
#' hyperplane which maximizes the minimum tropical distance from data points
#' to itself in order to separate the data points into sectors (half-spaces)
#' in the tropical projective torus.
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom Rfast eachrow
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom lpSolve lp
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y a response vector with one label for each row/component of x.
#' @param auto.assignment a logical value indicating if to provide an \code{assignment} manually.
#' If \code{FALSE}, an input is required, otherwise the function automatically
#' finds a good assignment.(default: FALSE)
#' @param assignment a numeric vector of length 4 indicating the sectors of tropical hyperplane that the
#' data will be assigned to. The first and third elements in the \code{assignment} are the coordinates of
#' an observed point in data matrix \code{x} believed from the first category where the maximum and second maximum
#' of the vector addition between the fitted optimal tropical hyperplane and the point itself are achieved.
#' The meanings for the second and the fourth element in the \code{assignment} are the same
#' but for the points in the second category. Namely, the first and second values in the \code{assignment}
#' are the indices of sectors where the two point clouds are assigned. Not needed when \code{auto.assignment = TRUE}. (default: NULL)
#' @param ind a numeric value or a numeric vector ranging from 1 to 70 indicating which classification method
#' to be used. There are 70 different classification methods. Details of a given method can be retrieved by \code{summary}.
#'  The different classification methods are proposed to resolve the issue when points fall on the intersections of sectors.
#'  Users can have personal choices if better knowledge is assumed. (default: 1)
#'
#' @return An object with S3 class \code{tropsvm} containing the fitted model, including:
#' \item{apex}{The negative apex of the fitted optimal tropical hyperplane.}
#' \item{assignment}{The user-input or auto-found \code{assignment}.}
#' \item{index}{The user-input classification method.}
#' \item{levels}{The name of each category, consistent with categories in \code{y}.}
#'
#' @seealso \code{predict}, \code{coef} and the \code{cv.tropsvm} function.
#'
#'
#' @examples
#'
#' # data generation
#' library(Rfast)
#' e <- 100
#' n <- 10
#' N <- 100
#' s <- 10
#' x <- rbind(
#'   rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
#'   rmvnorm(n, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
#' )
#' y <- as.factor(c(rep(1, n), rep(2, n)))
#' newx <- rbind(
#'   rmvnorm(N, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
#'   rmvnorm(N, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
#' )
#' newy <- as.factor(rep(c(1, 2), each = N))
#'
#' # train the tropical svm
#' tropsvm_fit <- tropsvm(x, y, auto.assignment = TRUE, ind = 1)
#'
#' coef(tropsvm_fit)
#'
#' # test with new data
#' pred <- predict(tropsvm_fit, newx)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(pred == newy) / length(newy)
#' @export
#' @export tropsvm
tropsvm <- function(x, y, auto.assignment = FALSE, assignment = NULL, ind = 1) {
  if (nrow(x) != length(y)) {
    stop("numbers of data and label don't match")
  }
  if (length(unique(y)) != 2) {
    stop("only two classes allowded")
  }
  if (is.data.frame(x)) {
    warning("input data not 'matrix'; set to 'Matrix'")
    x <- data.matrix(x)
  }
  classes <- unique(y)
  reorder_ind <- c(which(y == classes[1]), which(y == classes[2]))
  label <- y[reorder_ind]
  data <- x[reorder_ind, ]
  n1 <- sum(label == classes[1])
  n2 <- sum(label == classes[2])
  n <- n1 + n2

  if (auto.assignment) {
    assignment <- assignment_finder(data[1:n1, ], data[-c(1:n1), ])[1, ]
  }
  names(assignment) <- c("ip", "iq", "jp", "jq")
  ip <- assignment[1]
  jp <- assignment[3]
  iq <- assignment[2]
  jq <- assignment[4]
  f.obj <- c(1, rep(0, 4), c(
    rep(-1, n1), rep(-1, n1), rep(-1, n1), rep(-1, n1), rep(-1, n2),
    rep(-1, n2), rep(-1, n2), rep(-1, n2)
  ))
  f.conp <- rbind(
    cbind(rep(1, n1), rep(-1, n1), rep(1, n1), rep(0, n1), rep(0, n1)),
    cbind(rep(0, n1), rep(-1, n1), rep(1, n1), rep(0, n1), rep(0, n1)),
    cbind(rep(0, n1), rep(0, n1), rep(-1, n1), rep(1, n1), rep(0, n1)),
    cbind(rep(0, n1), rep(0, n1), rep(-1, n1), rep(0, n1), rep(1, n1))
  )
  f.conq <- rbind(
    cbind(rep(1, n2), rep(0, n2), rep(0, n2), rep(-1, n2), rep(1, n2)),
    cbind(rep(0, n2), rep(0, n2), rep(0, n2), rep(-1, n2), rep(1, n2)),
    cbind(rep(0, n2), rep(1, n2), rep(0, n2), rep(0, n2), rep(-1, n2)),
    cbind(rep(0, n2), rep(0, n2), rep(1, n2), rep(0, n2), rep(-1, n2))
  )
  f.con <- cbind(rbind(f.conp, f.conq), diag(-1, nrow = 4 * n, ncol = 4 * n))
  f.dir <- rep("<=", n)
  f.rhs <- c(
    rep(data[1:n1, ip] - data[1:n1, jp], 2),
    data[1:n1, jp] - data[1:n1, iq],
    data[1:n1, jp] - data[1:n1, jq],
    rep(data[-c(1:n1), iq] - data[-c(1:n1), jq], 2),
    data[-c(1:n1), jq] - data[-c(1:n1), ip],
    data[-c(1:n1), jq] - data[-c(1:n1), jp]
  )
  omega <- rep(0, ncol(data))
  sol <- lp("max", f.obj, f.con, f.dir, f.rhs)
  omega[c(ip, jp, iq, jq)] <- sol$solution[2:5]
  omega[-c(ip, jp, iq, jq)] <- colMins(-data[, -c(ip, jp, iq, jq)] + c(data[1:n1, jp] + omega[jp], data[-c(1:n1), jq] + omega[jq]), T)
  tropsvm.out <- list(
    "apex" = omega,
    "assignment" = assignment,
    "index" = ind,
    "levels" = as.factor(classes)
  )
  class(tropsvm.out) <- "tropsvm"
  tropsvm.out
}
