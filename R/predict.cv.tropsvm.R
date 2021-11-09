#' Predict Method for Tropical Support Vector Machines based on Cross-Validation
#'
#' Predicts values based upon a model trained by \code{cv.tropsvm}.
#' @importFrom RcppAlgos comboGeneral
#'
#' @param object a fitted \code{"cv.tropsvm"} object.
#' @param newx a data matrix, of dimension nobs x nvars used as testing data.
#' @param \dots Not used. Other arguments to predict.

#' @return A vector of predicted values of a vector of labels.
#'
#'
#' @seealso \code{summary}, \code{coef} and the \code{cv.tropsvm} function.
#'
#' @examples
#'
#' # data generation
#' library(Rfast)
#' e <- 20
#' n <- 10
#' N <- 10
#' s <- 5
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
#' cv_tropsvm_fit <- cv.tropsvm(x, y, parallel = FALSE)
#'
#' # test with new data
#' pred <- predict(cv_tropsvm_fit, newx)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(pred == newy) / length(newy)
#' @method predict cv.tropsvm
#'
#' @export
#' @export predict.cv.tropsvm
predict.cv.tropsvm <- function(object, newx, ...) {
  predict.tropsvm(object, newx)
}
