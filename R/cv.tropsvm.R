#' Cross-Validation for Tropical Support Vector Machines
#'
#' Conduct k-fold cross validation for tropsvm and return an object \code{"cv.tropsvm"}.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom parallel setDefaultCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom Rfast eachrow
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom lpSolve lp
#' @importFrom caret createFolds
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y a response vector with one label for each row/component of x.
#' @param parallel a logical value indicating if parallel computing should be used. (default: FALSE)
#' @param nfold a numeric value of the number of data folds for cross-validation. (default: 10)
#' @param nassignment a numeric value indicating the size of the parameter grid of assignments. (default: 10)
#' @param ncores a numeric value indicating the number of threads utilized for multi-cored CPUs. (default: 2)
#'
#' @return object with S3 class \code{cv.tropsvm} containing the fitted model, including:
#' \item{apex}{The negative apex of the fitted optimal tropical hyperplane.}
#' \item{assignment}{The best \code{assignment} tuned by cross-validation.}
#' \item{index}{The best classification method tuned by cross-validation.}
#' \item{levels}{The name of each category, consistent with categories in \code{y}.}
#' \item{accuracy}{The validation accuracy for each fold.}
#' \item{nfold}{The number of folds used in cross-validation.}
#'
#'
#' @seealso \code{summary}, \code{predict}, \code{coef} and the \code{tropsvm} function.
#'
#' @examples
#'
#' # data generation
#' library(Rfast)
#' set.seed(101)
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
#' summary(cv_tropsvm_fit)
#' coef(cv_tropsvm_fit)
#'
#' # test with new data
#' pred <- predict(cv_tropsvm_fit, newx)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(pred == newy) / length(newy)
#' @export
#' @export cv.tropsvm

cv.tropsvm <- function(x, y, parallel = FALSE, nfold = 10, nassignment = 10, ncores = 2) {
  if (nrow(x) != length(y)) {
    stop("numbers of data and label don't match")
  }
  if (nrow(x) <= nfold){
    stop("data set too small, please choose fewer folds to cross-validate")
  }
  if (length(unique(y)) != 2) {
    stop("only two classes allowded")
  }
  if (is.data.frame(x)) {
    warning("input data not 'matrix'; set to 'Matrix'")
    x <- data.matrix(x)
  }
  if (nassignment < 10) {
    warning("use more assignments when performance is bad.")
  }

  classes <- unique(y)
  reorder_ind <- c(which(y == classes[1]), which(y == classes[2]))
  np <- sum(y == classes[1])
  nq <- sum(y == classes[1])
  y <- y[reorder_ind]
  x <- x[reorder_ind, ]


  train_index <- createFolds(y, k = nfold, returnTrain = TRUE)
  all_assignment <- matrix(0, nrow = nfold * nassignment, ncol = 4)
  for (i in 1:length(train_index)) {
    P <- x[train_index[[i]] <= np, ]
    Q <- x[train_index[[i]] > np, ]
    all_assignment[((i - 1) * nassignment + 1):(i * nassignment), ] <- assignment_finder(P, Q)[1:nassignment, ]
  }
  all_assignment <- unique(all_assignment)
  all_assignment_list <- lapply(seq_len(nrow(all_assignment)), function(i) all_assignment[i, ])
  all_accuracy_list <- list()
  if (parallel) {cl <- makeCluster(ncores)}
  for (i in 1:length(train_index)) {
    # i = 1
    data <- x[train_index[[i]], ]
    label <- y[train_index[[i]]]
    n1 <- sum(label == classes[1])
    n2 <- sum(label == classes[2])
    n <- n1 + n2
    val_data <- x[-train_index[[i]], ]
    val_label <- y[-train_index[[i]]]
    val_n1 <- sum(val_label == classes[1])
    val_n2 <- sum(val_label == classes[2])
    val_n <- val_n1 + val_n2

    if (parallel) {
      all_accuracy <- parLapply(cl, all_assignment_list, function(assignment) {
        tropsvm_helper(x = data, y = label, assignment = assignment, ind = 1:70, newx = val_data, newy = val_label)
      })
    } else {
      all_accuracy <- lapply(all_assignment_list, function(assignment) {
        tropsvm_helper(data, label, assignment = assignment, ind = 1:70, newx = val_data, newy = val_label)
      })
    }
    accuracy_mat <- do.call("rbind", all_accuracy)
    all_accuracy_list[[i]] <- accuracy_mat
  }
  if (parallel) {stopCluster(cl)}
  all_accuracy_mat <- Reduce("+", all_accuracy_list)
  best_hyperparms <- matrix(which(all_accuracy_mat == max(all_accuracy_mat), arr.ind = T)[1, ], ncol = 2, byrow = TRUE)
  best_assignment <- all_assignment[best_hyperparms[1, 1], ]
  best_method_ind <- best_hyperparms[1, 2]
  best_accuracy <- sapply(all_accuracy_list, max)
  best_fold <- which.max(sapply(all_accuracy_list, function(x) {
    x[best_hyperparms]
  }))
  data <- x[train_index[[best_fold]], ]
  label <- y[train_index[[best_fold]]]
  n1 <- sum(label == classes[1])
  n2 <- sum(label == classes[2])
  n <- n1 + n2
  val_data <- x[-train_index[[best_fold]], ]
  val_label <- y[-train_index[[best_fold]]]
  val_n1 <- sum(val_label == classes[1])
  val_n2 <- sum(val_label == classes[2])
  val_n <- val_n1 + val_n2
  tropsvm.out <- tropsvm(data, label, assignment = best_assignment, ind = best_method_ind)
  cv.tropsvm.out <- list(
    "apex" = tropsvm.out$apex,
    "assignment" = best_assignment,
    "index" = best_method_ind,
    "levels" = as.factor(classes),
    "accuracy" = best_accuracy,
    "nfold" = nfold
  )
  class(cv.tropsvm.out) <- "cv.tropsvm"
  cv.tropsvm.out
}
