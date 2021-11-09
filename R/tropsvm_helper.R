#' Helper for \code{cv.tropsvm}
#' @keywords internal
#' Compute classification accuracy for a pair of assignment and classification method in cross validation
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom Rfast eachrow
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom lpSolve lp
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y a response vector with one label for each row/component of x.
#' @param assignment a numeric vector of length 4 indicating the sectors of tropical hyperplane that the
#' data will be assigned to. The first and third elements in the \code{assignment} are the coordinates of
#' an observed point in data matrix \code{x} believed from the first category where the maximum and second maximum
#' of the vector addition between the fitted optimal tropical hyperplane and the point itself are achieved.
#' The meanings for the second and the fourth element in the \code{assignment} are the same
#' but for the points in the second category. Namely, the first and second values in the \code{assignment}
#' are the indices of sectors where the two point cloud will be assigned. Not needed when \code{auto.assignment = TRUE}. (default: NULL)
#' @param ind a numeric value or a numeric vector ranging from 1 to 70 indicating which classification method
#' to be used. There are 70 different classification methods. Details of a given method can be retrieved by \code{summary}.
#'  The different classification methods are proposed to resolve the issue when points fall on the intersections of sectors.
#'  Users can have personal choices if better knowledge is assumed. (default: 1)
#' @param newx the same as "x" but only needed in \code{cv.tropsvm}, which is used as validation data. (default: \code{NULL})
#' @param newy the same as "y" but only needed in \code{cv.tropsvm}, which is used as validation labels (default: \code{NULL})
#'
#' @return Classification accuracy
#'
#'
#'
#'
tropsvm_helper <- function(x, y, assignment = NULL, ind = 1, newx = NULL, newy = NULL) {
  classes <- unique(y)
  reorder_ind <- c(which(y == classes[1]), which(y == classes[2]))
  label <- y[reorder_ind]
  data <- x[reorder_ind, ]
  n1 <- sum(label == classes[1])
  n2 <- sum(label == classes[2])
  n <- n1 + n2

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
  reorder_ind <- c(which(newy == classes[1]), which(newy == classes[2]))
  val_label <- newy[reorder_ind]
  val_data <- newx[reorder_ind, ]
  val_n1 <- sum(val_label == classes[1])
  omega <- rep(0, ncol(data))
  omega[c(ip, jp, iq, jq)] <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution[2:5]
  omega[-c(ip, jp, iq, jq)] <- colMins(-data[, -c(ip, jp, iq, jq)] + c(data[1:n1, jp] + omega[jp], data[-c(1:n1), jq] + omega[jq]), T)
  shifted_val_data <- eachrow(val_data, omega, "+")
  diff <- eachrow(t(shifted_val_data), rowMaxs(shifted_val_data, T), oper = "-")
  raw_classification <- lapply(lapply(seq_len(ncol(diff)), function(i) diff[, i]), function(x) {
    which(abs(x) < 1e-10)
  })

  if (length(unique(assignment)) == 2) {
    accuracy <- (sum(raw_classification[1:val_n1] == ip) + sum(raw_classification[-c(1:val_n1)] == iq)) / length(raw_classification)
  } else {
    all_method_ind <- comboGeneral(8, 4)
    # Algorithm 1
    if (length(unique(assignment)) == 4) {
      P_base <- matrix(c(
        1, 0, 0, 0,
        0, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 1, 1
      ), ncol = 4, byrow = T)
      Q_base <- matrix(c(
        0, 0, 1, 0,
        0, 0, 0, 1,
        0, 0, 1, 1,
        0, 0, 0, 0
      ), ncol = 4, byrow = T)
      PQ_com <- matrix(c(
        1, 0, 1, 0,
        1, 0, 0, 1,
        0, 1, 1, 0,
        0, 1, 0, 1,
        1, 1, 1, 0,
        1, 1, 0, 1,
        1, 0, 1, 1,
        0, 1, 1, 1
      ), ncol = 4, byrow = T)
    }
    if (length(unique(assignment)) == 3) {
      P_base <- c()
      Q_base <- c()
      PQ_com <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        1, 1, 0,
        1, 0, 1,
        0, 1, 1,
        1, 1, 1,
        0, 0, 0
      ), ncol = 3, byrow = T)
      if (ip == jq) {
        PQ_com <- PQ_com[, c(1, 2, 3, 1)]
      }
      if (iq == jp) {
        PQ_com <- PQ_com[, c(1, 2, 2, 3)]
      }
      if (jp == jq) {
        PQ_com <- PQ_com[, c(1, 2, 3, 2)]
      }
    }
    colnames(PQ_com) <- c("ip", "jp", "iq", "jq")

    accuracy <- sapply(ind, function(l) {
      P <- rbind(P_base, PQ_com[all_method_ind[l, ], ])
      Q <- rbind(Q_base, PQ_com[-all_method_ind[l, ], ])
      sum(c(sapply(raw_classification[1:val_n1], function(x) {
        v <- c(ip, jp, iq, jq) %in% x
        return(sum(colSums(t(P) == v) == ncol(P)))
      }), sapply(raw_classification[-c(1:val_n1)], function(x) {
        v <- c(ip, jp, iq, jq) %in% x
        return(sum(colSums(t(Q) == v) == ncol(Q)))
      }))) / length(raw_classification)
    })
  }
  accuracy
}
