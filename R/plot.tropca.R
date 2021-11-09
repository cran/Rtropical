#' Plot the Tropical Principal Components with Data Projections
#'
#' Visualize the second order tropical principle components in \code{troppca}
#' as a tropical triangle with projections on a two-dimensional plot via tropical isometry.
#'
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics par
#' @importFrom graphics legend
#'
#' @param x a fitted \code{troppca} object.
#' @param plab a vector of labels of all points in the given data matrix.
#' Not needed for unlabeled data. (default: NULL)
#' @param fw a logical variable to determine if to add
#' Fermat-Weber point of the data projection. (default: FALSE)
#' @param \dots Not used. Other arguments to plot
#'
#' @return \code{plot.troppca} does not return anything other than the plot.
#' @method plot troppca
#' @export
#' @export plot.troppca
plot.troppca <- function(x, plab = NULL, fw = FALSE, ...) {
  if (x$type == "linear space") {
    stop("Only principal component by tropical polytope is plottable.")
  }
  object <- x
  D <- eachrow(object$pc, object$pc[1, ], "-")[-1, ]
  if (fw){
    if (is.null(plab)) plab <- as.factor(c(rep(1, nrow(object$projection))))
    proj_points_plot <- t(apply(object$projection, 1, polytope_iso, D = object$pc))
    fw_point <- tropFW(proj_points_plot)
    proj_points_plot <- rbind(proj_points_plot, fw_point$fw)
    proj_2D_plot_m <- proj_points_plot - proj_points_plot[, 1]
  } else{
    if (is.null(plab)) plab <- as.factor(rep(1, nrow(object$projection)))
    proj_points_plot <- t(apply(object$projection, 1, polytope_iso, D = object$pc))
    proj_2D_plot_m <- proj_points_plot - proj_points_plot[, 1]
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(xpd = TRUE, mar = c(5.1, 4.1, 4.1, 5.1))
  k <- ncol(D)
  plot(D[1, ], D[2, ], xlab = "x1", ylab = "x2")
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      tseg1 <- tropseg(D[, i], D[, j])
      tseg2 <- tropseg(D[, i], D[, j], flag = TRUE)
      if (tseg1[[2]] < tseg2[[2]]) {
        tseg <- tseg1
      } else {
        tseg <- tseg2
      }
      segments(tseg[[1]][1, 1], tseg[[1]][2, 1], tseg[[1]][1, 2], tseg[[1]][2, 2], col = "black")
      segments(tseg[[1]][1, 2], tseg[[1]][2, 2], tseg[[1]][1, 3], tseg[[1]][2, 3], col = "black")
    }
  }
  for (i in 1:length(unique(plab))) {
    points(x = proj_2D_plot_m[plab == unique(plab)[i], 2], y = proj_2D_plot_m[plab == unique(plab)[i], 3], pch = 16, cex = .75, col = (i + 1))
  }
  coord <- par("usr")
  if (fw){
    plab <- c(plab, "FW")
    points(x = proj_2D_plot_m[nrow(proj_2D_plot_m), 2], y = proj_2D_plot_m[nrow(proj_2D_plot_m), 3], pch = 18, cex = 2, col = "black")
    legend(x = coord[2] * 1.05, y = coord[4], legend = unique(plab), pch = c(rep(16, length(unique(plab))-1), 18),
           col = c(2:(length(unique(plab))), "black"))
  }else{
    legend(x = coord[2] * 1.05, y = coord[4], legend = unique(plab), pch = rep(16, length(unique(plab))),
           col = c(2:(length(unique(plab)) + 1)))
  }
}
