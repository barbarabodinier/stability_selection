rm(list = ls())

library(fake)
library(sharp)

# Parameters
n <- 100
pk <- 50
seed <- 5
PFER_thr <- 5
ev <- 0.7
q <- 10

# Simulation
set.seed(seed)
simul <- SimulateRegression(
  n = n, pk = pk,
  theta = cbind(c(rep(1, q), rep(0, pk - q))),
  beta_abs = c(0.7, 1), beta_sign = c(-1, 1), continuous = TRUE,
  family = "gaussian", ev_xy = ev
)
table(simul$theta)

# Stability selection for variable selection
system.time({
  out <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian", lambda.min.ratio = 1e-2)
})

# Extracting selection proportions
A <- SelectedVariables(out)
SelectionPerformance(A, simul$theta)
selprop <- SelectionProportions(out)

# Re-ordering the selection proportions
ids <- sort.list(selprop, decreasing = TRUE)
Asum <- A + 2 * simul$theta
selprop <- selprop[ids]
Asum <- Asum[ids]

# Incremental analyses
incremental <- Incremental(
  xdata = simul$xdata, ydata = simul$ydata,
  tau = 0.5,
  stability = out, n_predictors = 50
)

# Figure parameters
mylwd <- 3
mycolours <- c("grey", "tomato", "forestgreen", "navy") # FP: red, TP: navy, FN: forestgreen
names(mycolours) <- c("TN", "FP", "FN", "TP")
mar <- c(1, 5, 6, 10)

# Saving figure
plotname <- "Figures/1-illustrations/Figure1_variable_selection.pdf"
{
  pdf(plotname, width = 10, height = 14)
  # Calibration plot
  par(mfrow = c(4, 1), mar = mar)
  CalibrationPlot(out)
  mtext(text = "A", cex = 2.5, side = 3, at = -6, line = 0)

  # Selection proportions
  plot(selprop,
    type = "h", col = mycolours[Asum + 1], lend = 1,
    las = 1, xaxt = "n", ylim = c(0, 1),
    xlab = "", ylab = "Selection Proportion", cex.lab = 1.5, lwd = mylwd, bty = "n"
  )
  abline(h = Argmax(out)[2], lty = 2, col = "darkred")
  abline(h = 1 - Argmax(out)[2], lty = 2, col = "orange")
  axis(side = 1, at = 1:length(selprop), labels = NA)
  for (k in 1:length(selprop)) {
    axis(
      side = 1, at = k, labels = names(selprop)[k], las = 2,
      col.axis = (mycolours[Asum + 1])[k],
      cex.axis = 1, tick = FALSE
    )
  }
  legend("right",
    lty = c(rep(1, 4), 2, 2), lwd = c(rep(mylwd, 4), 1, 1), bg = "white", bty = "n",
    col = c(mycolours[c(1, 4, 2, 3)], "darkred", "orange"), cex = 1.2,
    legend = c(
      paste0("True Negatives (N=", sum(Asum == 0), ")"),
      paste0("True Positives (N=", sum(Asum == 3), ")"),
      paste0("False Positives (N=", sum(Asum == 1), ")"),
      paste0("False Negatives (N=", sum(Asum == 2), ")"),
      eval(parse(text = "expression(hat(pi))")),
      eval(parse(text = "expression(1-hat(pi))"))
    )
  )
  mtext(text = "B", cex = 2.5, side = 3, at = -4, line = 0)

  # Incremental plot
  PlotIncremental(incremental,
    col = mycolours[Asum + 1], bty = "n",
    ylab = expression(R^2), ylim = c(0, 0.8),
    cex = 1.75, sfrac = 0.002
  )
  mtext(text = "C", cex = 2.5, side = 3, at = -4, line = 0)
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))

# Constrained calibration
system.time({
  out_constr <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    family = "gaussian", PFER_thr = PFER_thr, lambda.min.ratio = 1e-2
  )
})

# Saving figure
{
  pdf("Figures/1-illustrations/Constrained_calibration.pdf", width = 10, height = 11)
  par(mfrow = c(3, 1), mar = c(1, 6, 4.2, 10))
  b <- 1
  stability <- out_constr
  mat <- stability$S_2d
  # ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
  ids <- 1:nrow(mat)
  mat <- mat[ids, , drop = FALSE]
  mat <- mat[, , drop = FALSE]
  colnames(mat) <- stability$params$pi_list
  rownames(mat) <- formatC(stability$Lambda[, b], format = "e", digits = 2)[ids]

  # Extracting corresponding numbers of selected variables (q)
  Q <- stability$Q[, b]
  Q <- Q[ids]

  # Heatmap representation
  # CalibrationPlot(out_constr)
  Heatmap(t(mat[nrow(mat):1, ncol(mat):1]), axes = FALSE)

  # Identifying best pair of parameters
  graphics::abline(
    h = which.min(abs(as.numeric(colnames(mat)) - Argmax(stability)[b, 2])) - 0.5,
    lty = 3, col = "red"
  )
  graphics::abline(
    v = nrow(mat) - which(stability$Lambda[ids, b] == Argmax(stability)[b, 1]) + 0.5,
    lty = 3, col = "red"
  )

  # Including axes
  graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = 2, labels = rev(rownames(mat)))
  graphics::axis(
    side = 2, at = (1:ncol(mat)) - 0.5, las = 2,
    labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 2)
  )
  graphics::axis(
    side = 3, at = (1:nrow(mat)) - 0.5, las = 2,
    labels = rev(formatC(Q, format = "f", big.mark = ",", digits = 0))
  )

  # Including axis labels
  graphics::mtext(text = expression(pi), side = 2, line = 3.5, cex = 1.5)
  graphics::mtext(text = expression(italic(q)), side = 3, line = 2.5, cex = 1.5)
  mtext(text = expression(lambda), cex = 1.5, side = 1, at = -3, line = 2.3)
  mtext(text = "A", cex = 3.5, side = 3, at = -7, line = 0)

  # Uncontrained curve
  par(xaxs = "i", yaxs = "i")
  mat <- out$S_2d
  vect <- rev(apply(mat, 1, max, na.rm = TRUE))
  vect[is.infinite(vect)] <- NA
  plot(seq(1, length(vect)) - 0.5, vect,
    cex.lab = 1.5, pch = 19, cex = 0.25, col = "navy",
    ylim = range(vect, na.rm = TRUE) + c(-50, 50),
    xaxt = "n", xlab = "", ylab = "Stability Score", xlim = c(0, length(vect)), bty = "n"
  )
  axis(side = 3, at = seq(1, length(vect)) - 0.5, labels = NA)
  axis(side = 1, at = c(1, length(vect)) - 0.5, labels = NA)
  lines(seq(1, length(vect)) - 0.5, vect, col = "navy")
  abline(v = length(vect) - ArgmaxId(out)[1] + 0.5, col = "navy", lty = 3)
  mtext(text = "B", cex = 3.5, side = 3, at = -7, line = 0.5)

  # Constrained curve
  mat <- out_constr$S_2d
  ids <- which(apply(mat, 1, FUN = function(x) {
    any(!is.na(x))
  }))
  vect <- apply(mat[ids, ], 1, max, na.rm = TRUE)
  points(nrow(mat) - ids + 0.5, vect, col = "red", pch = 19, cex = 0.25)
  lines(nrow(mat) - ids + 0.5, vect, col = "red", lty = 2)
  abline(v = nrow(mat) - ArgmaxId(out_constr)[1] + 0.5, col = "red", lty = 3)
  legend("topleft",
    lty = c(1, 2), col = c("navy", "red"), pch = 19, pt.cex = 0.25,
    bty = "n", cex = 1.2,
    legend = c(
      "Unconstrained calibration",
      paste0("Constrained calibration PFER<", PFER_thr)
    )
  )

  # Extracting and re-ordering the selection proportions
  A <- SelectedVariables(out_constr)
  selprop <- SelectionProportions(out_constr)
  ids <- sort.list(selprop, decreasing = TRUE)
  Asum <- A + 2 * simul$theta
  selprop <- selprop[ids]
  Asum <- Asum[ids]

  # Selection proportions
  par(xaxs = "r", yaxs = "r")
  par(mar = c(4.2, 6, 1, 7))
  plot(selprop,
    type = "h", col = mycolours[Asum + 1], las = 1, xaxt = "n", ylim = c(0, 1),
    xlab = "", ylab = "Selection Proportion", cex.lab = 1.5, lwd = mylwd, bty = "n"
  )
  abline(h = Argmax(out_constr)[2], lty = 2, col = "darkred")
  abline(h = 1 - Argmax(out_constr)[2], lty = 2, col = "orange")
  axis(side = 1, at = 1:length(selprop), labels = NA)
  for (k in 1:length(selprop)) {
    axis(side = 1, at = k, labels = names(selprop)[k], las = 2, cex.axis = 1, tick = FALSE)
  }
  legend("right",
    lty = c(rep(1, 4), 2, 2), lwd = mylwd, bg = "white", bty = "n",
    col = c(mycolours[c(1, 4, 2, 3)], "darkred", "orange"), cex = 1.2,
    legend = c(
      paste0("True Negatives (N=", sum(Asum == 0), ")"),
      paste0("True Positives (N=", sum(Asum == 3), ")"),
      paste0("False Positives (N=", sum(Asum == 1), ")"),
      paste0("False Negatives (N=", sum(Asum == 2), ")"),
      eval(parse(text = "expression(hat(pi))")),
      eval(parse(text = "expression(1-hat(pi))"))
    )
  )
  mtext(text = "C", cex = 3.5, side = 3, at = -4.3, line = -2.7)
  dev.off()
}
