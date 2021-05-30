rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(igraph)

# Simulation parameters
n <- 100
pk <- 50
topology <- "scale-free"
if (topology != "scale-free") {
  N <- (sum(pk) * (sum(pk) - 1)) / 2
  print(paste("Number of edges:", nu * N))
} else {
  print(paste("Number of edges:", sum(pk) - 1))
}
seed <- 2

# Simulation of data with underlying graph structure
set.seed(seed)
simul <- SimulateGraphical(n = n, pk = pk, topology = topology, output_matrices = TRUE)

# Stability selection graphical LASSO
out <- GraphicalModel(xdata = simul$data, PFER_thr = 20, PFER_method = "MB")
A <- Adjacency(out)

# Evaluating the selection performances
perf <- SelectionPerformance(A, simul$theta)
perf
mygraph <- SelectionPerformanceGraph(A, simul$theta)

# Saving figure
{
  pdf("Figures/1-illustrations/Figure1_graphical_model.pdf", width = 15, height = 7)
  par(mfrow = c(1, 2), mar = c(7, 5, 5, 5))

  # Calibration plot
  CalibrationPlot(out)
  mtext(text = "C", cex = 4, side = 3, at = -9, line = 0)

  # Graph representation of the detected and missed edges
  par(mar = c(2, 2, 2, 0))
  set.seed(2)
  plot(mygraph, layout = layout_with_fr(mygraph))
  mtext(text = "D", cex = 4, side = 3, at = -1, line = -3)
  legend(
    x = -1.2, y = -0.5, col = c("navy", "tomato", "forestgreen"),
    lty = c(1, 2, 3), bty = "n", cex = 1.2,
    legend = c(
      paste0("True Positives (N=", perf$TP, ")"),
      paste0("False Positives (N=", perf$FP, ")"),
      paste0("False Negatives (N=", perf$FN, ")")
    )
  )
  dev.off()
}
