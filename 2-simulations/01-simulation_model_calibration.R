rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(sharp)
library(abind)

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
seed <- 1
u_list <- 10^-(seq(-1, 5, by = 0.1))

# Simulation
simul <- SimulateGraphical(
  n = n, pk = pk, topology = topology,
  v_within = 1, output_matrices = TRUE
)
omega <- simul$omega
diag(omega) <- 0
diag(omega) <- apply(abs(omega), 1, sum)
contrasts <- mycor <- NULL
for (k in 1:length(u_list)) {
  print(k)
  set.seed(seed)
  u <- u_list[k]
  tmpomega <- omega
  diag(tmpomega) <- diag(tmpomega) + u
  sigma <- cov2cor(solve(tmpomega))
  contrasts <- c(contrasts, Contrast(sigma))
  mycor <- abind(mycor, sigma, along = 3)
}

# Saving figure
plotname <- "Figures/2-simulations/Calibration_u_simulation_graphical_model.pdf"
{
  pdf(plotname, width = 12, height = 8)
  layout(mat = matrix(c(1, 1, 1, 2, 3, 4), ncol = 3, byrow = TRUE), heights = c(1.25, 1))
  par(mar = c(5, 5, 1, 5))

  # Contrast as a function of u
  par(xpd = FALSE)
  plot(log(u_list), contrasts,
    pch = 19, col = "navy", cex = 0.5,
    xlab = "log(u)", ylab = "Contrast", cex.lab = 1.5, las = 1, bty = "n"
  )
  abline(v = range(log(u_list)), lty = 2, col = "black")
  abline(v = axTicks(1), lty = 3, col = "grey")
  abline(h = axTicks(2), lty = 3, col = "grey")
  points(log(u_list), contrasts, pch = 19, col = "navy", cex = 0.5)
  points(log(u_list)[which.max(contrasts)], max(contrasts), pch = 19, col = "red")
  abline(v = log(u_list)[which.max(contrasts)], lty = 2, col = "red")
  abline(h = max(contrasts), lty = 2, col = "red")

  # Corresponding correlation matrices
  Heatmap(mycor[, , dim(mycor)[3]],
    col = c("darkblue", "white", "firebrick3"),
    legend_range = c(-1, 1), legend = FALSE, axes = FALSE
  )
  Heatmap(mycor[, , which.max(contrasts)],
    col = c("darkblue", "white", "firebrick3"),
    legend_range = c(-1, 1), legend = FALSE, axes = FALSE
  )
  Heatmap(mycor[, , 1],
    col = c("darkblue", "white", "firebrick3"),
    legend_range = c(-1, 1), legend_length = 10, legend = TRUE, axes = FALSE
  )
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))
