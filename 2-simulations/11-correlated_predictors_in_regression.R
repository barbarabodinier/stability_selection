rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(sharp)

# Reading arguments
seed <- 1
PFER_thr <- 5
pi_list <- seq(0.6, 0.9, by = 0.05)
pk <- rep(100, 10)
nu_xz <- 0.02

nu_within_list <- c(0, 0.02)
for (k in 1:length(nu_within_list)) {
  nu_within <- nu_within_list[k]
  print(nu_within)

  # Simulation
  set.seed(seed)
  adjacency <- SimulateAdjacency(
    pk = pk,
    nu_within = nu_within,
    nu_between = 0
  )
  theta <- sharp:::SamplePredictors(pk = pk, q = length(pk), nu = nu_xz, orthogonal = TRUE)
  theta <- apply(theta, 1, sum)
  simul <- SimulateRegression(
    n = 500, pk = sum(pk),
    adjacency_x = adjacency,
    theta_xz = cbind(theta),
    v_within = 1,
    eta_set = 1,
    family = "gaussian", ev_xz = 0.4
  )
  print(table(simul$theta))

  # Heatmap of correlations between predictors
  mycor <- cor(simul$xdata)
  plotname <- paste0("Figures/2-simulations/Heatmap_correlation_predictors_in_regression_", k, ".pdf")
  {
    pdf(plotname, width = 8, height = 8)
    par(mar = c(5, 5, 5, 7))
    Heatmap(mycor,
      col = c("darkblue", "white", "firebrick3"),
      legend_range = c(-1, 1), legend_length = 300, legend = TRUE, axes = FALSE
    )
    dev.off()
  }
}
