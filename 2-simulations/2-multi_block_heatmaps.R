rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(corpcor)


### Loading the data

# OMICS
cpg <- readRDS("Data/NOWAC_MTT_smoking_159.rds")
ttx <- readRDS("Data/NOWAC_TTX_smoking_208.rds")
ids <- intersect(rownames(cpg), rownames(ttx))
cpg <- cpg[ids, ]
ttx <- ttx[ids, ]

# Annotation
ttxannot <- readRDS("Data/NOWAC_TTX_smoking_208.rds")
cpgannot <- readRDS("Data/NOWAC_MTT_annot_159.rds")

# Extracting 200 observations, 50 random CpG sites and transcripts
set.seed(1)
ids <- sample(nrow(cpg), size = 200)
set.seed(1)
cpg <- cpg[ids, sample(ncol(cpg), size = 50)]
set.seed(1)
ttx <- ttx[ids, sample(ncol(ttx), size = 50)]


### Simulations

# Simulation parameters
n <- 200
pk <- c(50, 50)
topology <- "random"
nu <- 0.04
if (topology != "scale-free") {
  N <- (sum(pk) * (sum(pk) - 1)) / 2
  print(paste("Number of edges:", nu * N))
} else {
  print(paste("Number of edges:", p - 1))
}
seed <- 1
v_w <- 1
v_b <- 0.2

# Simulation of data with underlying graph structure
set.seed(seed)
simul <- SimulateGraphical(n = n, pk = pk, topology = topology, nu = nu, v_within = v_w, v_between = v_b)


### Estimation of (partial) correlations

# Saving figure
{
  pdf("Figures/2-simulations/Heatmaps_multi_block.pdf", width = 12, height = 11)
  par(mfrow = c(2, 2), mar = c(3, 3, 1, 5))

  # Real data
  mycor <- cor(cbind(cpg, ttx))
  mypcor <- cor2pcor(mycor)

  Heatmap(mycor,
    colours = c("darkblue", "white", "firebrick3"),
    legend_range = c(-1, 1), legend_length = 50, legend = FALSE, axes = FALSE
  )
  mtext(side = 2, at = nrow(mycor) * 1.025, text = "A", cex = 3, las = 1, line = 0.5)
  for (i in 1:2) {
    axis(side = i, at = c(0.5, 49.5), labels = NA)
    axis(
      side = i, at = mean(c(0.5, 49.5)), labels = ifelse(i == 1, yes = "DNA methylation", no = "Gene expression"),
      tick = FALSE, cex.axis = 1.5
    )
    axis(side = i, at = c(50.5, 99.5), labels = NA)
    axis(
      side = i, at = mean(c(50.5, 99.5)), labels = ifelse(i == 1, yes = "Gene expression", no = "DNA methylation"),
      tick = FALSE, cex.axis = 1.5
    )
  }

  Heatmap(mypcor,
    colours = c("darkblue", "white", "darkred"),
    legend_range = c(-1, 1), legend_length = 30, axes = FALSE
  )
  mtext(side = 2, at = nrow(mycor) * 1.025, text = "B", cex = 3, las = 1, line = 0.5)
  for (i in 1:2) {
    axis(side = i, at = c(0.5, 49.5), labels = NA)
    axis(
      side = i, at = mean(c(0.5, 49.5)), labels = ifelse(i == 1, yes = "DNA methylation", no = "Gene expression"),
      tick = FALSE, cex.axis = 1.5
    )
    axis(side = i, at = c(50.5, 99.5), labels = NA)
    axis(
      side = i, at = mean(c(50.5, 99.5)), labels = ifelse(i == 1, yes = "Gene expression", no = "DNA methylation"),
      tick = FALSE, cex.axis = 1.5
    )
  }

  # Simulations
  mycor <- cor(simul$data)
  mypcor <- cor2pcor(mycor)

  Heatmap(mycor,
    colours = c("darkblue", "white", "darkred"),
    legend_range = c(-1, 1), legend_length = 50, legend = FALSE, axes = FALSE
  )
  mtext(side = 2, at = nrow(mycor) * 1.025, text = "C", cex = 3, las = 1, line = 0.5)
  for (i in 1:2) {
    axis(side = i, at = c(0.5, 49.5), labels = NA)
    axis(
      side = i, at = mean(c(0.5, 49.5)), labels = ifelse(i == 1, yes = "Group 1", no = "Group 2"),
      tick = FALSE, cex.axis = 1.5
    )
    axis(side = i, at = c(50.5, 99.5), labels = NA)
    axis(
      side = i, at = mean(c(50.5, 99.5)), labels = ifelse(i == 1, yes = "Group 2", no = "Group 1"),
      tick = FALSE, cex.axis = 1.5
    )
  }

  Heatmap(mypcor,
    colours = c("darkblue", "white", "darkred"),
    legend_range = c(-1, 1), legend_length = 50, legend = FALSE, axes = FALSE
  )
  mtext(side = 2, at = nrow(mycor) * 1.025, text = "D", cex = 3, las = 1, line = 0.5)
  for (i in 1:2) {
    axis(side = i, at = c(0.5, 49.5), labels = NA)
    axis(
      side = i, at = mean(c(0.5, 49.5)), labels = ifelse(i == 1, yes = "Group 1", no = "Group 2"),
      tick = FALSE, cex.axis = 1.5
    )
    axis(side = i, at = c(50.5, 99.5), labels = NA)
    axis(
      side = i, at = mean(c(50.5, 99.5)), labels = ifelse(i == 1, yes = "Group 2", no = "Group 1"),
      tick = FALSE, cex.axis = 1.5
    )
  }
  dev.off()
}
