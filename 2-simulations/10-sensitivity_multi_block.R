rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(colorspace)
library(MASS)
library(openxlsx)

# Simulation parameters
simul_study_id <- 1
topology <- "random"
simul_id <- 1
PFER_thr <- Inf
niter <- 1000

# Loading the simulation results
simul_id <- 1
subtable <- readRDS(paste0("Results/2-simulations/2-multi_block/Sensitivity_", simul_study_id, "_", topology, "/Performances_multi_sensitivity_", simul_id, "_merged_PFER_thr_", PFER_thr, ".rds"))
subtable <- subtable[, , 1:niter]

# Computing the median performances
tmpmedian <- apply(subtable, c(1, 2), FUN = function(x) {
  median(as.numeric(x))
})
mymedian <- tmpmedian
mymedian[, c("precision", "recall", "F1_score")] <- formatC(tmpmedian[, c("precision", "recall", "F1_score")], format = "f", digits = 3)
mymedian[, c("TP", "FP", "FN", "time")] <- formatC(tmpmedian[, c("TP", "FP", "FN", "time")], format = "f", digits = 0)
mymedian <- mymedian[, c("TP", "FP", "FN", "precision", "recall", "F1_score", "time")]

# Computing the inter-quartile range of performances
tmpiqr <- apply(subtable, c(1, 2), FUN = function(x) {
  IQR(as.numeric(x))
})
myiqr <- tmpiqr
myiqr[, c("precision", "recall", "F1_score")] <- formatC(tmpiqr[, c("precision", "recall", "F1_score")], format = "f", digits = 3)
myiqr[, c("TP", "FP", "FN", "time")] <- formatC(tmpiqr[, c("TP", "FP", "FN", "time")], format = "f", digits = 0)
myiqr <- myiqr[, c("TP", "FP", "FN", "precision", "recall", "F1_score", "time")]

# Combining the metrics
mytable <- matrix(paste0(mymedian, " [", myiqr, "]"), ncol = ncol(mymedian))
colnames(mytable) <- colnames(mymedian)

# Re-formatting the tables
# mytable <- rbind(mytable_single, mytable_multi)
mytable[rep((2:4), 8) + rep(seq(0, 28, by = 4), each = 3), "time"] <- ""
tmp <- rep("", nrow(mytable))
tmp[seq(1, 32, by = 4)] <- c("", "", c(0, 0.001, 0.01, 0.1, 0.5, 1))
mytable <- cbind(rep("", nrow(mytable)), tmp, rep("", nrow(mytable)), mytable)
mytable[1, 1] <- "Single"
mytable[2, 1] <- "block"
mytable[5, 1] <- "Multi"
mytable[6, 1] <- "parameters"
mytable[9, 1] <- "Multi"
mytable[10, 1] <- "block"
mytable[,3]=c(rep(c("Overall", "Within 1", "Between", "Within 2"), 8))

# Saving table
write.xlsx(mytable, file = "Tables/2-simulations/Sensitivity_multi_block_table.xlsx")
write.table(mytable,
  file = "Tables/2-simulations/Sensitivity_multi_block_table.txt",
  quote = FALSE, eol = "@@ \n", sep = "&", row.names = FALSE
)
