rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(colorspace)
library(openxlsx)

# Simulation parameters
simul_study_id <- 1
topology <- "random"
# topology="scale-free"
PFER_thr <- 50

# Template design
pi_list <- seq(0.6, 0.9, by = 0.05)
pal_mb <- colorRampPalette(c("skyblue", "navy"))(length(pi_list))
pal_ss <- colorRampPalette(c("darkolivegreen2", "forestgreen"))(length(pi_list))
mycolours <- c("chocolate", "sienna", "pink", pal_mb, pal_ss, "darkred", "firebrick", "red", "tomato")
mypch <- c(rep(17, length(pi_list)), rep(17, length(pi_list)), 18, 18, 18, 18)
dimensionality <- c("Low", "Intermediate", "High")

# Saving table
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/2-simulations/1-graphical_model/Simulations_", simul_study_id, "_", topology, "/Performances_", simul_id, "_merged_PFER_thr_", PFER_thr, ".rds"))

  tmpmedian <- apply(performances, c(1, 2), FUN = function(x) {
    median(as.numeric(x))
  })
  mymedian <- tmpmedian
  mymedian[, c("precision", "recall", "F1_score")] <- formatC(tmpmedian[, c("precision", "recall", "F1_score")], format = "f", digits = 3)
  mymedian[, c("TP", "FP", "FN", "time")] <- formatC(tmpmedian[, c("TP", "FP", "FN", "time")], format = "f", digits = 0)
  mymedian <- mymedian[, c("pi", "TP", "FP", "FN", "precision", "recall", "F1_score", "time")]

  tmpiqr <- apply(performances[, -1, ], c(1, 2), FUN = function(x) {
    IQR(as.numeric(x))
  })
  myiqr <- tmpiqr
  myiqr[, c("precision", "recall", "F1_score")] <- formatC(tmpiqr[, c("precision", "recall", "F1_score")], format = "f", digits = 3)
  myiqr[, c("TP", "FP", "FN", "time")] <- formatC(tmpiqr[, c("TP", "FP", "FN", "time")], format = "f", digits = 0)
  myiqr <- myiqr[, c("TP", "FP", "FN", "precision", "recall", "F1_score", "time")]

  mytable <- matrix(paste0(mymedian[, -1], " [", myiqr, "]"), ncol = ncol(mymedian) - 1)
  mytable <- cbind(mymedian[, 1], mytable)
  colnames(mytable) <- colnames(mymedian)
  assign(paste0("mytable", simul_id), mytable)
}
mytable <- rbind(mytable1, mytable2, mytable3)
mytable <- cbind(rep(c("AIC", "BIC", "EBIC", "StARS", rep("", 3), "MB", rep("", 6), "SS", rep("", 3), "Subsampling", "CPSS", "MB", "SS"), 3), mytable)
mytable[is.na(mytable)] <- ""
mytable <- cbind(c(rep("", 10), "Low", rep("", 11), rep("", 10), "Intermediate", rep("", 11), rep("", 10), "High", rep("", 11)), mytable)

write.xlsx(mytable, paste0("Tables/2-simulations/Table_precision_recall_", simul_study_id, "_", topology, "_PFER_thr_", PFER_thr, ".xlsx"),
  row.names = FALSE, col.names = TRUE
)
write.table(mytable, paste0("Tables/2-simulations/Table_precision_recall_", simul_study_id, "_", topology, "_PFER_thr_", PFER_thr, ".txt"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)

# Saving figure
for (metric in c("F1_score")) {{ pdf(paste0("Figures/2-simulations/Boxplot_", metric, "_", simul_study_id, "_", topology, "_PFER_thr_", PFER_thr, ".pdf"),
  width = 12, height = 4.5
)
par(mar = c(7, 5, 3, 1), mfrow = c(1, 3))
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/2-simulations/1-graphical_model/Simulations_", simul_study_id, "_", topology, "/Performances_", simul_id, "_merged_PFER_thr_", PFER_thr, ".rds"))
  mylist <- list()
  for (k in 2:nrow(performances)) {
    mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
    assign(paste0("median", k), median(as.numeric(performances[k, metric, ])))
  }
  xseq <- c(1:7, 9:15, 17:20) + 4
  xseq <- c(1:3, xseq)
  myylab <- ifelse(simul_id == 1, yes = metric, no = "")
  if (myylab == "F1_score") {
    myylab <- eval(parse(text = "expression(F[1]*'-score')"))
  }
  boxplot(
    at = xseq, mylist, col = mycolours, boxcol = "white", whiskcol = mycolours, staplecol = mycolours,
    whisklty = 1, range = 0, las = 2, main = dimensionality[simul_id], cex.main = 1.5,
    ylab = myylab, cex.lab = 1.5, xaxt = "n", ylim = c(0, 1)
  )
  mtext(text = LETTERS[simul_id], side = 2, at = 1.1, line = 3, cex = 2, las = 1)
  abline(h = axTicks(2), lty = 3, col = "grey")
  boxplot(
    at = xseq, mylist, col = mycolours, boxcol = "white", whiskcol = mycolours, staplecol = mycolours,
    whisklty = 1, range = 0, las = 2, add = TRUE,
    ylab = myylab, cex.lab = 1.5, xaxt = "n"
  )
  abline(h = median19, col = "darkred", lty = 2)
  abline(h = median20, col = "firebrick", lty = 2)
  abline(h = median21, col = "red", lty = 2)
  abline(h = median22, col = "tomato", lty = 2)
  axis(side = 1, at = xseq[(1:3)], labels = c("BIC", "EBIC", "StARS"), las = 2)
  axis(side = 1, at = xseq[(1:(length(pi_list) * 2)) + 3], labels = rep(pi_list, 2), las = 2)
  axis(side = 1, at = xseq[((length(xseq) - 3):length(xseq))], labels = NA)
  axis(
    side = 1, at = xseq[((length(xseq) - 3):length(xseq))], las = 2, line = 0, tick = FALSE, # hadj=0.5,
    labels = c(
      "Subsampling", "CPSS",
      eval(parse(text = paste0("expression(PFER[MB]<", PFER_thr, ")"))),
      eval(parse(text = paste0("expression(PFER[SS]<", PFER_thr, ")")))
    )
  )
  axis(side = 1, at = c(1, 7) + 4, labels = NA, las = 2, line = 3.5)
  axis(
    side = 1, at = mean(c(1, 7)) + 4, labels = eval(parse(text = paste0("expression(PFER[MB]<", PFER_thr, ")"))),
    las = 1, line = 3.5, tick = FALSE
  )
  axis(side = 1, at = c(9, 15) + 4, labels = NA, las = 2, line = 3.5)
  axis(
    side = 1, at = mean(c(9, 15)) + 4, labels = eval(parse(text = paste0("expression(PFER[SS]<", PFER_thr, ")"))),
    las = 1, line = 3.5, tick = FALSE
  )
  abline(v = c(4, 12, 20), lty = 3, col = "black")
}
dev.off() }}
