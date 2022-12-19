rm(list = ls())

library(sharp)
library(igraph)
library(colorspace)
library(RColorBrewer)
library(data.table)


### Loading the data

# OMICS
cpg <- readRDS("Data/NOWAC_MTT_smoking_159.rds")
ttx <- readRDS("Data/NOWAC_TTX_smoking_208.rds")
ids <- intersect(rownames(cpg), rownames(ttx))
cpg <- cpg[ids, ]
ttx <- ttx[ids, ]
omic <- cbind(cpg, ttx)
pk <- c(ncol(cpg), ncol(ttx))

# Annotation
ttxannot <- data.frame(readRDS("Data/NOWAC_TTX_annot_208.rds"))
cpgannot <- readRDS("Data/NOWAC_MTT_annot_159.rds")
cpgannot <- cbind(cpgannot[, c("alt.name", "chr"), drop = FALSE], name = rownames(cpgannot))
ttxannot <- ttxannot[, c(3, 4, 2)]
colnames(ttxannot) <- c("alt.name", "chr", "name")
omicannot <- rbind(cpgannot, ttxannot)

# Covariates
covars <- readRDS("Data/Original/full_covariates_nowac.rds")
covars <- covars[!is.na(covars$sampleID), ]
rownames(covars) <- covars$sampleID
covars <- covars[, c("case_ctrl", "age.sample", "BMI", "smoking.status", "packyrsTF")]
covars$case <- ifelse(covars$case_ctrl == "case", yes = 1, no = 0)
covars <- covars[!is.na(covars$BMI), ]

# Merging the datasets
ids <- intersect(rownames(covars), rownames(cpg))
covars <- covars[ids, ]
cpg <- cpg[ids, ]
ttx <- ttx[ids, ]


### Logistic regression

# Stability selection
out <- VariableSelection(xdata = cpg, ydata = covars$case, family = "binomial")

# Calibration plot
plotname <- "Figures/3-applications/Calibration_methylation_lung_cancer_lasso.pdf"
{
  pdf(plotname,
    width = 8, height = 8
  )
  par(mar = c(7, 6, 7, 6))
  CalibrationPlot(out)
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))

# Incremental analyses
K_expl <- 100
tau_expl <- 0.7
n_predictors <- 50
incremental <- Incremental(
  xdata = cpg, ydata = covars$case,
  stability = out,
  tau = tau_expl, K = K_expl,
  n_predictors = n_predictors
)

# Explanatory performance of smoking (packyears)
expl_smoking <- ExplanatoryPerformance(
  xdata = covars[, "packyrsTF", drop = FALSE],
  ydata = covars$case,
  family = "binomial",
  tau = tau_expl, K = K_expl
)

# Figure parameters
mylwd <- 1.5
mycolours <- c(
  rep("navy", sum(SelectedVariables(out))),
  rep("grey", n_predictors - sum(SelectedVariables(out)))
)
mar <- c(3, 5, 0, 1)

# Saving figure
plotname <- "Figures/3-applications/Methylation_lung_cancer_lasso.pdf"
{
  pdf(plotname, width = 8, height = 8)
  par(mfrow = c(3, 1), mar = mar)
  # Selection proportions
  plot(NA,
    bty = "n", xlim = c(1, ncol(cpg)), ylim = c(0, 1),
    xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  abline(h = seq(0, 1, by = 0.2), lty = 3, col = "grey70")
  abline(h = Argmax(out)[2], lty = 2, col = "darkred")
  abline(v = c(1, n_predictors), lty = 3)
  par(new = TRUE)
  selprop <- sort(SelectionProportions(out), decreasing = TRUE)
  plot(selprop,
    type = "h", col = c(
      rep("navy", sum(SelectedVariables(out))),
      rep("grey", ncol(cpg) - sum(SelectedVariables(out)))
    ),
    lend = 1, las = 1, xaxt = "n", ylim = c(0, 1),
    xlab = "", ylab = "Selection Proportion", cex.lab = 1.5, lwd = mylwd, bty = "n"
  )
  axis(side = 1, at = 1:length(selprop), labels = NA)

  # Incremental plot
  plot(NA,
    bty = "n", xlim = c(1, n_predictors), ylim = c(0.45, 0.85),
    xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  abline(h = seq(0.5, 0.8, by = 0.1), lty = 3, col = "grey70")
  abline(h = median(expl_smoking$AUC), lty = 2, col = "darkolivegreen")
  par(new = TRUE)
  PlotIncremental(incremental,
    ylab = "AUC",
    ylim = c(0.45, 0.85), bty = "n",
    col = mycolours, cex = 1.25, sfrac = 0.0025
  )
  text(
    x = 1, y = median(expl_smoking$AUC), adj = c(0, -1),
    label = "Packyears",
    col = "darkolivegreen"
  )
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))
