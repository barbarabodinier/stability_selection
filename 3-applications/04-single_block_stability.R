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


### Single-block application

# Stability selection and calibration
PFER_thr <- 70
system.time({
  out <- GraphicalModel(
    xdata = cpg, 
    PFER_thr = PFER_thr, max_density = 0.2
  )
})
saveRDS(out, paste0("Results/3-applications/methylation_PFER_", PFER_thr, ".rds"))
out <- readRDS(paste0("Results/3-applications/methylation_PFER_", PFER_thr, ".rds"))

# Calibration plot
{
  pdf(paste0("Figures/3-applications/Calibration_methylation_graph_PFER_thr_", PFER_thr, ".pdf"), width = 7, height = 7)
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(out)
  dev.off()
}

# Numbers of block-specific edges
adjacency <- Adjacency(out)


### Graph

# Make igraph object
node_label <- paste(cpgannot$alt.name, "\n", cpgannot$name)
mygraph <- Graph(
  adjacency = adjacency, node_colour = lighten("skyblue", amount = 0.4),
  node_label = node_label, node_shape = "square")
V(mygraph)$size <- 0.7 * V(mygraph)$size

# Saving figure
myasp <- 1
myseed <- 1
{
  pdf(paste0("Figures/3-applications/Methylation_graph_platforms_PFER_", PFER_thr, ".pdf"),
      width = 14, height = myasp * 14
  )
  par(mar = rep(0, 4))
  set.seed(myseed)
  plot(mygraph, layout = layout_with_fr(mygraph), asp = myasp)
  dev.off()
}

