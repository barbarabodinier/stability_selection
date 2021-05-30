rm(list = ls())

library(focus)
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


### Multi-block application

# Heatmap of correlations
mycor <- cor(omic)
plotname <- "Figures/3-applications/Multi_omics_correlation_heatmap.pdf"
{
  pdf(plotname, width = 8, height = 8)
  par(mar = c(5, 5, 5, 5))
  Heatmap(mycor,
    colours = c("darkblue", "white", "firebrick3"),
    legend_range = c(-1, 1), legend_length = 50, legend = FALSE, axes = FALSE
  )
  axis(side = 1, at = c(0, ncol(cpg)), labels = NA)
  axis(side = 1, at = mean(c(0, ncol(cpg))), labels = "DNA methylation", tick = FALSE, cex.axis = 1.5)
  axis(side = 1, at = c(ncol(cpg), ncol(omic)), labels = NA)
  axis(side = 1, at = mean(c(ncol(cpg), ncol(omic))), labels = "Gene expression", tick = FALSE, cex.axis = 1.5)
  axis(side = 2, at = c(0, ncol(ttx)), labels = NA)
  axis(side = 2, at = mean(c(0, ncol(ttx))), labels = "Gene expression", tick = FALSE, cex.axis = 1.5)
  axis(side = 2, at = c(ncol(ttx), ncol(omic)), labels = NA)
  axis(side = 2, at = mean(c(ncol(ttx), ncol(omic))), labels = "DNA methylation", tick = FALSE, cex.axis = 1.5)
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))

# Stability selection and calibration
PFER_thr <- 150
system.time({
  out <- GraphicalModel(
    xdata = omic, pk = pk, lambda_other_blocks = 0.1, Lambda_cardinal = 30,
    PFER_thr = PFER_thr, max_density = 0.2
  )
})
saveRDS(out, paste0("Results/3-applications/multi_block_PFER_", PFER_thr, ".rds"))
out <- readRDS(paste0("Results/3-applications/multi_block_PFER_", PFER_thr, ".rds"))

# Calibration plot
{
  pdf(paste0("Figures/3-applications/Calibration_multi_omics_graph_PFER_thr_", PFER_thr, ".pdf"), width = 15, height = 5)
  par(mfrow = c(1, 3), mar = c(7, 5, 7, 6))
  CalibrationPlot(out)
  dev.off()
}

# Numbers of block-specific edges
adjacency <- Adjacency(out)
blockmat <- BlockMatrix(pk = pk)
sum(adjacency[blockmat == 1]) / 2 # cpg-cpg
sum(adjacency[blockmat == 2]) / 2 # cpg-transcript
sum(adjacency[blockmat == 3]) / 2 # transcript-transcript

# Comparison with BIC results
out_it <- readRDS("Results/3-applications/graphical_model_BIC.rds")
adjacency_bic <- out_it$path[, , which.min(out_it$BIC)]
diag(adjacency_bic) <- 0
sum(adjacency_bic) / 2
sum(adjacency_bic[blockmat == 2]) / 2 # cpg-transcript
table(adjacency + 2 * adjacency_bic) / 2
table((adjacency + 2 * adjacency_bic)[blockmat == 1]) / 2 # cpg-cpg
table((adjacency + 2 * adjacency_bic)[blockmat == 2]) / 2 # cpg-transcript
table((adjacency + 2 * adjacency_bic)[blockmat == 3]) / 2 # transcript-transcript


### Graph coloured by platform

# Define colours
colors <- lighten(colorRampPalette(brewer.pal(12, name = "Paired"))(23), amount = 0.4)
chr_number <- as.character(omicannot$chr)
chr_number[chr_number == "X"] <- 23
chr_number[chr_number == "Y"] <- 24
chr_number <- as.numeric(chr_number)

# Make igraph object
node_label <- paste(omicannot$alt.name, "\n", omicannot$name)
mygraph <- Graph(
  adjacency = adjacency, node_colour = c(rep("skyblue", ncol(cpg)), rep("lightsalmon", ncol(ttx))),
  node_label = node_label, node_shape = c(rep("square", ncol(cpg)), rep("circle", ncol(ttx)))
)
V(mygraph)$size <- 0.7 * V(mygraph)$size

# Saving figure
myasp <- 0.7
myseed <- 0
{
  pdf(paste0("Figures/3-applications/Multi_omics_graph_platforms_PFER_", PFER_thr, ".pdf"),
    width = 14, height = myasp * 14
  )
  par(mar = rep(0, 4))
  set.seed(myseed)
  plot(mygraph, layout = layout_with_fr(mygraph), asp = myasp)
  dev.off()
}


### Graph coloured by chromosome

# Define colours
colours <- lighten(colorRampPalette(brewer.pal(12, name = "Paired"))(23), amount = 0.4)
chr_number <- as.character(omicannot$chr)
chr_number[chr_number == "X"] <- 23
chr_number[chr_number == "Y"] <- 24
chr_number <- as.numeric(chr_number)

# Make igraph object
node_label <- paste(omicannot$alt.name, "\n", omicannot$name)
mygraph <- Graph(
  adjacency = adjacency, node_colour = colours[chr_number],
  node_label = node_label, node_shape = c(rep("square", ncol(cpg)), rep("circle", ncol(ttx)))
)
V(mygraph)$size <- 0.7 * V(mygraph)$size

# Saving figure
myasp <- 0.7
myseed <- 0
{
  pdf(paste0("Figures/3-applications/Multi_omics_graph_chromosomes_PFER_", PFER_thr, ".pdf"), width = 14, height = myasp * 14)
  par(mar = rep(0, 4))
  set.seed(myseed)
  plot(mygraph, layout = layout_with_fr(mygraph), asp = myasp)
  tmp <- sort(unique(chr_number))
  tmp[tmp == "23"] <- "X"
  legend("topleft",
    pch = 19, col = colors[sort(unique(chr_number))],
    legend = tmp, bty = "n", pt.cex = 1.2, cex = 1, title = "CHR"
  )
  dev.off()
}


### Annotation using Reactome pathways

g <- mygraph
annot <- data.frame(fread("Data/panther_annotation_208.txt"))
annot <- annot[, -1]
colnames(annot) <- c(
  "MappedID", "Protein_name", "Panther_family", "Protein_class",
  "Pathway", "BioProcess", "MolFunction", "CellComponent",
  "Reactome", "GOBioProcess", "GOMolFunction", "GOCellComponent"
)
table(annot$Protein_class)

genegroup <- c("Reactome", "GOBioProcess", "GOMolFunction", "GOCellComponent")
panther <- annot
filepath <- "Figures/3-applications/"
dir.create(paste0(filepath, "annotation/"), showWarnings = FALSE)

# for (i in 1:length(genegroup)){
for (i in 1) {
  print(genegroup[i])
  filepath_group <- paste0(filepath, "annotation/", genegroup[i], "/")
  dir.create(filepath_group, showWarnings = FALSE)

  # Create the unique list of pathways involved
  mylist_to_explore <- NULL
  for (x in eval(parse(text = paste0("panther$", genegroup[i])))) {
    tmp <- unlist(strsplit(x, ";"))
    tmp <- tmp[(tmp != "") & (tmp != ";")]
    mylist_to_explore <- c(mylist_to_explore, tmp)
  }
  mylist_to_explore <- gsub("->.*", "", mylist_to_explore)
  mylist_to_explore <- unique(mylist_to_explore)

  # Number of pathways with given number of transcripts
  number_of_hits <- sapply(mylist_to_explore,
    FUN = function(x) {
      sum(unique(panther$MappedID[grep(x, eval(parse(text = paste0("panther$", genegroup[i]))), fixed = TRUE)]) %in% omicannot[V(g)$name, 1])
    }
  )

  # Choose the pathways to investigate further (i.e. if > 3 transcripts in it)
  mylist_to_investigate <- names(which(number_of_hits >= 5))

  for (k in 1:length(mylist_to_investigate)) {
    print(mylist_to_investigate[k])
    pathway_ttx <- panther$MappedID[grep(mylist_to_investigate[k], eval(parse(text = paste0("panther$", genegroup[i]))), fixed = TRUE)]
    {
      pdf(paste0(
        filepath_group, "Multi_omics_graph_PFER_", sum(PFER_thr), "_", genegroup[i], "_",
        gsub("/", "_", gsub(" ", "_", mylist_to_investigate[k])), ".pdf"
      ), width = 14, height = myasp * 14)
      par(mar = rep(0, 4))
      set.seed(myseed)
      V(g)$color <- ifelse(omicannot[V(g)$name, 1] %in% pathway_ttx, yes = "tomato", no = "grey90")
      V(g)$frame.color <- V(g)$color
      plot(g, layout = layout_with_fr(g), asp = myasp)
      dev.off()
    }
  }
  cat("\n")
}

# Extracting targeted pathways
k <- grep("Signal Transduction", mylist_to_investigate)
pathway_signal <- panther$MappedID[grep(mylist_to_investigate[k], eval(parse(text = paste0("panther$", genegroup[i]))), fixed = TRUE)]
k <- grep("Translation", mylist_to_investigate)[2]
pathway_translation <- panther$MappedID[grep(mylist_to_investigate[k], eval(parse(text = paste0("panther$", genegroup[i]))), fixed = TRUE)]
k <- grep("lipids", mylist_to_investigate)
pathway_lipids <- panther$MappedID[grep(mylist_to_investigate[k], eval(parse(text = paste0("panther$", genegroup[i]))), fixed = TRUE)]

# Make annotated igraph object
node_label <- paste(omicannot$alt.name, "\n", omicannot$name)
mynode_colours <- rep("grey", ncol(adjacency))
names(mynode_colours) <- omicannot[colnames(adjacency), 1]
pathway_signal <- intersect(pathway_signal, names(mynode_colours))
pathway_translation <- intersect(pathway_translation, names(mynode_colours))
pathway_lipids <- intersect(pathway_lipids, names(mynode_colours))
mynode_colours[pathway_signal] <- "tomato"
mynode_colours[pathway_translation] <- "seagreen"
mynode_colours[pathway_lipids] <- "darkgoldenrod1"
mynode_colours <- lighten(mynode_colours, amount = 0.05)
mynode_colours[1:ncol(cpg)] <- "skyblue"
mygraph <- Graph(
  adjacency = adjacency, node_colour = mynode_colours,
  node_label = node_label, node_shape = c(rep("square", ncol(cpg)), rep("circle", ncol(ttx)))
)
V(mygraph)$size <- 0.7 * V(mygraph)$size

# Saving figure
myasp <- 0.7
myseed <- 1
{
  pdf(paste0("Figures/3-applications/Multi_omics_graph_PFER_", PFER_thr, ".pdf"), width = 14, height = myasp * 14)
  par(mar = rep(0, 4))
  set.seed(myseed)
  plot(mygraph, layout = layout_with_fr(mygraph), asp = myasp)
  legend("bottomleft",
    pch = c(15, rep(19, 5)), pt.cex = 2, cex = 1.5,
    col = c("skyblue", "black", lighten(c("tomato", "seagreen", "darkgoldenrod1", "grey"), amount = 0.05)),
    legend = c("DNA methylation", "Gene expression", "Signal transduction", "Translation", "Metabolism of lipids", "Other"), bty = "n"
  )
  dev.off()
}
