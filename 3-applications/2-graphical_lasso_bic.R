rm(list=ls())

library(focus)
library(igraph)
library(colorspace)
library(RColorBrewer)
source("Scripts/additional_functions_specific_to_comparisons.R")

# Loading OMICS
cpg=readRDS('Data/NOWAC_MTT_smoking_159.rds')
ttx=readRDS('Data/NOWAC_TTX_smoking_208.rds')
ids=intersect(rownames(cpg), rownames(ttx))
cpg=cpg[ids,]
ttx=ttx[ids,]
omic=cbind(cpg, ttx)
pk=c(ncol(cpg), ncol(ttx))

# Loading Annotation
ttxannot=data.frame(readRDS('Data/NOWAC_TTX_annot_208.rds'))
cpgannot=readRDS('Data/NOWAC_MTT_annot_159.rds')
cpgannot=cbind(cpgannot[,c("alt.name"),drop=FALSE], name=rownames(cpgannot))
ttxannot=ttxannot[,c(3,2)]
colnames(ttxannot)=c("alt.name", "name")
omicannot=rbind(cpgannot, ttxannot)

# Calibration using the BIC
Lambda=LambdaGridNetwork(data=omic)
out_it=CalibrateInformationTheory(x=omic, Lambda=Lambda)
saveRDS(out_it, "Results/3-applications/graphical_model_BIC.rds")

# Adjacency matrix
adjacency=out_it$path[,,which.min(out_it$BIC)]
rownames(adjacency)=colnames(adjacency)=colnames(omic)
diag(adjacency)=0

# Number of block-specific edges
blockmat=GetBlockMatrix(pk=pk)
sum(adjacency[blockmat==1])/2 # cpg-cpg
sum(adjacency[blockmat==2])/2 # cpg-transcript
sum(adjacency[blockmat==3])/2 # transcript-transcript

# Preparing igraph object 
node_label=rownames(adjacency)
node_label=paste(omicannot$alt.name, "\n", omicannot$name)
g=Graph(adjacency=adjacency, node_colour=c(rep("skyblue",ncol(cpg)), rep("lightsalmon", ncol(ttx))), 
        node_label=node_label, node_shape=c(rep("square",ncol(cpg)), rep("circle", ncol(ttx))))

# Saving figure
{plotname=paste0("Figures/3-applications/Multi_omics_graph_BIC.pdf")
  pdf(plotname, height=15, width=10)
  layout(mat=matrix(c(1,2),ncol=1), height=c(1,2))
  par(mar=c(5,5,1,1))
  
  # Calibration plot
  plot(Lambda, out_it$BIC, pch=19, col="navy", cex=0.5,
       xlab=expression(lambda), ylab="BIC", cex.lab=1.5)
  mtext(text="A", side=2, line=3, at=max(out_it$BIC), las=1, cex=3)
  abline(v=axTicks(1),col="grey",lty=2)
  abline(h=axTicks(2),col="grey",lty=2)
  points(Lambda, out_it$BIC, pch=19, col="navy", cex=0.5)
  lines(Lambda, out_it$BIC, col="navy")
  abline(h=min(out_it$BIC),lty=3,col="red")
  abline(v=Lambda[which.min(out_it$BIC)],lty=3,col="red")
  
  # Calibrated graph
  set.seed(2)
  par(mar=rep(0,4))
  plot(g, layout=layout_with_fr(g), xlim=c(-1.2,1.2))
  mtext(text="B", side=2, line=-2, at=1, las=1, cex=3)
  legend("bottomright", legend=c("DNA methylation", "Gene expression"), 
         col=c("skyblue", "lightsalmon"), pch=c(15,19), bty="n", cex=1.2)
  # legend("left", pch=19, col=colors[sort(unique(chr_number))], legend=sort(unique(chr_number)), bty="n", cex=0.7, title="CHR")
  dev.off()
  system(paste("pdfcrop --margin 10",plotname,plotname))}

