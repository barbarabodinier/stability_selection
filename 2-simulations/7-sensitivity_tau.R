rm(list=ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(plotrix)
library(colorspace)

# Simulation parameters
simul_study_id=1
topology="random"
PFER_thr=30

# Template design
mycolours=c(colorRampPalette(c("darkred", lighten("red", amount=0.5)))(9), "navy", "darkolivegreen")
dimensionality=c("Low", "Intermediate", "High")

# Saving figure
ylabs=c("Precision", "Recall", "F1_score")
names(ylabs)=c("precision", "recall", "F1_score")
plotname=paste0("Figures/2-simulations/Sensitivity_tau_", simul_study_id, "_", topology, ".pdf")
{pdf(plotname, width=8, height=12)
  par(mar=c(0,5,2,1), mfrow=c(5,1))
  plot.new()
  for (metric_id in 1:length(ylabs)){
    metric=names(ylabs)[metric_id]
    for (simul_id in 1:3){
      performances=readRDS(paste0("Results/2-simulations/1-graphical_model/Sensitivity_tau_",simul_study_id,"_",topology,"/Performances_",simul_id,"_merged.rds"))
      mylist=list()
      for (k in 1:nrow(performances)){
        mylist=c(mylist, list(as.numeric(performances[k,metric,])))
        assign(paste0("median",k), median(as.numeric(performances[k,metric,])))
      }
      nrows=nrow(performances)
      xseq=c(1:nrows)+(simul_id-1)*(nrows+1)
      myylab=ylabs[metric]
      if (myylab=="F1_score"){
        myylab=eval(parse(text="expression(F[1]*'-score')"))
      }
      if (simul_id==1){
        plotCI(x=xseq, y=sapply(mylist,median), ui=sapply(mylist, quantile, probs=0.95), li=sapply(mylist, quantile, probs=0.05),
               col=mycolours, ylim=c(0,1), xlim=c(1,3*nrows+2), pch=18, xaxt="n", xlab="", ylab=myylab, cex.lab=1.5, cex=1.5, sfrac=0.005, las=1)
        abline(v=c(1:nrows,(nrows+2):(2*nrows+1),(2*nrows+3):(3*nrows+2)), lty=3, col=mycolours)
        abline(v=c(nrows+1,2*nrows+2), lty=2)
        abline(h=axTicks(2), lty=3, col="grey", lwd=0.5)
        mtext(text=LETTERS[metric_id], side=2, at=1, line=2.7, las=1, cex=2)
      }
      if (metric=="precision"){
        par(xpd=TRUE)
        text(x=mean(xseq), y=1.1, labels=c("Low", "Intermediate", "High")[simul_id], cex=1.5)
        par(xpd=FALSE)
      }
      plotCI(x=xseq, y=sapply(mylist,median), ui=sapply(mylist, quantile, probs=0.95), li=sapply(mylist, quantile, probs=0.05),
             col=mycolours, pch=18, add=TRUE, cex=1.5, sfrac=0.005)
      if (metric=="F1_score"){
        axis(side=1, at=xseq[(1:nrows)], labels=c(seq(0.1,0.9,by=0.1), "CPSS", "Bootstrap"), las=2)
        axis(side=1, at=xseq[c(2,6)], line=3.5, labels=NA)
        axis(side=1, at=mean(xseq[c(2,6)]), line=3.5, labels="Subsampling", tick=FALSE)
      }
    }
  }
  dev.off()}
system(paste("pdfcrop --margin 10",plotname,plotname))


