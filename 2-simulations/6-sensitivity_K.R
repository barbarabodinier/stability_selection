rm(list=ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(plotrix)
library(colorspace)

# Simulation parameters
simul_study_id=1
topology="random"
# topology="scale-free"
PFER_thr=30

# Template design
pi_list=seq(0.6,0.9,by=0.05)
mycolours=c(colorRampPalette(c(lighten("royalblue", amount=0.6), darken("navy", amount=0.6)))(9))
mypch=18
dimensionality=c("Low", "Intermediate", "High")

# Saving figure
metric="F1_score"
{pdf(paste0("Figures/2-simulations/Sensitivity_K_", metric, "_", simul_study_id, "_", topology, ".pdf"), 
     width=12, height=4)
  par(mar=c(5,5,3,1), mfrow=c(1,3))
  for (simul_id in 1:3){
    performances=readRDS(paste0("Results/2-simulations/1-graphical_model/Sensitivity_K_",simul_study_id,"_",topology,"/Performances_",simul_id,"_merged.rds"))
    performances=performances[,,1:min(dim(performances)[3],1000)]
    
    # Extracting the performance metric
    mylist=list()
    for (k in 1:nrow(performances)){
      mylist=c(mylist, list(as.numeric(performances[k,metric,])))
      assign(paste0("median",k), median(as.numeric(performances[k,metric,])))
    }
    
    # Extracting the computation times
    mylist_time=list()
    for (k in 1:nrow(performances)){
      mylist_time=c(mylist_time, list(as.numeric(performances[k,"time.user.self",])))
      assign(paste0("median",k), median(as.numeric(performances[k,"time.user.self",])))
    }
    
    # Figure parameters
    xseq=c(1:9)+(simul_id-1)*10
    myylab=ifelse(simul_id==1, yes=metric, no="")
    if (myylab=="F1_score"){
      myylab=eval(parse(text="expression(F[1]*'-score')"))
    }
    
    # Making boxplots
    plotCI(x=log(sapply(mylist_time,median)), y=sapply(mylist,median), ui=sapply(mylist, quantile, probs=0.95), li=sapply(mylist, quantile, probs=0.05),
           col=mycolours, pch=18, xlab="Time (s)", ylab=myylab, cex.lab=1.5, cex=1, sfrac=0.005, las=1, 
           xlim=c(1, log(max(sapply(mylist_time, quantile, probs=0.95)))), ylim=c(0,1), xaxt="n", main=c("Low", "Intermediate", "High")[simul_id], cex.main=1.5)
    mtext(text=LETTERS[simul_id], side=2, at=1, line=3, cex=2, las=1)
    xticks=c(1, 10, 100, 1000, 10000, 100000)
    axis(side=1, at=log(xticks), labels=xticks)
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
    lines(x=log(sapply(mylist_time,median)), y=sapply(mylist,median), col="grey")
    plotCI(x=log(sapply(mylist_time,median)), y=sapply(mylist,median), ui=sapply(mylist, quantile, probs=0.95), li=sapply(mylist, quantile, probs=0.05),
           col=mycolours, pch=18, xlab="", ylab=myylab, cex.lab=1.5, cex=1, sfrac=0.005, las=1, 
           xlim=c(1, log(max(sapply(mylist_time, quantile, probs=0.95)))), add=TRUE)
    plotCI(x=log(sapply(mylist_time,median)), y=sapply(mylist,median), ui=log(sapply(mylist_time, quantile, probs=0.95)), li=log(sapply(mylist_time, quantile, probs=0.05)),
           col=mycolours, pch=18, ylab=myylab, cex.lab=1.5, cex=1, sfrac=0.005, las=1, add=TRUE, err="x")
    if (simul_id==1){
      legend("bottomright", col=mycolours, pch=18, lty=1, title="K", bty="n",
             legend=formatC(c(10, 20, 50, 100, 500, 1000, 2000, 5000, 10000), 
                            format="f", big.mark=",", digits=0))
    }
  }
  dev.off()}


