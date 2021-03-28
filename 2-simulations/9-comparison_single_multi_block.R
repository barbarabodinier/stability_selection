rm(list=ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(colorspace)
library(MASS)
library(openxlsx)

# Simulation parameters
simul_study_id=1
topology="random"
simul_id=1
PFER_thr=Inf
niter=1000

# Saving table
for (simul_id in 1:3){
  print(simul_id)
  
  # Loading the simulation results 
  perf=readRDS(paste0("Results/2-simulations/3-multi_block/Simulations_",simul_study_id,"_",topology,"/Performances_single_",simul_id,"_merged_PFER_thr_",PFER_thr,".rds"))
  perf_block=readRDS(paste0("Results/2-simulations/3-multi_block/Simulations_",simul_study_id,"_",topology,"/Performances_multi_",simul_id,"_merged_PFER_thr_",PFER_thr,".rds"))
  
  # Summarising the simulation results
  mynames=c("perf", "perf_block")
  names(mynames)=c("single", "multi")
  for (k in 1:2){
    subtable=eval(parse(text=mynames[k]))
    
    # Computing the median performances
    tmpmedian=apply(subtable, c(1,2), FUN=function(x){median(as.numeric(x))})
    mymedian=tmpmedian
    mymedian[,c("precision", "recall", "F1_score")]=formatC(tmpmedian[,c("precision", "recall", "F1_score")], format="f", digits=3)
    mymedian[,c("TP", "FP", "FN", "time")]=formatC(tmpmedian[,c("TP", "FP", "FN", "time")], format="f", digits=0)
    mymedian=mymedian[,c("TP", "FP", "FN", "precision", "recall", "F1_score", "time")]
    
    # Computing the inter-quartile range of performances
    tmpiqr=apply(subtable, c(1,2), FUN=function(x){IQR(as.numeric(x))})
    myiqr=tmpiqr
    myiqr[,c("precision", "recall", "F1_score")]=formatC(tmpiqr[,c("precision", "recall", "F1_score")], format="f", digits=3)
    myiqr[,c("TP", "FP", "FN", "time")]=formatC(tmpiqr[,c("TP", "FP", "FN", "time")], format="f", digits=0)
    myiqr=myiqr[,c("TP", "FP", "FN", "precision", "recall", "F1_score", "time")]
    
    # Combining the metrics
    mytable=matrix(paste0(mymedian, " [", myiqr, "]"), ncol=ncol(mymedian))
    colnames(mytable)=colnames(mymedian)
    
    # Storing the output
    assign(paste0("mytable_",names(mynames)[k]), mytable)
  }
  
  mytable=rbind(mytable_single, mytable_multi)
  mytable[c(2:4,6:8),"time"]=""
  assign(paste0("mytable",simul_id), mytable)
}
mytable=rbind(mytable1, mytable2, mytable3)
mytable=cbind(c(rep("",3),"Low",rep("",7),"Intermediate",rep("",7),"High",rep("",4)),
              rep(c("Single-block",rep("",3),"Multi-block",rep("",3)),3),
              rep(c("Overall","Within-block 1","Between-block","Within-block 2"),6),mytable)
write.xlsx(mytable, file="Tables/2-simulations/Single_multi_block_table.xlsx")
write.table(mytable, file="Tables/2-simulations/Single_multi_block_table.txt", 
            quote=FALSE, eol="@@ \n", sep="&", row.names=FALSE)

# Template design
block_name=c("Within-block 1", "Between-blocks", "Within-block 2")
mybandwidth=0.3
myngrid=50
mylwd=0.5

# Saving figure
plotname=paste0("Figures/2-simulations/Precision_recall_contour_",simul_study_id,"_",topology,"_",simul_id,"_PFER_thr_",PFER_thr,".pdf")
{pdf(plotname, width=11, height=12)
  par(mar=c(0,5,3,1), mfrow=c(4,3))
  for (simul_id in 1:3){
    print(simul_id)
    
    # Loading the simulation results 
    perf=readRDS(paste0("Results/2-simulations/3-multi_block/Simulations_",simul_study_id,"_",topology,"/Performances_single_",simul_id,"_merged_PFER_thr_",PFER_thr,".rds"))
    perf_block=readRDS(paste0("Results/2-simulations/3-multi_block/Simulations_",simul_study_id,"_",topology,"/Performances_multi_",simul_id,"_merged_PFER_thr_",PFER_thr,".rds"))
    perf=perf[,,1:niter]
    perf_block=perf_block[,,1:min(dim(perf)[3], dim(perf_block)[3])]
    
    # For loop over the different blocks
    for (k in 1:3){
      # Contour plot
      z <- kde2d(perf[k+1,"precision",], perf[k+1,"recall",], 
                 h=mybandwidth, n=myngrid, lims = c(0,1,0,1))
      contour(z, las=1, cex.main=1.5, xlim=c(0,1), ylim=c(0,1),
              main=ifelse(simul_id==1, yes=block_name[k], no=""), 
              xlab="", ylab="", bty="n", 
              panel.first=c(abline(h=seq(0,1,by=0.2), lty=3, col="grey"),
                            abline(v=seq(0,1,by=0.2), lty=3, col="grey")),
              lwd=mylwd, cex.lab=1.5, col=lighten("navy", amount=0),
              pch=19, cex=0.5, las=1, lty=1, drawlabels=FALSE)
      
      # Axis labels
      if (k==1){
        mtext(LETTERS[simul_id], side=2, at=1, las=1, cex=2, line=3)
        mtext(text="Recall", side=2, line=3)
      }
      if (simul_id==3){
        par(xpd=TRUE)
        mtext(text="Precision", side=1, line=3)
        par(xpd=FALSE)
      }
      
      # Simulation points
      points(perf[k+1,"precision",], perf[k+1,"recall",], 
             col=lighten("skyblue", amount=0.1), pch=19, cex=0.5)
      z1=z
      if ((length(unique(perf_block[k+1,"recall",]))>1)&(length(unique(perf_block[k+1,"precision",]))>1)){
        z <- kde2d(perf_block[k+1,"precision",], perf_block[k+1,"recall",], 
                   h=mybandwidth, n=myngrid)
        contour(z, las=1,main=block_name[k], cex.main=2, add = TRUE,
                lwd=mylwd, cex.lab=1.5, col=lighten("red", amount=0),
                pch=17, cex=0.5, las=1, lty=5, drawlabels=FALSE)
      }
      points(perf_block[k+1,"precision",], perf_block[k+1,"recall",], 
             col=lighten("tomato", amount=0.3), pch=17, cex=0.7)
      contour(z1, las=1, main=block_name[k], cex.main=2, xlim=c(0,1), ylim=c(0,1),
              xlab=ifelse(k==1, yes="Precision", no=""),
              ylab=ifelse(k==1, yes="Recall", no=""),
              lwd=mylwd, cex.lab=1.5, col=lighten("navy", amount=0),
              pch=19, cex=0.5, las=1, lty=1, add=TRUE, drawlabels=FALSE)
      contour(z, las=1,main=block_name[k], cex.main=2, add = TRUE,
              lwd=mylwd, cex.lab=1.5, col=lighten("red", amount=0),
              pch=17, cex=0.5, las=1, lty=5, drawlabels=FALSE)
    }
  }
  dev.off()}
system(paste("pdfcrop --margin 10",plotname,plotname))

