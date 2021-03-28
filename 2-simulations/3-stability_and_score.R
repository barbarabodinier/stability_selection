rm(list=ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(glassoFast)
library(colorspace)

source("Scripts/additional_functions_specific_to_comparisons.R")


### Simulations

# Simulation parameters
n=200
pk=100
topology="random"
nu=0.02
seed=0

# Stability selection parameters
PFER_thr=10
K=100

# Data simulation
set.seed(seed)
simul=SimulateGraph(n=n, pk=pk, topology=topology, nu=nu, output_matrices=TRUE)


### Graphical models

# Definition of the grid of penalty parameters
Lambda=LambdaGridNetwork(data=simul$data, max_density=0.9)
Lambda_sub=LambdaGridNetwork(data=simul$data, max_density=0.5) # to subset for stability selection

# Graphical LASSO (no stability)
out_it=CalibrateInformationTheory(x=simul$data, Lambda=Lambda)
perf_original=NULL
for (k in 1:length(Lambda)){
  print(k)
  A=out_it$path[,,k]
  perf_original=rbind(perf_original, data.frame(SelectionPerformance(theta=A, theta_star=simul$theta)))
}
perf_original=perf_original[perf_original$precision+perf_original$recall>0,]
perf_bic=data.frame(SelectionPerformance(theta=out_it$path[,,which.min(out_it$BIC)], theta_star=simul$theta))
perf_ebic=data.frame(SelectionPerformance(theta=out_it$path[,,which.min(out_it$EBIC)], theta_star=simul$theta))

# Stability selection
system.time({out_visited=GraphicalModel(data=simul$data, Lambda=Lambda[Lambda>=min(Lambda_sub)], K=K, n_cores=2)})
perf_visited=NULL
for (i in 1:nrow(out_visited$S_2d)){
  for (j in 1:ncol(out_visited$S_2d)){
    A=Adjacency(out_visited, argmax_id=rbind(c(i,j)))
    perf_visited=rbind(perf_visited, 
                       data.frame(c(SelectionPerformance(theta=A, theta_star=simul$theta), 
                                    lambda=as.numeric(out_visited$Lambda[i,1]), 
                                    pi=out_visited$params$pi_list[j], 
                                    stability=out_visited$S_2d[i,j])))
  }
}
perf_visited=perf_visited[perf_visited$precision+perf_visited$recall>0,]
myperf=data.frame(SelectionPerformance(theta=Adjacency(out_visited), theta_star=simul$theta))

# Constrained calibration
system.time({out_constr=GraphicalModel(data=simul$data, Lambda=out_visited$Lambda[1:max(which(apply(out_visited$PFER_2d,1,min)<=PFER_thr)),], 
                                       K=K, PFER_thr=PFER_thr)})
myperf_constr=data.frame(SelectionPerformance(theta=Adjacency(out_constr), theta_star=simul$theta))


### Visualisation

# Figure parameters
mycolours=c("slategrey", "navy",
            "chocolate", "sienna", "darkred", "red")

# Saving figure
{pdf("Figures/2-simulations/Stability_selection_stability_score.pdf", width=14, height=7)
  par(mar=c(5,5,1,1), mfrow=c(1,2))
  plot(perf_visited$precision, perf_visited$recall, 
       las=1, cex.main=2, xlim=c(0,1), ylim=c(0,1),
       xlab="Precision", ylab="Recall", 
       cex.lab=1.5, col=mycolours[1],
       pch=17, cex=0.5, las=1, lty=1, bty="n")
  abline(h=axTicks(2), lty=2, col="lightgrey", lwd=0.5)
  abline(v=axTicks(1), lty=2, col="lightgrey", lwd=0.5)
  abline(h=perf_bic$recall, lty=3, col=mycolours[3])
  abline(v=perf_bic$precision, lty=3, col=mycolours[3])
  abline(h=perf_ebic$recall, lty=3, col=mycolours[4])
  abline(v=perf_ebic$precision, lty=3, col=mycolours[4])
  id=which.max(perf_visited$stability)
  abline(h=perf_visited$recall[id], lty=3, col=mycolours[5])
  abline(v=perf_visited$precision[id], lty=3, col=mycolours[5])
  abline(h=myperf_constr$recall, col=mycolours[6], lty=3)
  abline(v=myperf_constr$precision, col=mycolours[6], lty=3)
  points(perf_visited$precision, perf_visited$recall, 
         col=mycolours[1], pch=17, cex=0.5)
  points(perf_original$precision, perf_original$recall, 
         col=mycolours[2], pch=16, cex=0.7)
  lines(perf_original$precision, perf_original$recall, 
        col=mycolours[2],lty=1)
  points(perf_bic$precision, perf_bic$recall, pch=18, cex=1.2, col=mycolours[3])
  points(perf_ebic$precision, perf_ebic$recall, pch=18, cex=1.2, col=mycolours[4])
  points(perf_visited$precision[id], perf_visited$recall[id], pch=18, cex=1.2, col=mycolours[5])
  points(myperf_constr$precision, myperf_constr$recall, pch=18, cex=1.2, col=mycolours[6])
  legend("bottomleft", pch=c(16,rep(18,2),17,rep(18,2)), col=mycolours[c(2,3,4,1,5,6)], 
         pt.cex=c(0.7,rep(1.2,2),0.5,rep(1.2,2)), bg="white", lty=c(1, 3, 3, NA, 3, 3), bty="n",
         legend=c("Graphical LASSO", "BIC", expression(EBIC~(gamma*"="*0.5)),
                  "Stability selection", "Unconstrained", 
                  eval(parse(text="expression('Constrained ('*PFER[MB]*'<10)')"))))
  mtext("A", side=2, line=3, at=1, cex=3, las=1)
  
  plot(perf_visited$stability, perf_visited$F1_score, las=1, 
       col=mycolours[1], pch=17, cex=0.5, cex.lab=1.5, ylim=c(0,1),
       xlab="Stability score", ylab=expression(F[1]-score), bty="n")
  abline(h=axTicks(2), lty=2, col="lightgrey", lwd=0.5)
  abline(h=perf_bic$F1_score, lty=3, col=mycolours[3])
  abline(h=perf_ebic$F1_score, lty=3, col=mycolours[4])
  abline(v=axTicks(1), lty=2, col="lightgrey", lwd=0.5)
  abline(h=perf_visited$F1_score[id], col=mycolours[5], lty=3)
  abline(v=perf_visited$stability[id], col=mycolours[5], lty=3)
  abline(h=myperf_constr$F1_score, col=mycolours[6], lty=3)
  abline(v=max(out_constr$S, na.rm=TRUE), col=mycolours[6], lty=3)
  points(perf_visited$stability[id], perf_visited$F1_score[id], pch=18, cex=1.2, col=mycolours[5])
  points(max(out_constr$S, na.rm=TRUE), myperf_constr$F1_score, pch=18, cex=1.2, col=mycolours[6])
  mtext("B", side=2, line=3, at=1, cex=3, las=1)
  dev.off()}

