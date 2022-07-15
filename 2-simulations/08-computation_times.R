rm(list = ls())
setwd("~/Dropbox/Stability_selection/")

library(focus)
library(plotrix)
library(colorspace)

# Simulation parameters
simul_study_id <- 1
simul_id=1
topology <- "random"

for (model in c("graphical", "regression")){
  if (model=="graphical"){
    performances <- readRDS(paste0("Results/2-simulations/1-graphical_model/Computation_times_", simul_study_id, "_", topology, "/Computation_times_", simul_id, "_merged.rds"))
  } else {
    performances <- readRDS(paste0("Results/2-simulations/3-regression_model/Computation_times_", simul_study_id, "/Computation_times_", simul_id, "_merged.rds"))
  }
  
  # Extracting the computation times
  mylist_time <- list()
  for (k in 1:nrow(performances)) {
    mylist_time <- c(mylist_time, list(as.numeric(performances[k, "time.user.self", ])))
    assign(paste0("median", k), median(as.numeric(performances[k, "time.user.self", ])))
  }
  
  mytable=matrix(paste0(formatC(sapply(mylist_time, median), format="f", digits=0, big.mark = ","), 
                        rep(" [", length(mylist_time)), 
                        formatC(sapply(mylist_time, IQR), format="f", digits=0, big.mark = ","), 
                        rep("]", length(mylist_time))), 
                 ncol=ifelse(model=="graphical", yes=2, no=1), byrow=TRUE)
  if (model=="graphical"){
    mytable=cbind(c(100, 250, 500, 750, 1000), mytable)
    colnames(mytable)=c("p", "Time (cold)", "Time (warm)")
  } else {
    mytable=cbind(c(1000, 2500, 5000, 7500, 10000), mytable)
    colnames(mytable)=c("p", "Time")
  }
  
  write.xlsx(as.data.frame(mytable), paste0("Tables/2-simulations/Table_computation_times_", model, "_", simul_study_id, "_", topology, ".xlsx"),
             rowNames = FALSE, colNames = TRUE
  )
  write.table(mytable, paste0("Tables/2-simulations/Table_computation_times_", model, "_", simul_study_id, "_", topology, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
  )
}

