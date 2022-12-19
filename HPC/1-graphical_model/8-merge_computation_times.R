library(abind)

args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- 1
topology <- "random"
i <- 1

print(paste("ID of the simulation study:", simul_study_id))
print(paste("Topology:", topology))
print(paste("Simulation ID:", i))

setwd(paste0("../../Results/1-graphical_model/Computation_times_", simul_study_id, "_", topology))

myfiles <- list.files(pattern = paste0("Computation_times_", i, "_"))
myfiles <- myfiles[1:min(length(myfiles), 1000)]
print(length(myfiles))

results <- NULL
for (j in 1:length(myfiles)) {
  tmp <- readRDS(myfiles[j])
  results <- abind(results, tmp, along = 3)
}

saveRDS(results, paste0("Computation_times_", i, "_merged.rds"))
