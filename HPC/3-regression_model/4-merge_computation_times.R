library(abind)

simul_study_id <- 1
i <- 1

print(paste("ID of the simulation study:", simul_study_id))
print(paste("Simulation ID:", i))

setwd(paste0("../../Results/3-regression_model/Computation_times_", simul_study_id))

myfiles <- list.files(pattern = paste0("Computation_times_", i, "_"))
myfiles <- myfiles[1:1000]
print(length(myfiles))

results <- NULL
for (j in 1:length(myfiles)) {
  tmp <- readRDS(myfiles[j])
  results <- abind(results, tmp, along = 3)
}

saveRDS(results, paste0("Computation_times_", i, "_merged.rds"))
