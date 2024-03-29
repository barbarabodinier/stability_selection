library(abind)

args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
PFER_thr <- as.character(args[2])
i <- 3
data <- ""
if (simul_study_id == 1) {
  nu_within <- 0
} else {
  nu_within <- 0.02
}

print(paste("ID of the simulation study:", simul_study_id))
print(paste("Simulation ID:", i))
print(data)

setwd(paste0("../../Results/3-regression_model/Simulations_", simul_study_id))

myfiles <- list.files(pattern = paste0("Performances", data))
myfiles <- myfiles[grep(paste0("Performances_nu.*_", nu_within, "_", i, "_.*_PFER_thr_", PFER_thr, ".rds"), myfiles)]
print(length(myfiles))
results <- NULL
for (j in 1:min(1000, length(myfiles))) {
  tmp <- readRDS(myfiles[j])
  results <- abind(results, tmp, along = 3)
}
saveRDS(results, paste0("Performances_", i, "_merged_PFER_thr_", PFER_thr, ".rds"))
