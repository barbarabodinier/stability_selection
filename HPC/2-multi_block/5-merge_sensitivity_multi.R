library(abind)

args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
topology <- as.character(args[2])
i <- as.numeric(args[3])
data <- ""

print(paste("ID of the simulation study:", simul_study_id))
print(paste("Topology:", topology))
print(paste("Simulation ID:", i))
print(data)

setwd(paste0("../../Results/2-multi_block/Sensitivity_", simul_study_id, "_", topology))

myfiles <- list.files(pattern = paste0("Performances", data))
myfiles <- myfiles[grep(paste0("Performances_multi_sensitivity_nu.*_", i, "_.*_PFER_thr_Inf.rds"), myfiles)]
print(length(myfiles))
results <- NULL
for (j in 1:length(myfiles)) {
  tmp <- readRDS(myfiles[j])
  results <- abind(results, tmp, along = 3)
}
saveRDS(results, paste0("Performances_multi_sensitivity_", i, "_merged_PFER_thr_Inf.rds"))
