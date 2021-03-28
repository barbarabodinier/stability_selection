library(abind)

args=commandArgs(trailingOnly=TRUE)
simul_study_id=as.numeric(args[1])
i=as.numeric(args[2])
PFER_thr=as.character(args[3])
data=""

print(paste("ID of the simulation study:", simul_study_id))
print(paste("Simulation ID:",i))
print(data)

setwd(paste0("../../Results/2-variable_selection/Simulations_",simul_study_id,"_linear/"))

myfiles=list.files(pattern=paste0("Performances",data))
myfiles=myfiles[grep(paste0("Performances_",data,"nu0.05_",i,"_.*_PFER_thr_", PFER_thr, ".rds"), myfiles)]
print(length(myfiles))
print(myfiles[1])
results=NULL
for (j in 1:length(myfiles)){
  tmp=readRDS(myfiles[j])
  results=abind(results, tmp, along=3)
}
saveRDS(results, paste0("Performances",data,"_",i,"_merged_PFER_thr_",PFER_thr,".rds"))
