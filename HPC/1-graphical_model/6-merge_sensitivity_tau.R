library(abind)

args=commandArgs(trailingOnly=TRUE)
simul_study_id=as.numeric(args[1])
topology=as.character(args[2])
i=as.numeric(args[3])
PFER_thr=as.character(args[4])
data=""

print(paste("ID of the simulation study:", simul_study_id))
print(paste("Topology:",topology))
print(paste("Simulation ID:",i))
print(data)

setwd(paste0("../../Results/1-graphical_model/Sensitivity_tau_",simul_study_id,"_",topology))

myfiles=list.files(pattern=paste0("Performances",data))
myfiles=myfiles[grep(paste0("Performances_",i,"_.*.rds"), myfiles)]
print(length(myfiles))

results=NULL
for (j in 1:length(myfiles)){
  tmp=readRDS(myfiles[j])
  results=abind(results, tmp, along=3)
}

saveRDS(results, paste0("Performances_",i,"_merged.rds"))
