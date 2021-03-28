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

setwd(paste0("../../Results/1-graphical_model/Simulations_",simul_study_id,"_",topology))

myfiles=list.files(pattern=paste0("Performances",data))
if (topology=="random"){
  myfiles=myfiles[grep(paste0("Performances_",data,"nu0.02_",i,"_.*_PFER_thr_", PFER_thr, ".rds"), myfiles)]
} else {
  myfiles=myfiles[grep(paste0("Performances_",data,i,"_.*_PFER_thr_", PFER_thr, ".rds"), myfiles)]
}
print(length(myfiles))
results=NULL
for (j in 1:min(1000, length(myfiles))){
  #print(myfiles[j])
  tmp=readRDS(myfiles[j])
  results=abind(results, tmp, along=3)
}
saveRDS(results, paste0("Performances",data,"_",i,"_merged_PFER_thr_",PFER_thr,".rds"))
