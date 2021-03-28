args=commandArgs(trailingOnly=TRUE)
simul_study_id=as.numeric(args)

n=c(200,100,50)
p1=50
p2=50
nu=0.02
v_between=0.2

mylist=expand.grid(n=n, p1=p1, p2=p2, nu=nu, v_between=v_between)
mylist=cbind(simulation_id=1:nrow(mylist), mylist)

dir.create("Simulation_parameters", showWarnings=FALSE)
write.table(mylist, paste0("Simulation_parameters/Simulation_parameters_list_",simul_study_id,".txt"), sep="\t", row.names=FALSE)
