library(focus)

# Reading arguments
args=commandArgs(trailingOnly=TRUE)
simul_study_id=as.numeric(args[1])
topology=as.character(args[2])
do_exp=as.character(args[3])
params_id=as.numeric(args[4])
seed=as.numeric(args[5])
PFER_thr=as.numeric(args[6])
simulation_id=paste0(params_id, "_", seed)

# Extracting simulation parameters
params_list=read.table(paste0("Simulation_parameters/Simulation_parameters_list_",simul_study_id,".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
n=params_list[params_id,"n"]
p1=params_list[params_id,"p1"]
p2=params_list[params_id,"p2"]
pk=c(p1,p2)
nu=params_list[params_id,"nu"]
v_between=params_list[params_id,"v_between"]

# Printing
print(paste("ID of simulation study:", simul_study_id))
print(paste("Number of observations:", n))
print(paste("Number of variables:", pk))
print(paste("Network topology:", topology))
print(paste("Network density:", nu))
print(paste("v_between:", v_between))
print(paste("Exponential:", do_exp))
print(paste("Simulation ID:", simulation_id))

# Creating folder of simulation study
filepath=paste0("../../Results/3-multi_block/Sensitivity_",simul_study_id,"_",topology,"/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Simulation
set.seed(seed)
simul=SimulateGraph(n=n, pk=pk, topology=topology, nu=nu, v_within=1, v_between=v_between)

# Lambda path
Lambda_single=LambdaGridNetwork(data=simul$data, Lambda_cardinal=30)
Lambda_blocks=expand.grid(Lambda_single,Lambda_single)
Lambda_blocks=as.matrix(cbind(Lambda_blocks,Lambda_blocks[,1]))

# Multi-parameters
foo_unconstr=function(){
out=GraphicalModel(data=simul$data, pk=pk, Lambda=Lambda_blocks, start="cold", PFER_method="MB")
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo_unconstr())
print(tmptime)
print(GetArgmax(out))
A=Adjacency(out)
nperf=SelectionPerformance(theta=A, theta_star=simul$theta, pk=pk)
rownames(nperf)=c("Overall",1:3)
perf_full=cbind(nperf, time=as.numeric(tmptime[1]))

# Multi-block (different lambda_dense)
myperfs=perf_full
for (lambda_dense in c(0,0.001,0.01,0.1,0.5,1)){
print(lambda_dense)
foo_unconstr=function(){
out=GraphicalModel(data=simul$data, pk=pk, Lambda=Lambda_single, lambda_other_blocks=lambda_dense, start="cold", PFER_method="MB")
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo_unconstr())
print(tmptime)
print(GetArgmax(out))
A=Adjacency(out)
nperf=SelectionPerformance(theta=A, theta_star=simul$theta, pk=pk)
rownames(nperf)=c("Overall",1:3)
perf_full=cbind(nperf, time=as.numeric(tmptime[1]))
myperfs=rbind(myperfs, perf_full)
}

# Saving output object
saveRDS(myperfs, paste0(filepath, "Performances_multi_sensitivity_nu", nu, "_", simulation_id, "_PFER_thr_", PFER_thr, ".rds"))

