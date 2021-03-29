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
filepath=paste0("../../Results/3-multi_block/Simulations_",simul_study_id,"_",topology,"/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Simulation
set.seed(seed)
simul=SimulateGraph(n=n, pk=pk, topology=topology, nu=nu, v_within=1, v_between=v_between)

# Lambda path
Lambda_single=LambdaGridNetwork(data=simul$data)

# Unconstrained calibration (subsampling)
foo_unconstr=function(){
out=GraphicalModel(data=simul$data, pk=pk, Lambda=Lambda_single, lambda_other_blocks=0.1, start="cold", PFER_method="MB")
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo_unconstr())
A=Adjacency(out)
nperf=SelectionPerformance(theta=A, theta_star=simul$theta, pk=pk)
rownames(nperf)=c("Overall",1:3)
perf_full=cbind(nperf, time=as.numeric(tmptime[1]))

# Saving output object
saveRDS(perf_full, paste0(filepath, "Performances_multi_nu", nu, "_", simulation_id, "_PFER_thr_", PFER_thr, ".rds"))

