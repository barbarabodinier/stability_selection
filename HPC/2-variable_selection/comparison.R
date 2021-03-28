library(focus)
library(glmnet)
library(stabs)

source("../additional_functions_specific_to_comparisons.R")

# Reading arguments
args=commandArgs(trailingOnly=TRUE)
simul_study_id=as.numeric(args[1])
do_exp=as.character(args[2])
sd_pred_error=as.numeric(args[3])
params_id=as.numeric(args[4])
seed=as.numeric(args[5])
PFER_thr=as.numeric(args[6])
simulation_id=paste0(params_id, "_", seed)
pi_list=seq(0.6, 0.9, by=0.05)

# Extracting simulation parameters
params_list=read.table(paste0("Simulation_parameters/Simulation_parameters_list_",simul_study_id,".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
n=params_list[params_id,"n"]
pk=params_list[params_id,"p"]
nu=params_list[params_id,"nu"]
p=sum(pk)

# Printing
print(paste("ID of simulation study:", simul_study_id))
print(paste("Number of observations:", n))
print(paste("Number of variables:", pk))
print(paste("Network density:", nu))
print(paste("Exponential:", do_exp))
print(paste("Simulation ID:", simulation_id))

# Creating folder of simulation study
filepath=paste0("../../Results/2-variable_selection/Simulations_",simul_study_id,"_linear/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Simulation
set.seed(seed)
simul=SimulateXY(n=n, pk=pk, nu_pred=nu, family="gaussian", sd_pred_error=sd_pred_error, continuous=TRUE)
X=simul$X
Y=simul$Y

# Lambda path
Lambda=LambdaGridRegression(xdata=simul$X, ydata=simul$Y)

# Univariate regressions
foo=function(){
Summary=matrix(NA, ncol=2, nrow=ncol(X))
colnames(Summary)=c("beta", "pval")
for (k in 1:ncol(X)){
  mymodel=lm(Y~X[,k])
  Summary[k,1]=summary(mymodel)$coefficients[2,1]
  Summary[k,2]=summary(mymodel)$coefficients[2,4]
}
Summary=data.frame(Summary)
assign("Summary", Summary, envir=.GlobalEnv)
}
tmptime=system.time(foo())

# Bonferroni correction
A=ifelse(p.adjust(Summary$pval, method="bonf")<0.05, yes=1, no=0)
perf_bonf=c(pi=NA, SelectionPerformance(theta=A, theta_star=simul$theta_pred), time=as.numeric(tmptime[1]))

# FDR correction
A=ifelse(p.adjust(Summary$pval, method="BH")<0.05, yes=1, no=0)
perf_bh=c(pi=NA, SelectionPerformance(theta=A, theta_star=simul$theta_pred), time=as.numeric(tmptime[1]))

# LASSO (no stability)
foo_lasso=function(){
mylasso=cv.glmnet(x=X, y=Y, lambda=Lambda)
assign("mylasso", mylasso, envir=.GlobalEnv)
}
tmptime=system.time(foo_lasso())

# Lambda min
A=ifelse(t(as.matrix(coef(mylasso, s="lambda.min")))!=0, yes=1, no=0)[-1]
perf_min=c(pi=NA, SelectionPerformance(theta=A, theta_star=simul$theta_pred), time=as.numeric(tmptime[1]))

# Lambda 1se
A=ifelse(t(as.matrix(coef(mylasso, s="lambda.1se")))!=0, yes=1, no=0)[-1]
perf_1se=c(pi=NA, SelectionPerformance(theta=A, theta_star=simul$theta_pred), time=as.numeric(tmptime[1]))

# Calibration by PFER (MB)
foo=function(){
  set.seed(1)
  stab_mb <- stabsel(x=X, y=Y, fitfun=glmnet.lasso_model, sampling.type="MB", mc.cores=1,
                     cutoff=pi_list[k], PFER=PFER_thr, args.fitfun=list(lambda=Lambda))
  assign("stab_mb", stab_mb, envir=.GlobalEnv)
}
perf_mb=NULL
for (k in 1:length(pi_list)){
  print(k)
  tmptime=system.time(foo())
  edgelist=cbind(gsub(" :.*", "", names(stab_mb$selected)), gsub(".*: ", "", names(stab_mb$selected)))
  A=rep(0, ncol(X))
  names(A)=colnames(X)
  A[stab_mb$selected]=1
  perf=data.frame(c(pi=pi_list[k], SelectionPerformance(theta=A, theta_star=simul$theta)))
  perf_mb=rbind(perf_mb, c(perf, time=as.numeric(tmptime[1])))
}

# Calibration by PFER (SS)
foo=function(){
  set.seed(1)
  stab_ss <- stabsel(x=X, y=Y, fitfun=glmnet.lasso_model, sampling.type="SS", assumption="unimodal", mc.cores=1,
                     cutoff=pi_list[k], PFER=PFER_thr, args.fitfun=list(lambda=Lambda))
  assign("stab_ss", stab_ss, envir=.GlobalEnv)
}
perf_ss=NULL
for (k in 1:length(pi_list)){
  print(k)
  tmptime=system.time(foo())
  edgelist=cbind(gsub(" :.*", "", names(stab_mb$selected)), gsub(".*: ", "", names(stab_mb$selected)))
  A=rep(0, ncol(X))
  names(A)=colnames(X)
  A[stab_ss$selected]=1
  perf=data.frame(c(pi=pi_list[k], SelectionPerformance(theta=A, theta_star=simul$theta)))
  perf_ss=rbind(perf_ss, c(perf, time=as.numeric(tmptime[1])))
}

# Unconstrained calibration (subsampling)
foo=function(){
out=VariableSelection(xdata=X, ydata=Y, family="gaussian", Lambda=Lambda, PFER_method="MB")
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo())
print(GetArgmaxId(out))
A=SelectedVariables(out)
perf=data.frame(c(pi=as.numeric(GetArgmax(out)[2]), SelectionPerformance(theta=A, theta_star=simul$theta)))
nperf=NULL
nperf=rbind(nperf, c(perf, time=as.numeric(tmptime[1])))

# Unconstrained calibration (CPSS)
foo=function(){
out=VariableSelection(xdata=X, ydata=Y, family="gaussian", Lambda=Lambda, PFER_method="SS")
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo())
print(GetArgmaxId(out))
A=SelectedVariables(out)
perf=data.frame(c(pi=as.numeric(GetArgmax(out)[2]), SelectionPerformance(theta=A, theta_star=simul$theta)))
nperf=rbind(nperf, c(perf, time=as.numeric(tmptime[1])))

# Constrained calibration (MB)
foo=function(){
out=VariableSelection(xdata=X, ydata=Y, family="gaussian", Lambda=Lambda, PFER_thr=PFER_thr)
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo())
print(GetArgmaxId(out))
A=SelectedVariables(out)
perf=data.frame(c(pi=as.numeric(GetArgmax(out)[2]), SelectionPerformance(theta=A, theta_star=simul$theta)))
nperf=rbind(nperf, c(perf, time=as.numeric(tmptime[1])))

# Constrained calibration (SS)
foo=function(){
out=VariableSelection(xdata=X, ydata=Y, family="gaussian", Lambda=Lambda, PFER_thr=PFER_thr, PFER_method="SS")
assign("out", out, envir=.GlobalEnv)
}
tmptime=system.time(foo())
print(GetArgmaxId(out))
A=SelectedVariables(out)
perf=data.frame(c(pi=as.numeric(GetArgmax(out)[2]), SelectionPerformance(theta=A, theta_star=simul$theta)))
nperf=rbind(nperf, c(perf, time=as.numeric(tmptime[1])))
rownames(nperf)=c("Subsampling", "CPSS", "MB", "SS")

# Re-formatting output object
perf_full=rbind(Bonferroni=perf_bonf, BH=perf_bh, lambda_min=perf_min, lambda_1se=perf_1se, perf_mb, perf_ss, nperf)

# Saving output object
saveRDS(perf_full, paste0(filepath, "Performances_nu", nu, "_", simulation_id, "_PFER_thr_", PFER_thr, ".rds"))

