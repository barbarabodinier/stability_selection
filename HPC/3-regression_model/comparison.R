library(sharp)
library(stabs)
library(glmnet)

source("../additional_functions_specific_to_comparisons.R")


# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
params_id <- as.numeric(args[2])
seed <- as.numeric(args[3])
PFER_thr <- as.numeric(args[4])
simulation_id <- paste0(params_id, "_", seed)
pi_list <- seq(0.6, 0.9, by = 0.05)

# Extracting simulation parameters
params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
n <- params_list[params_id, "n"]
nu_within <- params_list[params_id, "nu_within"]
nu_xz <- params_list[params_id, "nu_xz"]
ev <- params_list[params_id, "ev"]
pk <- rep(100, 10)
p <- sum(pk)

# Printing
print(packageVersion("sharp"))
print(paste("ID of simulation study:", simul_study_id))
print(paste("Number of observations:", n))
print(paste("Number of variables:", pk))
print(paste("Network density:", nu_within))
print(paste("Proportion of contributing predictors:", nu_xz))
print(paste("Proportion of explained variance:", ev))
print(paste("Simulation ID:", simulation_id))

# Creating folder of simulation study
filepath <- paste0("../../Results/3-regression_model/Simulations_", simul_study_id, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Simulation
set.seed(seed)
adjacency <- SimulateAdjacency(
  pk = pk,
  nu_within = nu_within,
  nu_between = 0
)
theta <- sharp:::SamplePredictors(pk = pk, nu = nu_xz)
theta <- apply(theta, 1, sum)
simul <- SimulateRegression(
  n = n, pk = sum(pk),
  adjacency_x = adjacency,
  theta_xz = cbind(theta),
  eta_set = 1,
  v_within = 1,
  family = "gaussian",
  ev_xz = ev
)
print(table(simul$theta))

# Lambda path
Lambda_single <- LambdaGridRegression(xdata = simul$xdata, ydata = simul$ydata)

# Cross-validation
tmptime <- system.time({
  out <- cv.glmnet(x = simul$xdata, y = simul$ydata, lambda = Lambda_single, family = "gaussian")
})
out_min <- ifelse(as.vector(coef.glmnet(out, s = "lambda.min"))[-1] != 0, yes = 1, no = 0)
out_1se <- ifelse(as.vector(coef.glmnet(out, s = "lambda.1se"))[-1] != 0, yes = 1, no = 0)
perf_min <- c(data.frame(c(pi = NA, SelectionPerformance(theta = out_min, theta_star = simul$theta))),
  time = as.numeric(tmptime[1])
)
perf_1se <- c(data.frame(c(pi = NA, SelectionPerformance(theta = out_1se, theta_star = simul$theta))),
  time = as.numeric(tmptime[1])
)
nperf <- rbind(
  perf_min,
  perf_1se
)

# Calibration by PFER (MB)
foo_mb <- function() {
  set.seed(1)
  stab_mb <- stabsel(
    x = simul$xdata, y = simul$ydata, fitfun = glmnet.lasso_model, sampling.type = "MB", mc.cores = 1,
    cutoff = pi_list[k], PFER = PFER_thr, args.fitfun = list(lams = Lambda_single)
  )
  assign("stab_mb", stab_mb, envir = .GlobalEnv)
}
perf_mb <- NULL
for (k in 1:length(pi_list)) {
  print(k)
  tmptime <- system.time(foo_mb())
  selected <- rep(0, ncol(simul$xdata))
  names(selected) <- colnames(simul$xdata)
  selected[names(stab_mb$selected)] <- 1
  perf <- data.frame(c(pi = pi_list[k], SelectionPerformance(theta = selected, theta_star = simul$theta)))
  perf_mb <- rbind(perf_mb, c(perf, time = tmptime[1]))
}
nperf <- rbind(nperf, perf_mb)

# Calibration by PFER (SS)
foo_mb <- function() {
  set.seed(1)
  stab_mb <- stabsel(
    x = simul$xdata, y = simul$ydata, fitfun = glmnet.lasso_model, sampling.type = "SS", assumption = "unimodal", mc.cores = 1,
    cutoff = pi_list[k], PFER = PFER_thr, args.fitfun = list(lams = Lambda_single)
  )
  assign("stab_mb", stab_mb, envir = .GlobalEnv)
}
perf_mb <- NULL
for (k in 1:length(pi_list)) {
  print(k)
  tmptime <- system.time(foo_mb())
  selected <- rep(0, ncol(simul$xdata))
  names(selected) <- colnames(simul$xdata)
  selected[names(stab_mb$selected)] <- 1
  perf <- data.frame(c(pi = pi_list[k], SelectionPerformance(theta = selected, theta_star = simul$theta)))
  perf_mb <- rbind(perf_mb, c(perf, time = tmptime[1]))
}
nperf <- rbind(nperf, perf_mb)

# Unconstrained calibration (subsampling)
foo_unconstr <- function() {
  out <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, Lambda = Lambda_single, PFER_method = "MB", verbose = FALSE)
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = out, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Unconstrained calibration (simultaneous selection in complementary pairs)
foo_unconstr <- function() {
  out <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, Lambda = Lambda_single, PFER_method = "SS", verbose = FALSE)
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = out, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Constrained calibration (MB)
foo_unconstr <- function() {
  out <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata, Lambda = Lambda_single,
    PFER_method = "MB", PFER_thr = PFER_thr, verbose = FALSE
  )
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = out, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Constrained calibration (SS)
foo_unconstr <- function() {
  out <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata, Lambda = Lambda_single,
    PFER_method = "SS", PFER_thr = PFER_thr, verbose = FALSE
  )
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = out, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Saving output object
saveRDS(nperf, paste0(filepath, "Performances_nu", nu_xz, "_", nu_within, "_", simulation_id, "_PFER_thr_", PFER_thr, ".rds"))
