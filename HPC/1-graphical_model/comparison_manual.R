library(fake)
library(sharp)
library(stabs)
library(pulsar)
library(glassoFast)

source("../additional_functions_specific_to_comparisons.R")


# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
topology <- as.character(args[2])
do_exp <- as.character(args[3])
params_id <- as.numeric(args[4])
seed <- as.numeric(args[5])
PFER_thr <- as.numeric(args[6])
simulation_id <- paste0(params_id, "_", seed)
pi_list <- seq(0.6, 0.9, by = 0.05)

# Extracting simulation parameters
params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
n <- params_list[params_id, "n"]
pk <- params_list[params_id, "p"]
nu <- params_list[params_id, "nu"]
p <- sum(pk)

# Printing
print(paste0("fake: ", packageVersion("fake")))
print(paste0("glmnet: ", packageVersion("glmnet")))
print(paste0("sharp: ", packageVersion("sharp")))
print(paste("ID of simulation study:", simul_study_id))
print(paste("Number of observations:", n))
print(paste("Number of variables:", pk))
print(paste("Network topology:", topology))
print(paste("Network density:", nu))
print(paste("Exponential:", do_exp))
print(paste("Simulation ID:", simulation_id))

# Creating folder of simulation study
filepath <- paste0("../../Results/1-graphical_model/Simulations_", simul_study_id, "_", topology, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Simulation
set.seed(seed)
simul <- SimulateGraphical(n = n, pk = pk, topology = topology, v_within = 1, nu_within = nu)

# Lambda path
Lambda_single <- LambdaGridGraphical(xdata = simul$data)

# Information theory
foo_it <- function() {
  out_it <- CalibrateInformationTheory(x = simul$data, Lambda = Lambda_single)
  assign("out_it", out_it, envir = .GlobalEnv)
}
tmptime <- system.time(foo_it())
print(which.min(out_it$AIC))
print(which.min(out_it$BIC))
print(which.min(out_it$EBIC))

# AIC
A <- out_it$path[, , which.min(out_it$AIC)]
perf <- data.frame(c(pi = NA, SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- c(perf, time = as.numeric(tmptime[1]))

# BIC
A <- out_it$path[, , which.min(out_it$BIC)]
perf <- data.frame(c(pi = NA, SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = as.numeric(tmptime[1])))

# EBIC
A <- out_it$path[, , which.min(out_it$EBIC)]
perf <- data.frame(c(pi = NA, SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = as.numeric(tmptime[1])))

# StARS
foo_stars <- function() {
  hugeargs <- list(lambda = Lambda_single)
  pulsar <- pulsar(data = simul$data, fun = glasso.pulsar, criterion = "stars", fargs = hugeargs, subsample.ratio = 0.5, rep.num = 100, ncores = 1)
  assign("pulsar", pulsar, envir = .GlobalEnv)
}
tmptime <- system.time(foo_stars())
opt_pulsar <- opt.index(pulsar, "stars")
A <- ifelse(as.matrix(pulsar$est$path[[opt_pulsar]]) != 0, yes = 1, no = 0)
perf <- data.frame(c(pi = NA, SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = as.numeric(tmptime[1])))

# Stability selection (subsampling)
time_unconstr_sub <- system.time({
  out_unconstr_sub <- GraphicalModel(xdata = simul$data, 
                                     Lambda = Lambda_single, 
                                     start = "cold", 
                                     PFER_method = "MB",
                                     verbose=FALSE)
})

# Stability selection (CPSS)
time_unconstr_cpss <- system.time({
  out_unconstr_cpss <- GraphicalModel(xdata = simul$data, 
                                      Lambda = Lambda_single, 
                                      start = "cold", 
                                      PFER_method = "SS",
                                      verbose=FALSE)
})

# Calibration by error control
perf_mb <- ErrorControl(out_unconstr_sub, simul = simul, time = as.numeric(time_unconstr_sub[1]))
perf_ss <- ErrorControl(out_unconstr_cpss, simul = simul, time = as.numeric(time_unconstr_cpss[1]))
nperf <- rbind(nperf, perf_mb, perf_ss)

# Unconstrained calibration
perf_sub <- data.frame(c(
  pi = as.numeric(Argmax(out_unconstr_sub)[2]),
  SelectionPerformance(theta = out_unconstr_sub, theta_star = simul$theta),
  time = as.numeric(time_unconstr_sub[1])
))
perf_cpss <- data.frame(c(
  pi = as.numeric(Argmax(out_unconstr_cpss)[2]),
  SelectionPerformance(theta = out_unconstr_cpss, theta_star = simul$theta),
  time = as.numeric(time_unconstr_cpss[1])
))

# Constrained stability selection (MB)
time_constr_sub <- system.time({
  out_constr_sub <- GraphicalModel(xdata = simul$data, 
                                   Lambda = Lambda_single, 
                                   start = "cold", 
                                   PFER_method = "MB",
                                   PFER_thr = PFER_thr,
                                   verbose=FALSE)
})

# Constrained stability selection (SS)
time_constr_cpss <- system.time({
  out_constr_cpss <- GraphicalModel(xdata = simul$data, 
                                    Lambda = Lambda_single, 
                                    start = "cold", 
                                    PFER_method = "SS",
                                    PFER_thr = PFER_thr,
                                    verbose=FALSE)
})

# Constrained calibration
perf_constr_mb <- data.frame(c(
  pi = as.numeric(Argmax(out_constr_sub)[2]),
  SelectionPerformance(theta = out_constr_sub, theta_star = simul$theta),
  time = as.numeric(time_constr_sub[1])
))
perf_constr_ss <- data.frame(c(
  pi = as.numeric(Argmax(out_constr_cpss)[2]),
  SelectionPerformance(theta = out_constr_cpss, theta_star = simul$theta),
  time = as.numeric(time_constr_cpss[1])
))
nperf <- rbind(nperf, perf_sub, perf_cpss, perf_constr_mb, perf_constr_ss)

# Re-formatting output object
rownames(nperf)[1:4] <- c("AIC", "BIC", "EBIC", "StARS")
rownames(nperf)[(nrow(nperf) - 2):nrow(nperf)] <- c("Unconstrained", "MB", "SS")
colnames(nperf)[ncol(nperf)] <- "time"
perf_full <- nperf

# Saving output object
if (topology == "scale-free") {
  saveRDS(perf_full, paste0(filepath, "Performances_", simulation_id, "_PFER_thr_", PFER_thr, ".rds"))
} else {
  saveRDS(perf_full, paste0(filepath, "Performances_nu", nu, "_", simulation_id, "_PFER_thr_", PFER_thr, ".rds"))
}
