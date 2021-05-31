library(focus)
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
simul <- SimulateGraphical(n = n, pk = pk, topology = topology, nu = nu)
C_hat <- cor(simul$data)

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

# Calibration by PFER (MB)
foo_mb <- function() {
  set.seed(1)
  stab_mb <- stabsel(
    x = simul$data, fitfun = glasso.graphical_model, sampling.type = "MB", mc.cores = 1,
    cutoff = pi_list[k], PFER = PFER_thr, args.fitfun = list(lams = Lambda_single)
  )
  assign("stab_mb", stab_mb, envir = .GlobalEnv)
}
perf_mb <- NULL
for (k in 1:length(pi_list)) {
  print(k)
  tmptime <- system.time(foo_mb())
  edgelist <- cbind(gsub(" :.*", "", names(stab_mb$selected)), gsub(".*: ", "", names(stab_mb$selected)))
  p <- ncol(simul$data)
  A <- matrix(0, p, p)
  colnames(A) <- rownames(A) <- colnames(simul$data)
  A[edgelist] <- 1
  A <- A + t(A)
  perf <- data.frame(c(pi = pi_list[k], SelectionPerformance(theta = A, theta_star = simul$theta)))
  perf_mb <- rbind(perf_mb, c(perf, time = tmptime[1]))
}
nperf <- rbind(nperf, perf_mb)

# Calibration by PFER (SS)
foo_mb <- function() {
  set.seed(1)
  stab_mb <- stabsel(
    x = simul$data, fitfun = glasso.graphical_model, sampling.type = "SS", assumption = "unimodal", mc.cores = 1,
    cutoff = pi_list[k], PFER = PFER_thr, args.fitfun = list(lams = Lambda_single)
  )
  assign("stab_mb", stab_mb, envir = .GlobalEnv)
}
perf_mb <- NULL
for (k in 1:length(pi_list)) {
  print(k)
  tmptime <- system.time(foo_mb())
  edgelist <- cbind(gsub(" :.*", "", names(stab_mb$selected)), gsub(".*: ", "", names(stab_mb$selected)))
  p <- ncol(simul$data)
  A <- matrix(0, p, p)
  colnames(A) <- rownames(A) <- colnames(simul$data)
  A[edgelist] <- 1
  A <- A + t(A)
  perf <- data.frame(c(pi = pi_list[k], SelectionPerformance(theta = A, theta_star = simul$theta)))
  perf_mb <- rbind(perf_mb, c(perf, time = tmptime[1]))
}
nperf <- rbind(nperf, perf_mb)

# Unconstrained calibration (subsampling)
foo_unconstr <- function() {
  out <- GraphicalModel(xdata = simul$data, Lambda = Lambda_single, start = "cold", PFER_method = "MB")
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
A <- Adjacency(out)
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Unconstrained calibration (simultaneous selection in complementary pairs)
foo_unconstr <- function() {
  out <- GraphicalModel(xdata = simul$data, Lambda = Lambda_single, start = "cold", PFER_method = "SS")
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
A <- Adjacency(out)
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Constrained calibration (MB)
foo_unconstr <- function() {
  out <- GraphicalModel(xdata = simul$data, Lambda = Lambda_single, start = "cold", PFER_method = "MB", PFER_thr = PFER_thr)
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
A <- Adjacency(out)
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

# Constrained calibration (SS)
foo_unconstr <- function() {
  out <- GraphicalModel(xdata = simul$data, Lambda = Lambda_single, start = "cold", PFER_method = "SS", PFER_thr = PFER_thr)
  assign("out", out, envir = .GlobalEnv)
}
tmptime <- system.time(foo_unconstr())
A <- Adjacency(out)
perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = A, theta_star = simul$theta)))
nperf <- rbind(nperf, c(perf, time = tmptime[1]))

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
