library(focus)

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
topology <- as.character(args[2])
do_exp <- as.character(args[3])
params_id <- as.numeric(args[4])
seed <- as.numeric(args[5])
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
filepath <- paste0("../../Results/1-graphical_model/Sensitivity_K_", simul_study_id, "_", topology, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Simulation
set.seed(seed)
simul <- SimulateGraphical(n = n, pk = pk, topology = topology, nu = nu)
C_hat <- cor(simul$data)

# Lambda path
Lambda_single <- LambdaGridGraphical(data = simul$data, lambda_path_refined_cardinal = 100)

# Unconstrained calibration
foo_unconstr <- function() {
  out <- GraphicalModel(xdata = simul$data, Lambda = Lambda_single, start = "cold", K = k)
  assign("out", out, envir = .GlobalEnv)
}

# Running for different numbers of iterations K
nperf <- NULL
for (k in c(10, 20, 50, 100, 500, 1000, 2000, 5000, 10000)) {
  print(k)
  tmptime <- system.time(foo_unconstr())
  A <- Adjacency(out)
  perf <- SelectionPerformance(theta = A, theta_star = simul$theta)
  nperf <- rbind(nperf, c(perf, time = tmptime[1]))
}

# Saving output
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))
