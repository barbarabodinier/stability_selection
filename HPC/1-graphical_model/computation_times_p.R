library(fake)
library(sharp)

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
topology <- as.character(args[2])
seed <- as.numeric(args[3])
params_id <- 1
simulation_id <- paste0(params_id, "_", seed)
pi_list <- seq(0.6, 0.9, by = 0.05)

# Extracting simulation parameters
params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
nu <- params_list[params_id, "nu"]

# Defining other simulation parameters
n <- 500
p_list <- c(100, 250, 500, 750, 1000)

# Printing
print(paste0("fake: ", packageVersion("fake")))
print(paste0("glmnet: ", packageVersion("glmnet")))
print(paste0("sharp: ", packageVersion("sharp")))
print(paste("ID of simulation study:", simul_study_id))
print(paste("Number of observations:", n))
print(paste("Numbers of variables:"))
print(p_list)
print(paste("Network topology:", topology))
print(paste("Network density:", nu))
print(paste("Simulation ID:", simulation_id))

# Creating folder of simulation study
filepath <- paste0("../../Results/1-graphical_model/Computation_times_", simul_study_id, "_", topology, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Running for different numbers of iterations K
nperf <- NULL
for (p in p_list) {
  print(p)

  # Simulation
  set.seed(seed)
  simul <- SimulateGraphical(n = n, pk = p, topology = topology, v_within = 1, nu_within = nu)

  # Lambda path
  Lambda_single <- LambdaGridGraphical(xdata = simul$data, lambda_path_refined_cardinal = 100)

  for (start in c("cold", "warm")) {
    # Unconstrained calibration
    foo_unconstr <- function() {
      out <- GraphicalModel(xdata = simul$data, start = start, Lambda = Lambda_single)
      assign("out", out, envir = .GlobalEnv)
    }

    # Selection performance
    tmptime <- system.time(foo_unconstr())
    A <- Adjacency(out)
    perf <- SelectionPerformance(theta = A, theta_star = simul$theta)
    nperf <- rbind(nperf, c(perf, time = tmptime[1]))
  }

  # Saving intermediate output
  saveRDS(nperf, paste0(filepath, "Computation_times_up_to_", p, "_", simulation_id, ".rds"))
}

# Saving output
saveRDS(nperf, paste0(filepath, "Computation_times_", simulation_id, ".rds"))
