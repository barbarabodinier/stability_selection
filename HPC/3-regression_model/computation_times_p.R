library(sharp)

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
seed <- as.numeric(args[2])
params_id <- 1
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

# Defining other simulation parameters
n <- 500
p_list <- c(1000, 2500, 5000, 7500, 10000)

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
filepath <- paste0("../../Results/3-regression_model/Computation_times_", simul_study_id, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Running for different numbers of variables p
nperf <- NULL
for (p in p_list) {
  print(p)

  # Simulation
  set.seed(seed)
  simul <- SimulateRegression(
    n = n, pk = p, nu_xz = nu_xz,
    eta_set = 1, nu_zy = 1,
    v_within = 1, nu_within = nu_within,
    family = "gaussian", ev_xz = ev
  )
  simul$ydata <- simul$ydata[, 1]
  simul$theta <- simul$theta[, 1, drop = FALSE]
  print(table(simul$theta))

  # Lambda path
  Lambda_single <- LambdaGridRegression(xdata = simul$xdata, ydata = simul$ydata)

  # Unconstrained calibration
  foo_unconstr <- function() {
    out <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, Lambda = Lambda_single, PFER_method = "MB", verbose = FALSE)
    assign("out", out, envir = .GlobalEnv)
  }

  # Selection performance
  tmptime <- system.time(foo_unconstr())
  perf <- data.frame(c(pi = as.numeric(Argmax(out)[2]), SelectionPerformance(theta = out, theta_star = simul$theta)))
  nperf <- rbind(nperf, c(perf, time = tmptime[1]))

  # Saving intermediate output
  saveRDS(nperf, paste0(filepath, "Computation_times_up_to_", p, "_", simulation_id, ".rds"))
}

# Saving output
saveRDS(nperf, paste0(filepath, "Computation_times_", simulation_id, ".rds"))
