args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args)

n <- c(2000, 1000, 500)
p <- 1000
nu_within <- 0.02
nu_xz <- 0.02
ev <- 0.4

mylist <- expand.grid(n = n, p = p, nu_within = nu_within, nu_xz = nu_xz, ev = ev)
mylist <- cbind(simulation_id = 1:nrow(mylist), mylist)

dir.create("Simulation_parameters", showWarnings = FALSE)
write.table(mylist, paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"), sep = "\t", row.names = FALSE)
