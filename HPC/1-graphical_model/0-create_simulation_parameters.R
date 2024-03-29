args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args)

n <- c(200, 100, 50)
p <- 100
nu <- 0.02

mylist <- expand.grid(n = n, p = p, nu = nu)
mylist <- cbind(simulation_id = 1:nrow(mylist), mylist)

dir.create("Simulation_parameters", showWarnings = FALSE)
write.table(mylist, paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"), sep = "\t", row.names = FALSE)
