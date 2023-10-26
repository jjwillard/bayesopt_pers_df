## Example simulation used in "Bayesian optimization for personalized dose-finding trials with 
## combination therapies"
## 2023-10-26

pacman::p_load(laGP, ggdist, dplyr, tidyr, purrr, foreach, doParallel)
source("./simulation_functions.R")

#------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------#

n_simulations <- 1100
standard <- tibble(filename = "./Data/bimodal_resp_het_regular.RDS",
                   iterations = 15,
                   reps = 4,
                   delta = 0,
                   min_stop_iteration = 15,
                   std_dev = 0.2)
personalized <- tibble(filename = "./Data/bimodal_resp_het_personalized.RDS",
                   iterations = 15,
                   reps = 2,
                   delta = 0,
                   min_stop_iteration = 15,
                   std_dev = 0.2)

objective_function <- function(x_1, x_2, x_3){
  data_x <- as.matrix(tibble(x_1 = x_1, x_2 = x_2))
  
  # same for scenario 3
  if(x_3 == 0){
    f_x <-  -mvtnorm::dmvnorm(data_x,
                              mean = c(0.3, 0.7),
                              sigma = matrix(c(0.2, 0.05, 0.05, 0.1), ncol = 2))
  } else {
    f_x <-  -mvtnorm::dmvnorm(data_x,
                              mean = c(0.7, 0.3),
                              sigma = matrix(c(0.2, 0.05, 0.05, 0.1), ncol = 2))
  }
  tibble(f_x = f_x)
  
}

optimal_value_x_3_0 <- pull(objective_function(0.3, 0.7, 0))
optimal_value_x_3_1 <- pull(objective_function(0.7, 0.3, 1))
optimal_dose_combo_x_3_0 <- c(0.3, 0.7)
optimal_dose_combo_x_3_1 <- c(0.7, 0.3)

#------------------------------------------------------------------------#
# Standard Algorithm
#------------------------------------------------------------------------#

#tictoc::tic()
registerDoParallel(cores = parallel::detectCores())
seed <-  100000:(100000 + n_simulations - 1)
sim_res_reg <- foreach(j=1:n_simulations,
                       .errorhandling = "remove",
                       .combine = 'rbind') %dopar% run_sim_reg(j, 
                                                               iterations = standard$iterations, 
                                                               reps = standard$reps, 
                                                               delta = standard$delta,
                                                               min_stop_iteration = standard$min_stop_iteration, #allows you stop AFTER the specified iteration
                                                               optimal_dose_combo_x_3_0 = optimal_dose_combo_x_3_0,
                                                               optimal_dose_combo_x_3_1 = optimal_dose_combo_x_3_1,
                                                               optimal_value_x_3_0 = optimal_value_x_3_0,
                                                               optimal_value_x_3_1 = optimal_value_x_3_1,
                                                               std_dev = standard$std_dev,
                                                               seed = seed[j])
stopImplicitCluster()

saveRDS(sim_res_reg, standard$filename)

#------------------------------------------------------------------------#
# Personalized Algorithm
#------------------------------------------------------------------------#

#tictoc::tic()
registerDoParallel(cores = parallel::detectCores())
sim_res_pers <- foreach(j=1:n_simulations,
                        .errorhandling = "remove",
                        .combine = 'rbind') %dopar% {
                          run_sim_pers(j, 
                                       iterations = personalized$iterations, 
                                       reps = personalized$reps, 
                                       delta = personalized$delta,
                                       min_stop_iteration = personalized$min_stop_iteration,
                                       optimal_dose_combo_x_3_0 = optimal_dose_combo_x_3_0,
                                       optimal_dose_combo_x_3_1 = optimal_dose_combo_x_3_1,
                                       optimal_value_x_3_0 = optimal_value_x_3_0,
                                       optimal_value_x_3_1 = optimal_value_x_3_1,
                                       std_dev = personalized$std_dev,
                                       seed = seed[j])
                        }
  
stopImplicitCluster()

saveRDS(sim_res_pers, personalized$filename)

