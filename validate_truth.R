source(here("Survival-simulation-master/code_simulation/simulate_data.R"))

test <- simulate_data(n_sim = 1e6)
k_grid <- 1:10

mc_surv_est <- sapply(k_grid-1,function(k)mean(k<test$dat$T1))
wrong_truth <- test$true_surv1(k_grid-1) # this is what we're using in the sim now
right_truth <- test$true_surv1(k_grid-0.5)
mean((wrong_truth-mc_surv_est)^2)
mean((right_truth-mc_surv_est)^2)
