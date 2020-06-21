source(here("moss-simulation-master/code_simulation/simulate_data.R"))
library(devtools)
load_all(here("MOSS"))
test <- simulate_data(n_sim = 1e6)
k_grid <- 1:10

mc_surv_est <- sapply(k_grid-1,function(k)mean(k<test$dat$T1))
mc_haz_est <- sapply(k_grid-1,function(k)mean(test$dat$T1==k)/mean(test$dat$T1>=k))

wrong_truth <- test$true_surv1(k_grid-1) # this is what we're using in the sim now
right_truth <- test$true_surv1(k_grid-0.5)
mean((wrong_truth-mc_surv_est)^2)
mean((right_truth-mc_surv_est)^2)
surv_obj <- survival_curve$new(t = c(1,k_grid+1), survival=c(1,right_truth))
haz_truth <- surv_obj$survival_to_hazard()$hazard

