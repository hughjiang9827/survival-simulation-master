# assertthat::assert_that(packageVersion("MOSS") >= "1.1.2")
# TODO: check path
devtools::load_all("Moss")
devtools::load_all("tmle3")
library(sl3)

library(survival)
# library(MOSS)
library(ggpubr)
library(tidyverse)
library(here)
library(foreach)
# source(here("fit_survtmle.R"))
source(here("Moss-simulation-master/fit_survtmle.R"))
# simulate data
# source(here("./code_simulation/simulate_data.R"))
source(here("Moss-simulation-master/code_simulation/simulate_data.R"))
# TODO: check
source(here("Moss-simulation-master/do_once.R"))

# TODO: check
set.seed(42)
N_SIMULATION <- 100
n_sim_grid <- c(1e2)
# n_sim_grid <- c(1e3, 1e2)
# n_sim_grid <- c(1e3, 5e2, 1e2)
df_metric <- foreach(
  n_sim = n_sim_grid,
  .combine = rbind,
  # TODO: check
  # .packages = c("R6", "MOSS", "survtmle", "survival"),
  .packages = c("R6", "survtmle", "survival"),
  .inorder = FALSE,
  .verbose = TRUE
) %:%
  foreach(
    it2 = 1:N_SIMULATION, .combine = rbind, .errorhandling = "remove"
  ) %dopar% {
    df <- do_once(n_sim = n_sim)
    df$id_mcmc <- it2
    df$n <- n_sim
    return(df)
  }
table(df_metric$id_mcmc)

save(df_metric, file = "Moss-simulation-master/code_simulation/df_metric.rda")


# shut down for memory
# closeCluster(cl)
# mpi.quit()
stopCluster(cl)
