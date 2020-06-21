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


perf <- df_metric[t!=1,list(mean_est=mean(estimate),
                        true_surv=mean(true_surv),
                        mse = mean((estimate-true_surv)^2),
                        mse_se = sd((estimate-true_surv)^2)/sqrt(.N),
                        bias = mean(estimate-true_surv),
                        bias2 = mean(estimate-true_surv)^2,
                        var = var(estimate),
                        mean_n_t = mean(n_t),
                        n_sim = .N,
                        EED = mean(abs(ED)),
                        EED2 = mean(ED2)),by=list(method,t)]


ggplot(perf,aes(x=t,y=mse,color=method))+geom_line()
ggplot(perf,aes(x=t,y=n_sim))+geom_line()

long <- melt(perf,id=c("method","t","true_surv"),measure=c("mse","bias2","var"))
ggplot(long[t<=10],aes(x=t,y=value/true_surv^2,color=variable))+geom_line()+facet_wrap(~method)


# shut down for memory
# closeCluster(cl)
# mpi.quit()
stopCluster(cl)
