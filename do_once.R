moss_estimates <- function(moss_fit){
  psi <- colMeans(moss_fit$density_failure$survival)
  IC <- eic$new(
    A = moss_fit$A,
    T_tilde = moss_fit$T_tilde,
    Delta = moss_fit$Delta,
    density_failure = moss_fit$density_failure,
    density_censor = moss_fit$density_censor,
    g1W = moss_fit$g1W,
    psi = psi,
    A_intervene = moss_fit$A_intervene
  )$all_t(k_grid = moss_fit$k_grid)
  result <- list(psi=psi,IC=IC)
  return(result)
}

format_results <- function(method,estimates, k_grid, n_t){
  ED <- colMeans(estimates$IC)
  ED2 <- colMeans(estimates$IC^2)
  results_df <- data.table(method = method,
                           t = k_grid,
                           estimate = estimates$psi,
                           ED = ED,
                           ED2 = ED2,
                           n_t=n_t)
}

do_once <- function(n_sim = 2e2) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1
  W_names <- c("W", "W1")

  # TODO: check
  while (all(df$Delta == 1)) {
    simulated <- simulate_data(n_sim = n_sim)
    df <- simulated$dat
    true_surv <- simulated$true_surv1
  }

  # TODO: check
  sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
  sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam")
  sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam")
  range(df$T.tilde)
  df$T.tilde <- df$T.tilde + 1
  k_grid <- 1:max(df$T.tilde)
  n_t <- sapply(k_grid, function(k)sum(df$T.tilde>=k))
  # message("KM")
  # n_sample <- nrow(df)
  # km_fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
  # surv1_km <- tail(km_fit$surv, km_fit$strata["A=1"])
  # time1_km <- tail(km_fit$time, km_fit$strata["A=1"])
  # surv0_km <- tail(km_fit$surv, km_fit$strata["A=0"])
  # time0_km <- tail(km_fit$time, km_fit$strata["A=0"])
  # library(zoo)
  # impute_KM <- function(time, km) {
  #   surv1_km_final <- rep(NA, max(df$T.tilde))
  #   surv1_km_final[time] <- km
  #   surv1_km_final <- na.locf(surv1_km_final, na.rm = FALSE)
  #   surv1_km_final[is.na(surv1_km_final)] <- 1
  #   surv1_km_final <- c(1, surv1_km_final)
  #   surv1_km_final <- surv1_km_final[-length(surv1_km_final)]
  #   return(surv1_km_final)
  # }
  # surv1_km_final <- impute_KM(time = time1_km, km = surv1_km)
  # surv0_km_final <- impute_KM(time = time0_km, km = surv0_km)
  # km_fit_1 <- survival_curve$new(t = k_grid, survival = surv1_km_final)
  # km_fit_0 <- survival_curve$new(t = k_grid, survival = surv0_km_final)

  message("SL")
  sl_fit <- initial_sl_fit(
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    A = df$A,
    W = data.frame(df[, W_names]),
    t_max = max(df$T.tilde),
    sl_treatment = sl_lib_g,
    sl_censoring = sl_lib_censor,
    sl_failure = sl_lib_failure
  )
  
  sl_fit$density_failure_1$hazard_to_survival()
  sl_fit$density_failure_0$hazard_to_survival()
  # WILSON hack no data is t_tilde = 2
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid
  
  
  
  
  message("moss")

  message("moss with l2 submodel")
  moss_hazard_l2 <- MOSS_hazard$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  
  sl_estimates <- moss_estimates(moss_hazard_l2) 
  sl_results <- format_results("super learner",sl_estimates,k_grid, n_t)
  
  
  moss_hazard_l1 <- moss_hazard_l2$clone(deep = TRUE)
  
  # TODO: check
  moss_hazard_l2$iterate_onestep(
    method = "l2", epsilon = 1e-1 / sqrt(n_sim), max_num_interation = 5e1, verbose = FALSE
  )
  
  moss_l2_estimates <- moss_estimates(moss_hazard_l2)
  moss_l2_results <- format_results("moss_l2", moss_l2_estimates, k_grid, n_t)
  
  # moss_hazard_l1$iterate_onestep(
  #   method = "l1", epsilon = 1e-1 / sqrt(n_sim), max_num_interation = 1e0, verbose = FALSE
  # )
  # 
  # moss_l1_estimates <- moss_estimates(moss_hazard_l1)
  # moss_l1_results <- format_results("moss_l1", moss_l1_estimates, k_grid, n_t)
  # 
  ################################################################################
  # TODO: check
  message("tlverse 1-dimensional ulfm")

  tmax <- max(df$T.tilde)
  all_times <- lapply(seq_len(tmax), function(t_current){
    df_time <- copy(df)
    # TODO: check
    df_time$N <- as.numeric(t_current == df$T.tilde & df$Delta == 1)
    df_time$A_c <- as.numeric(t_current == df$T.tilde & df$Delta == 0)
    df_time$pre_failure <- as.numeric(t_current<=df$T.tilde)
    df_time$t <- t_current

    return(df_time)
  })

  df_long <- rbindlist(all_times)

  node_list <- list(W = c("W", "W1"), A = "A", T_tilde = "T.tilde", Delta = "Delta", 
    time = "t", N = "N", A_c = "A_c", id ="ID", pre_failure = "pre_failure")

  
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_gam <- make_learner(Lrnr_gam)
  sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
  learner_list <- list(A = sl_A, N = sl_A, A_c = sl_A)

  # TODO: check
  var_types <- list(T_tilde = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), 
    Delta = Variable_Type$new("binomial"))
  survival_spec <- tmle_survival(treatment_level = 1, control_level = 0, 
                                 target_times = intersect(1:10, k_grid),
                                 variable_types = var_types)
  survival_task <- survival_spec$make_tmle_task(df_long, node_list)

  likelihood <- survival_spec$make_initial_likelihood(survival_task, learner_list)

  initial_likelihood <- likelihood
  # TODO: check
  # up <- tmle3_Update_survival$new(maxit = 1e2, clipping = 1e-1 / sqrt(n_sim))
  # up <- tmle3_Update_survival$new(
  #   one_dimensional = TRUE, constrain_step = TRUE,
  #   maxit = 1e2, cvtmle = TRUE,
  #   convergence_type = "sample_size",
  #   delta_epsilon = 1e-2,
  #   fit_method = "classic"
  # )
  # up <- tmle3_Update_survival$new(
  #     constrain_step = TRUE, one_dimensional = TRUE, 
  #     delta_epsilon = 3e-2, 
  #     verbose = TRUE,
  #     convergence_type = "scaled_var", 
  #     maxit = 5e1, 
  #     fit_method = "classic")
  up <- tmle3_Update$new(constrain_step = TRUE, one_dimensional = TRUE, 
                       delta_epsilon = 3e-2, verbose = FALSE,
                       convergence_type = "scaled_var", 
                       maxit = 5e1, use_best = TRUE)
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
  # targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  tmle_task <- survival_task
  tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)

  # TODO: initial
  ps <- tmle_params[[1]]
  tl_initial_estimates <- ps$estimates(tmle_task)
  tl_initial_results <- format_results("sl3", tl_initial_estimates, k_grid + 1, n_t)
  
  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  
  tl_ulfs_results <- format_results("tmle3 ulfs", tmle_fit_manual$estimates[[1]], k_grid + 1, n_t)


  message("tlverse l2")
  initial_likelihood <- likelihood
  

  up <- tmle3_Update_survival$new(
    # TODO: check
    # one_dimensional = TRUE, constrain_step = TRUE,
    maxit = 5e1, 
    # cvtmle = TRUE,
    convergence_type = "scaled_var",
    # delta_epsilon = 1e-2,
    fit_method = "l2",
    clipping = 1e-2,
    use_best = TRUE
  )
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
  # targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  tmle_task <- survival_task
  tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)

  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  
  tl_l2_results <- format_results("tmle3 l2", tmle_fit_manual$estimates[[1]], k_grid + 1, n_t)
  
  
  # evaluate against truth
  # TODO: fix the truth
  true_indx <- c(0, k_grid - 0.5)
  true_indx <- true_indx[seq(1, length(true_indx) - 1)]
  true_surv <- simulated$true_surv1(true_indx)
  true_surv_dt <- data.table(t=k_grid, true_surv=true_surv)
  
  results <- rbindlist(list(
    sl_results,
    # moss_l1_results,
    moss_l2_results,
    tl_initial_results,
    tl_ulfs_results,
    tl_l2_results
  ))
  
  results <- merge(results, true_surv_dt, by="t")
  # TODO: think about time 0, what these survivals are actually indicating
  
  
  
  return(results)
}

get_mean_eic <- function(name, df_metric) {
  idx <- which(df_metric$method == name)
  eic <- unique(df_metric$eic[idx])
  return(mean(eic, na.rm=TRUE))
}
