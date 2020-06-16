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

  # # ipcw
  # message("ipcw + ee")
  # ipcw_fit_1_all <- repeat_t_grid$new(
  #   method = ipcw,
  #   A = df$A,
  #   T_tilde = df$T.tilde,
  #   Delta = df$Delta,
  #   density_failure = sl_fit$density_failure_1,
  #   density_censor = sl_fit$density_censor_1,
  #   g1W = sl_fit$g1W,
  #   A_intervene = 1
  # )$fit(k_grid = k_grid)
  # ipcw_fit_0_all <- repeat_t_grid$new(
  #   method = ipcw,
  #   A = df$A,
  #   T_tilde = df$T.tilde,
  #   Delta = df$Delta,
  #   density_failure = sl_fit$density_failure_0,
  #   density_censor = sl_fit$density_censor_0,
  #   g1W = sl_fit$g1W,
  #   A_intervene = 0
  # )$fit(k_grid = k_grid)
  # ee_fit_1_all <- repeat_t_grid$new(
  #   method = ee,
  #   A = df$A,
  #   T_tilde = df$T.tilde,
  #   Delta = df$Delta,
  #   density_failure = sl_fit$density_failure_1,
  #   density_censor = sl_fit$density_censor_1,
  #   g1W = sl_fit$g1W,
  #   A_intervene = 1
  # )$fit(k_grid = k_grid)
  # ee_fit_0_all <- repeat_t_grid$new(
  #   method = ee,
  #   A = df$A,
  #   T_tilde = df$T.tilde,
  #   Delta = df$Delta,
  #   density_failure = sl_fit$density_failure_0,
  #   density_censor = sl_fit$density_censor_0,
  #   g1W = sl_fit$g1W,
  #   A_intervene = 0
  # )$fit(k_grid = k_grid)
  # ipcw_fit_1 <- survival_curve$new(t = k_grid, survival = ipcw_fit_1_all)
  # ipcw_fit_0 <- survival_curve$new(t = k_grid, survival = ipcw_fit_0_all)
  # ee_fit_1 <- survival_curve$new(t = k_grid, survival = ee_fit_1_all)
  # ee_fit_0 <- survival_curve$new(t = k_grid, survival = ee_fit_0_all)

  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)
  message("moss")
  moss_fit <- MOSS$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  # TODO: check
  psi_moss_1 <- moss_fit$onestep_curve(
    epsilon = 1e-1 / n_sim,
    max_num_interation = 1e0,
    verbose = FALSE
  )
  moss_fit <- MOSS$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0,
    k_grid = k_grid
  )
  # TODO: check
  psi_moss_0 <- moss_fit$onestep_curve(
    epsilon = 1e-1 / n_sim,
    max_num_interation = 1e0,
    verbose = FALSE
  )
  moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
  moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)

  # # tmle
  # message("tmle")
  # tmle_fit <- tryCatch({
  #   tmle_fit <- fit_survtmle(
  #     T.tilde = df$T.tilde,
  #     Delta = df$Delta,
  #     A = df$A,
  #     W_df = data.frame(df[, W_names]),
  #     SL.trt = sl_lib_g,
  #     SL.ctime = sl_lib_censor,
  #     SL.ftime = sl_lib_failure
  #   )
  # },
  # error = function(cond) {
  #   message("tmle error")
  #   NULL
  # }
  # )
  # if (is.null(tmle_fit)) {
  #   tmle_fit_1 <- sl_density_failure_1_marginal$clone(deep = TRUE)
  #   tmle_fit_0 <- sl_density_failure_0_marginal$clone(deep = TRUE)
  #   is_tmle1_converge <- FALSE
  # } else {
  #   s_1 <- c(1, tmle_fit$s_1)
  #   s_1 <- s_1[-length(s_1)]
  #   s_0 <- c(1, tmle_fit$s_0)
  #   s_0 <- s_0[-length(s_0)]
  #   tmle_fit_1 <- survival_curve$new(t = k_grid, survival = s_1)
  #   tmle_fit_0 <- survival_curve$new(t = k_grid, survival = s_0)
  #   is_tmle1_converge <- TRUE
  # }
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
  moss_hazard_l1 <- moss_hazard_l2$clone(deep = TRUE)
  # TODO: check
  moss_l2_1_rs <- moss_hazard_l2$iterate_onestep(
    method = "l2", epsilon = 1e-1 / sqrt(n_sim), max_num_interation = 5e1, verbose = FALSE
  )
  psi_moss_l2_1 <- moss_l2_1_rs$psi_n
  mean_eic_inner_product_moss_l2_1 <- min(moss_l2_1_rs$eic_list)

  moss_hazard_l2_1 <- survival_curve$new(t = k_grid, survival = psi_moss_l2_1)

  # TODO: check
  moss_l1_1_rs <- moss_hazard_l1$iterate_onestep(
    method = "l1", epsilon = 1e-1 / sqrt(n_sim), max_num_interation = 1e0, verbose = FALSE
  )
  psi_moss_hazard_l1_1 <- moss_l1_1_rs$psi_n
  mean_eic_inner_product_moss_l1_1 <- min(moss_l1_1_rs$eic_list)

  moss_hazard_l1_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_l1_1)

  ################################################################################
  # TODO: check
  message("tlverse 1-dimensional ulfm")
  # tlverse
  # tmax <- max(df$T.tilde)
  # all_times <- lapply(seq_len(tmax), function(t_current){
  #   df_time <- copy(df)
  #   # TODO: check
  #   df_time$N <- ifelse(t_current == df$T.tilde & df$Delta == 1, 1, 0)
  #   df_time$A_c <- ifelse(t_current == df$T.tilde & df$Delta == 0, 1, 0)
  #   df_time$t <- t_current

  #   return(df_time)
  # })
  # df_long <- rbindlist(all_times)

  # node_list <- list(W = c("W", "W1"), A = "A", T_tilde = "T.tilde", Delta = "Delta", 
  #   t = "t", N = "N", A_c = "A_c")
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

  # TODO: check
  # lrnr_xgb <- make_learner(Lrnr_xgboost)
  # learner_list <- list(A = lrnr_xgb, N = lrnr_xgb, A_c = lrnr_xgb)
  # lrnr_glm <- make_learner(Lrnr_glm)
  # learner_list <- list(A = lrnr_glm, N = lrnr_glm, A_c = lrnr_glm)
  # lrnr_mean <- make_learner(Lrnr_mean)
  # learner_list <- list(A = lrnr_mean, N = lrnr_mean, A_c = lrnr_mean)
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_gam <- make_learner(Lrnr_gam)
  # lrnr_earth <- make_learner(Lrnr_earth)
  # sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam, lrnr_earth))
  sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
  learner_list <- list(A = sl_A, N = sl_A, A_c = sl_A)

  # TODO: check
  var_types <- list(T_tilde = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), 
    Delta = Variable_Type$new("binomial"))
  survival_spec <- tmle_survival(treatment_level = 1, control_level = 0, variable_types = var_types)
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
                       convergence_type = "scaled_var", maxit = 5e1)
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
  # targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  tmle_task <- survival_task
  tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)

  # TODO: initial
  ps <- tmle_params[[1]]
  cf_task <- ps$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
  pN1 <- ps$observed_likelihood$get_likelihoods(cf_task, "N")
  time <- tmle_task$time
  id <- tmle_task$id
  pN1_mat <- ps$long_to_mat(pN1,id,time)
  SN1_mat <- ps$hm_to_sm(pN1_mat)
  psi_tl_initial <- colMeans(SN1_mat)
  # TODO : check
  psi_tl_initial <- c(1, psi_tl_initial[seq(1, length(psi_tl_initial) - 1)])
  psi_tl_initial <- survival_curve$new(t = k_grid, survival = psi_tl_initial)

  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  rs <- tmle_fit_manual$estimates[[1]]
  psi1_tl <- rs$psi
  # TODO : check
  psi1_tl <- c(1, psi1_tl[seq(1, length(psi1_tl) - 1)])
  psi1_tl <- survival_curve$new(t = k_grid, survival = psi1_tl)

  # TODO : check
  eic_tl <- rs$IC
  mean_eic_tl <- colMeans(eic_tl)
  mean_eic_inner_product_tl <- abs(sqrt(sum(mean_eic_tl ^ 2)))

  # TODO: l2 update
  initial_likelihood <- likelihood
  # TODO: check
  # up <- tmle3_Update_survival$new(maxit = 1e2, clipping = 1e-1 / sqrt(n_sim))
  up <- tmle3_Update_survival$new(
    # TODO: check
    # one_dimensional = TRUE, constrain_step = TRUE,
    maxit = 5e1, 
    # cvtmle = TRUE,
    # convergence_type = "sample_size",
    # delta_epsilon = 1e-2,
    fit_method = "l2",
    clipping = 1e-2
  )
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
  # targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  tmle_task <- survival_task
  tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)

  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  # TODO : check
  print("--------------")
  print("l2 step number")
  print(tmle_fit_manual$steps)
  print("--------------")
  rs <- tmle_fit_manual$estimates[[1]]
  psi_tl_l2 <- rs$psi
  # TODO : check
  psi_tl_l2 <- c(1, psi_tl_l2[seq(1, length(psi_tl_l2) - 1)])
  psi_tl_l2 <- survival_curve$new(t = k_grid, survival = psi_tl_l2)

  eic_tl_l2 <- rs$IC
  mean_eic_tl_l2 <- colMeans(eic_tl_l2)
  mean_eic_inner_product_tl_2 <- abs(sqrt(sum(mean_eic_tl_l2 ^ 2)))
  ################################################################################

  # evaluate against truth
  survival_truth_1 <- survival_curve$new(t = k_grid, survival = simulated$true_surv1(k_grid - 1))
  survival_truth_0 <- survival_curve$new(t = k_grid, survival = simulated$true_surv0(k_grid - 1))

  # TODO: check
  evaluate_tlverse_initial <- evaluate_metric$new(
    survival = psi_tl_initial, survival_truth = survival_truth_1
  )
  df_entropy_tlverse_initial <- evaluate_tlverse_initial$evaluate_cross_entropy()
  df_entropy_tlverse_initial$metric_name <- "cross_entropy"
  df_mse_tlverse_initial <- evaluate_tlverse_initial$evaluate_mse()
  df_mse_tlverse_initial$metric_name <- "mse"

  # TODO: check
  evaluate_tlverse_l2 <- evaluate_metric$new(
    survival = psi_tl_l2, survival_truth = survival_truth_1
  )
  df_entropy_tlverse_l2 <- evaluate_tlverse_l2$evaluate_cross_entropy()
  df_entropy_tlverse_l2$metric_name <- "cross_entropy"
  # TODO: check
  df_mse_tlverse_l2 <- evaluate_tlverse_l2$evaluate_mse(mean_eic_inner_product_tl_2)
  df_mse_tlverse_l2$metric_name <- "mse"

  # TODO: check
  evaluate_tlverse <- evaluate_metric$new(
    survival = psi1_tl, survival_truth = survival_truth_1
  )
  df_entropy_tlverse_1 <- evaluate_tlverse$evaluate_cross_entropy()
  df_entropy_tlverse_1$metric_name <- "cross_entropy"
  # TODO: check
  df_mse_tlverse_1 <- evaluate_tlverse$evaluate_mse(mean_eic_inner_product_tl)
  df_mse_tlverse_1$metric_name <- "mse"

  evaluate_moss <- evaluate_metric$new(
    survival = moss_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_1 <- evaluate_moss$evaluate_cross_entropy()
  df_entropy_moss_1$metric_name <- "cross_entropy"
  df_mse_moss_1 <- evaluate_moss$evaluate_mse()
  df_mse_moss_1$metric_name <- "mse"

  evaluate_moss_l2 <- evaluate_metric$new(
    survival = moss_hazard_l2_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_l2_1 <- evaluate_moss_l2$evaluate_cross_entropy()
  df_entropy_moss_l2_1$metric_name <- "cross_entropy"
  df_mse_moss_l2_1 <- evaluate_moss_l2$evaluate_mse(mean_eic_inner_product_moss_l2_1)
  df_mse_moss_l2_1$metric_name <- "mse"

  evaluate_moss_l1 <- evaluate_metric$new(
    survival = moss_hazard_l1_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_l1_1 <- evaluate_moss_l1$evaluate_cross_entropy()
  df_entropy_moss_l1_1$metric_name <- "cross_entropy"
  df_mse_moss_l1_1 <- evaluate_moss_l1$evaluate_mse(mean_eic_inner_product_moss_l1_1)
  df_mse_moss_l1_1$metric_name <- "mse"

  evaluate_sl <- evaluate_metric$new(
    survival = sl_density_failure_1_marginal, survival_truth = survival_truth_1
  )
  df_entropy_sl_1 <- evaluate_sl$evaluate_cross_entropy()
  df_entropy_sl_1$metric_name <- "cross_entropy"
  df_mse_sl_1 <- evaluate_sl$evaluate_mse()
  df_mse_sl_1$metric_name <- "mse"

  # evaluate_ipcw <- evaluate_metric$new(
  #   survival = ipcw_fit_1, survival_truth = survival_truth_1
  # )
  # df_entropy_ipcw_1 <- evaluate_ipcw$evaluate_cross_entropy()
  # df_entropy_ipcw_1$metric_name <- "cross_entropy"
  # df_mse_ipcw_1 <- evaluate_ipcw$evaluate_mse()
  # df_mse_ipcw_1$metric_name <- "mse"

  # evaluate_ee <- evaluate_metric$new(
  #   survival = ee_fit_1, survival_truth = survival_truth_1
  # )
  # df_entropy_ee_1 <- evaluate_ee$evaluate_cross_entropy()
  # df_entropy_ee_1$metric_name <- "cross_entropy"
  # df_mse_ee_1 <- evaluate_ee$evaluate_mse()
  # df_mse_ee_1$metric_name <- "mse"

  # evaluate_tmle <- evaluate_metric$new(
  #   survival = tmle_fit_1, survival_truth = survival_truth_1
  # )
  # df_entropy_tmle_1 <- evaluate_tmle$evaluate_cross_entropy()
  # df_entropy_tmle_1$metric_name <- "cross_entropy"
  # df_mse_tmle_1 <- evaluate_tmle$evaluate_mse()
  # df_mse_tmle_1$metric_name <- "mse"

  # evaluate_km <- evaluate_metric$new(
  #   survival = km_fit_1, survival_truth = survival_truth_1
  # )
  # df_entropy_km_1 <- evaluate_km$evaluate_cross_entropy()
  # df_entropy_km_1$metric_name <- "cross_entropy"
  # df_mse_km_1 <- evaluate_km$evaluate_mse()
  # df_mse_km_1$metric_name <- "mse"

  # TODO: check
  df_mse_tlverse_initial$method <- "sl3"
  df_mse_tlverse_l2$method <- "tmle3 l2"
  df_mse_tlverse_1$method <- "tmle3 1-dim ulfm"
  df_mse_moss_1$method <- "MOSS_classic"
  df_mse_moss_l2_1$method <- "MOSS_l2"
  df_mse_moss_l1_1$method <- "MOSS_l1"
  df_mse_sl_1$method <- "super learner"
  # df_mse_ipcw_1$method <- "IPCW"
  # df_mse_ee_1$method <- "EE"
  # df_mse_tmle_1$method <- "TMLE"
  # df_mse_km_1$method <- "KM"
  # TODO: check
  df_entropy_tlverse_initial$method <- "sl3"
  df_entropy_tlverse_l2$method <- "tmle3 l2"
  df_entropy_tlverse_1$method <- "tmle3 1-dim ulfm"
  df_entropy_moss_1$method <- "MOSS_classic"
  df_entropy_moss_l2_1$method <- "MOSS_l2"
  df_entropy_moss_l1_1$method <- "MOSS_l1"
  df_entropy_sl_1$method <- "super learner"
  # df_entropy_ipcw_1$method <- "IPCW"
  # df_entropy_ee_1$method <- "EE"
  # df_entropy_tmle_1$method <- "TMLE"
  # df_entropy_km_1$method <- "KM"
  df_plot <- plyr::rbind.fill(
    # TODO: check
    df_mse_tlverse_initial,
    df_mse_tlverse_l2,
    df_mse_tlverse_1,
    df_mse_moss_1,
    df_mse_moss_l2_1,
    df_mse_moss_l1_1,
    df_mse_sl_1,
    # df_mse_ipcw_1,
    # df_mse_ee_1,
    # df_mse_tmle_1,
    # df_mse_km_1,
    # TODO: check
    df_entropy_tlverse_initial,
    df_entropy_tlverse_l2,
    df_entropy_tlverse_1,
    df_entropy_moss_1,
    df_entropy_moss_l2_1,
    df_entropy_moss_l1_1,
    df_entropy_sl_1
    # df_entropy_ipcw_1,
    # df_entropy_ee_1,
    # df_entropy_tmle_1,
    # df_entropy_km_1
  )
  # track if the estimators are monotone
  # df_plot$is_monotone_tlverse1 <- all(diff(as.numeric(psi1_tl)) <= 0)
  # df_plot$is_monotone_tmle1 <- all(diff(as.numeric(tmle_fit_1$survival)) <= 0)
  # df_plot$is_monotone_ee1 <- all(diff(as.numeric(ee_fit_1$survival)) <= 0)
  # df_plot$is_monotone_tmle0 <- all(diff(as.numeric(tmle_fit_0$survival)) <= 0)
  # df_plot$is_monotone_ee0 <- all(diff(as.numeric(ee_fit_0$survival)) <= 0)
  # df_plot$is_tmle1_converge <- is_tmle1_converge
  return(df_plot)
}

get_mean_eic <- function(name, df_metric) {
  idx <- which(df_metric$method == name)
  eic <- unique(df_metric$eic[idx])
  return(mean(eic, na.rm=TRUE))
}
