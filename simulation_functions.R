## Simulation functions used in "Bayesian optimization for personalized dose-finding
## trials with combination therapies"
## 2023-10-26

#-----------------------------------------------------------------------------#
# Standard Algorithm
#-----------------------------------------------------------------------------#

get_fits_reg <- function(data_x, data_y,
                         d = darg(NULL, data_x), g = garg(list(mle=TRUE), data_y), dK = TRUE,
                         tmax = c(1, var(data_y)), xnew){
  # note that this fits a separable GP
  eps <- sqrt(.Machine$double.eps)
  gpi <- newGPsep(data_x,
                  data_y, d=rep(d$start, ncol(data_x)), g=g$start, dK=dK)
  mle <- mleGPsep(gpi, param="both", tmin=c(d$min, g$min), tmax=c(d$max, g$max))
  p <- predGPsep(gpi, xnew)
  deleteGPsep(gpi)
  bind_cols(xnew,
            mu = p$mean,
            sigma_hat = sqrt(diag(p$Sigma)),
            df = p$df)
}

gen_ei_reg <- function(fits){
  y_best <- min(fits$mu) #heuristic/plug-in approach
  
  tibble(
    x_1 = fits$x_1,
    x_2 = fits$x_2,
    ei = map2_dbl(
      fits$mu,
      fits$sigma_hat,
      function(m, s) {
        if (s == 0) return(0)
        gamma <- (y_best - m) / s
        phi <- pnorm(gamma)
        return(s * (gamma * phi + dnorm(gamma)))
      }
    )
  ) %>%
    bind_cols(y_best = y_best)
}

gen_next_eval_reg <- function(ei_data){
  ei_data %>%
    top_n(1, ei) %>%
    select(x_1, x_2, ei, y_best) %>%
    slice_sample(n = 1)
}


run_sim_reg <- function(simulation, iterations, reps,
                        delta = 0, min_stop_iteration = 15,
                        optimal_dose_combo_x_3_0,
                        optimal_dose_combo_x_3_1,
                        optimal_value_x_3_0,
                        optimal_value_x_3_1,
                        std_dev,
                        seed = 20190416){
  set.seed(seed)
  
  .std_dev <- std_dev
  .reps <- reps # parts per dose combo
  
  # randomly assign x_3
  .initial_x <- tibble(x_1 = c(0, 0, 0.5, 1, 1),
                       x_2 = c(0, 1, 0.5, 0, 1),
                       iteration = 0) %>%
    uncount(reps) %>%
    mutate(x_3 = rbinom(nrow(.), 1, 0.5))
  
  .initial_objective_function <- foreach(i = 1:nrow(.initial_x), .combine = 'rbind') %do% objective_function(.initial_x$x_1[i],
                                                                                                             .initial_x$x_2[i],
                                                                                                             .initial_x$x_3[i])
  
  .initial_evals <- tibble(y = rnorm(nrow(.initial_x), mean = as_vector(.initial_objective_function), sd = .std_dev)) %>%
    bind_cols(.initial_x)
  
  
  .new_x <- expand_grid(x_1 = seq(0, 1, 0.01),
                        x_2 = seq(0, 1, 0.01))
  
  ## add in all variables which will be added eventually
  .evals <- bind_cols(.initial_evals, ei = NA_real_, y_best = NA_real_,
                      ei_stop = NA_real_)
  
  .delta <- delta # expected improvement stopping threshold
  
  ## Get fits from single GP
  .yhat <- get_fits_reg(.evals %>% filter(!is.na(y)) %>% select(x_1, x_2) %>% as.matrix(),
                        .evals %>% filter(!is.na(y)) %>% select(y) %>% as.matrix(),
                        xnew = .new_x)
  
  # yields recommended dose and value of mean/sd at this point at end of fitting ignoring x_3
  # set recommended dose as minimizer of posterior mean
  .recommended_dose_0 <- .yhat %>%
    filter(mu == min(mu)) %>%
    distinct() %>%
    slice_sample(n = 1) %>%
    mutate(lb_95_rec = mu - qstudent_t(0.975, df = df)*sigma_hat,
           ub_95_rec = mu + qstudent_t(0.975, df = df)*sigma_hat) %>% 
    rename(x_1_rec = x_1,
           x_2_rec = x_2,
           mu_rec = mu,
           sigma_hat_rec = sigma_hat)
  
  # collect posterior samples of f at recommended dose
  
  .post_samps_f <- rstudent_t(3000,
                              df = .recommended_dose_0$df,
                              mu = .recommended_dose_0$mu_rec,
                              sigma = .recommended_dose_0$sigma_hat_rec)
  
  # calculate posterior squared error loss
  .rmse_x_3_0 <- sqrt(mean((.post_samps_f - optimal_value_x_3_0)^2))
  .rmse_x_3_1 <- sqrt(mean((.post_samps_f - optimal_value_x_3_1)^2))
  
  # calculate euclidean distance between recommended dose and optimal here
  .euc_dist_dose_x_3_0 <- sqrt((.recommended_dose_0$x_1_rec - optimal_dose_combo_x_3_0[1])^2 +
                                 (.recommended_dose_0$x_2_rec - optimal_dose_combo_x_3_0[2])^2)
  .euc_dist_dose_x_3_1 <- sqrt((.recommended_dose_0$x_1_rec - optimal_dose_combo_x_3_1[1])^2 +
                                 (.recommended_dose_0$x_2_rec - optimal_dose_combo_x_3_1[2])^2)

  .performance_data <- bind_rows(tibble(x_3 = 0,
                                        rmse = .rmse_x_3_0,
                                        euc_dist_dose = .euc_dist_dose_x_3_0),
                                        #post_samps_f = list(.post_samps_f)), #can keep posterior samples if desired
                                 tibble(x_3 = 1,
                                        rmse = .rmse_x_3_1,
                                        euc_dist_dose = .euc_dist_dose_x_3_1))
                                        #post_samps_f = list(.post_samps_f)))
  
  .evals <- .recommended_dose_0 %>%
    uncount(2) %>%
    mutate(x_3 = c(0,1)) %>%
    right_join(.evals, by = "x_3") %>%
    left_join(.performance_data, by = "x_3")
  
  for(i in 1:iterations){
    .ei_data <- gen_ei_reg(.yhat %>% distinct())
    .next_eval <- gen_next_eval_reg(.ei_data) %>%
      uncount(.reps) %>%
      bind_cols(tibble(x_3 = rbinom(.reps, 1, 0.5)))
    
    .previous_ei_stop <- .evals %>% filter(iteration == i - 1) %>%
      pull(ei_stop) %>% unique()
    
    
    if(!is.na(.previous_ei_stop) & .previous_ei_stop == 1){
      .next_eval$ei <- 0
    }
    
    # allows stopping AFTER min iteration
    .ei_stop <- unique(.next_eval$ei)  < .delta  & i > min_stop_iteration
    
    if(.ei_stop == TRUE){
      .next_y_obs <- tibble(y = NA_real_,
                            ei_stop = as.numeric(.ei_stop))
    } else {
      .next_f_vals <- foreach(j=1:nrow(.next_eval), .combine = 'rbind') %do% {
        objective_function(.next_eval$x_1[j],
                           .next_eval$x_2[j],
                           .next_eval$x_3[j])
      }
      
      .next_y_obs <- tibble(y = rnorm(.reps, mean = as_vector(.next_f_vals),
                                      sd = .std_dev),
                            ei_stop = 0)
    }
    
    .new_eval_rows <- bind_cols(tibble(iteration = rep(i, .reps)),
                                .next_eval,
                                .next_y_obs)
    
    .evals_temp <- bind_rows(.evals, .new_eval_rows)
    
    .continue <- .new_eval_rows  %>% select(ei_stop) %>%
      distinct() %>% pull() == 0
    
    ## performs model fitting only if exploration should continue
    # otherwise skips to end and records the stop
    if(.continue == TRUE){
      
      # yields recommended dose and value of mean/sd at this point at end of fitting
      .yhat <- get_fits_reg(.evals_temp %>% filter(!is.na(y)) %>% select(x_1, x_2) %>% as.matrix(),
                            .evals_temp %>% filter(!is.na(y)) %>% select(y) %>% as.matrix(),
                            xnew = .new_x) 
      
      .new_eval_rows <- .yhat %>%
        filter(mu == min(mu)) %>%
        distinct() %>%
        slice_sample(n = 1) %>%
        mutate(lb_95_rec = mu - qstudent_t(0.975, df = df)*sigma_hat,
               ub_95_rec = mu + qstudent_t(0.975, df = df)*sigma_hat) %>% 
        rename(x_1_rec = x_1,
               x_2_rec = x_2,
               mu_rec = mu,
               sigma_hat_rec = sigma_hat) %>%
        uncount(2) %>%
        mutate(x_3 = c(0,1)) %>%
        right_join(.new_eval_rows, by = "x_3")
      
      .post_samps_f <- rstudent_t(3000,
                                   df = unique(.new_eval_rows$df),
                                   mu = unique(.new_eval_rows$mu_rec),
                                   sigma = unique(.new_eval_rows$sigma_hat_rec))
      
      # calculate posterior squared error loss (can call it rmse) for f here
      .rmse_x_3_0 <- sqrt(mean((.post_samps_f - optimal_value_x_3_0)^2))
      .rmse_x_3_1 <- sqrt(mean((.post_samps_f - optimal_value_x_3_1)^2))
      
      # calculate euclidean distance between recommended dose and optimal here
      .euc_dist_dose_x_3_0 <- sqrt((unique(.new_eval_rows$x_1_rec) - optimal_dose_combo_x_3_0[1])^2 +
                                     (unique(.new_eval_rows$x_2_rec) - optimal_dose_combo_x_3_0[2])^2)
      .euc_dist_dose_x_3_1 <- sqrt((unique(.new_eval_rows$x_1_rec) - optimal_dose_combo_x_3_1[1])^2 +
                                     (unique(.new_eval_rows$x_2_rec) - optimal_dose_combo_x_3_1[2])^2)
      
      #do for x_3_1
      
      .performance_data <- bind_rows(tibble(x_3 = 0,
                                            rmse = .rmse_x_3_0,
                                            euc_dist_dose = .euc_dist_dose_x_3_0),
                                            #post_samps_f = list(.post_samps_f)),
                                     tibble(x_3 = 1,
                                            rmse = .rmse_x_3_1,
                                            euc_dist_dose = .euc_dist_dose_x_3_1))
                                            #post_samps_f = list(.post_samps_f)))
      
      
      .new_eval_rows <- .new_eval_rows %>%
        left_join(.performance_data, by = "x_3")
      
      
    }
    
    .evals <- .evals %>%
      bind_rows(.new_eval_rows)
    
    if(.continue == FALSE){
      break
    }
    
  }

  bind_cols(tibble(simulation = rep(simulation, nrow(.evals)),
                   delta = rep(delta, nrow(.evals)),
                   reps = rep(reps, nrow(.evals)),
                   max_iterations = rep(iterations, nrow(.evals)),
                   min_stop_iteration = rep(min_stop_iteration, nrow(.evals)),
                   optimal_value_x_3_0 = rep(optimal_value_x_3_0, nrow(.evals)),
                   optimal_value_x_3_1 = rep(optimal_value_x_3_1, nrow(.evals)),
                   seed = rep(seed, nrow(.evals))),
            .evals)
  
}


#-----------------------------------------------------------------------------#
# Personalized Algorithm
#-----------------------------------------------------------------------------#

get_fits_pers <- function(data_x, data_y,
                          d = darg(NULL, data_x), g = garg(list(mle=TRUE), data_y), dK = TRUE,
                          tmax = c(1, var(data_y)), xnew){
  # note that this fits a separable GP
  eps <- sqrt(.Machine$double.eps)
  gpi <- newGPsep(data_x,
                  data_y, d=rep(d$start, ncol(data_x)), g=g$start, dK=dK)
  mle <- mleGPsep(gpi, param="both", tmin=c(d$min, g$min), tmax=c(d$max, g$max))
  p <- predGPsep(gpi, xnew)
  deleteGPsep(gpi)
  bind_cols(xnew,
            mu = p$mean,
            sigma_hat = sqrt(diag(p$Sigma)),
            df = p$df)
}


gen_ei_pers <- function(fits){
  y_best <- min(fits$mu) #heuristic approach
  
  tibble(
    x_1 = fits$x_1,
    x_2 = fits$x_2,
    x_3 = fits$x_3,
    ei = map2_dbl(
      fits$mu,
      fits$sigma_hat,
      function(m, s) {
        if (s == 0) return(0)
        gamma <- (y_best - m) / s
        phi <- pnorm(gamma)
        return(s * (gamma * phi + dnorm(gamma)))
      }
    )
  ) %>%
    bind_cols(y_best = y_best)
}

gen_next_eval_pers <- function(ei_data){
  ei_data %>%
    top_n(1, ei) %>%
    select(x_1, x_2, x_3, ei, y_best) %>%
    slice_sample(n = 1)
  
}

run_sim_pers <- function(simulation, iterations, reps,
                         delta = 0, min_stop_iteration = 15,
                         optimal_dose_combo_x_3_0,
                         optimal_dose_combo_x_3_1,
                         optimal_value_x_3_0,
                         optimal_value_x_3_1,
                         std_dev,
                         seed = 20190416){
  set.seed(seed)
  
  .std_dev <- std_dev
  .reps <- reps # parts per dose combo
  
  .initial_x <- tibble(x_1 = rep(c(0, 0, 0.5, 1, 1),2),
                       x_2 = rep(c(0, 1, 0.5, 0, 1),2),
                       x_3 = c(rep(0,5), rep(1,5)),
                       iteration = 0) %>% 
    uncount(.reps)
  
  .initial_objective_function <- foreach(i = 1:nrow(.initial_x), .combine = 'rbind') %do% objective_function(.initial_x$x_1[i],
                                                                                                             .initial_x$x_2[i],
                                                                                                             .initial_x$x_3[i])
  
  .initial_evals <- tibble(y = rnorm(nrow(.initial_x), mean = as_vector(.initial_objective_function), sd = .std_dev)) %>%
    bind_cols(.initial_x)
  
  .new_x <- expand_grid(x_1 = seq(0, 1, 0.01),
                        x_2 = seq(0, 1, 0.01),
                        x_3 = c(0,1))
  
  ## add in all variables which will be added eventually
  .evals <- bind_cols(.initial_evals, ei = NA_real_, y_best = NA_real_,
                      ei_stop = NA_real_)
  .delta <- delta
  
  ## Initial fit of data at iteration 0
  ## Get fits from single GP to borrow strength across strata
  ## have to split predictions at xnew to reduce computational burden 
  .yhat_x_3_0 <- get_fits_pers(.evals %>% filter(!is.na(y)) %>% select(x_1, x_2, x_3) %>% as.matrix(),
                               .evals %>% filter(!is.na(y)) %>% select(y) %>% as.matrix(),
                               xnew = .new_x %>% filter(x_3 == 0))
  .yhat_x_3_1 <- get_fits_pers(.evals %>% filter(!is.na(y)) %>% select(x_1, x_2, x_3) %>% as.matrix(),
                               .evals %>% filter(!is.na(y)) %>% select(y) %>% as.matrix(),
                               xnew = .new_x %>% filter(x_3 == 1))
  
  .yhat <- .yhat_x_3_0 %>%
    bind_rows(.yhat_x_3_1)
  
  # yields recommended dose and value of mean/sd at this point at end of fitting
  .recommended_dose_0 <- .yhat %>%
    group_by(x_3) %>%
    filter(mu == min(mu)) %>%
    distinct() %>%
    slice_sample(n = 1) %>%
    mutate(lb_95_rec = mu - qstudent_t(0.975, df = df)*sigma_hat,
           ub_95_rec = mu + qstudent_t(0.975, df = df)*sigma_hat) %>% 
    rename(x_1_rec = x_1,
           x_2_rec = x_2,
           mu_rec = mu,
           sigma_hat_rec = sigma_hat) %>%
    ungroup()

  .post_samps_f_x_3_0 <- rstudent_t(3000,
                                    df = .recommended_dose_0 %>% filter(x_3 == 0) %>% pull(df),
                                    mu = .recommended_dose_0 %>% filter(x_3 == 0) %>% pull(mu_rec),
                                    sigma = .recommended_dose_0 %>% filter(x_3 == 0) %>% pull(sigma_hat_rec))
  .post_samps_f_x_3_1 <- rstudent_t(3000,
                                    df = .recommended_dose_0 %>% filter(x_3 == 1) %>% pull(df),
                                    mu = .recommended_dose_0 %>% filter(x_3 == 1) %>% pull(mu_rec),
                                    sigma = .recommended_dose_0 %>% filter(x_3 == 1) %>% pull(sigma_hat_rec))
  
  # calculate posterior squared error loss 
  .rmse_x_3_0 <- sqrt(mean((.post_samps_f_x_3_0 - optimal_value_x_3_0)^2))
  .rmse_x_3_1 <- sqrt(mean((.post_samps_f_x_3_1 - optimal_value_x_3_1)^2))
  
  # calculate euclidean distance between recommended dose and optimal here
  .opt_dose_x_3_0 <- .recommended_dose_0 %>%
    filter(x_3 == 0)
  .opt_dose_x_3_1 <- .recommended_dose_0 %>%
    filter(x_3 == 1)
  .euc_dist_dose_x_3_0 <- sqrt((.opt_dose_x_3_0$x_1_rec - optimal_dose_combo_x_3_0[1])^2 +
                                 (.opt_dose_x_3_0$x_2_rec - optimal_dose_combo_x_3_0[2])^2)
  .euc_dist_dose_x_3_1 <- sqrt((.opt_dose_x_3_1$x_1_rec - optimal_dose_combo_x_3_1[1])^2 +
                                 (.opt_dose_x_3_1$x_2_rec - optimal_dose_combo_x_3_1[2])^2)
  
  
  .performance_data <- bind_rows(tibble(x_3 = 0,
                                        rmse = .rmse_x_3_0,
                                        euc_dist_dose = .euc_dist_dose_x_3_0),
                                        #post_samps_f = list(.post_samps_f_x_3_0)),
                                 tibble(x_3 = 1,
                                        rmse = .rmse_x_3_1,
                                        euc_dist_dose = .euc_dist_dose_x_3_1))
                                        #post_samps_f = list(.post_samps_f_x_3_1)))
  
  .evals <- .recommended_dose_0 %>%
    left_join(.evals, by = "x_3") %>%
    left_join(.performance_data, by = "x_3")

  for(i in 1:iterations){
    .ei_data_x3_0 <- gen_ei_pers(.yhat %>% filter(x_3 == 0) %>% distinct())
    .ei_data_x3_1 <- gen_ei_pers(.yhat %>% filter(x_3 == 1) %>% distinct())
    .next_eval_x3_0 <- gen_next_eval_pers(.ei_data_x3_0)
    .next_eval_x3_1 <- gen_next_eval_pers(.ei_data_x3_1)
    
    
    .previous_ei_stop_x3_0 <- .evals %>% filter(x_3 == 0, iteration == i - 1) %>%
      pull(ei_stop) %>% unique()
    .previous_ei_stop_x3_1 <- .evals %>% filter(x_3 == 1, iteration == i - 1) %>%
      pull(ei_stop) %>% unique()
    
    if(!is.na(.previous_ei_stop_x3_0) & .previous_ei_stop_x3_0 == 1){
      .next_eval_x3_0$ei <- 0
    }
    
    if(!is.na(.previous_ei_stop_x3_1) & .previous_ei_stop_x3_1 == 1){
      .next_eval_x3_1$ei <- 0
    }
    
    # allows stopping AFTER min iteration
    .ei_stop_x3_0 <- .next_eval_x3_0$ei  < .delta  & i > min_stop_iteration
    
    if(.ei_stop_x3_0 == TRUE){
      .next_y_obs_x3_0 <- tibble(y = NA_real_,
                                 ei_stop = as.numeric(.ei_stop_x3_0))
    } else {
      .next_y_obs_x3_0 <- tibble(y = rnorm(.reps, mean =
                                             as_vector(objective_function(.next_eval_x3_0$x_1,
                                                                          .next_eval_x3_0$x_2,
                                                                          .next_eval_x3_0$x_3)),
                                           sd = .std_dev),
                                 ei_stop = 0)
    }
    
    .ei_stop_x3_1 <- .next_eval_x3_1$ei < .delta & i > min_stop_iteration
    if(.ei_stop_x3_1 == TRUE) {
      .next_y_obs_x3_1 <- tibble(y = NA_real_,
                                 ei_stop = as.numeric(.ei_stop_x3_1))
    } else {
      .next_y_obs_x3_1 <- tibble(y = rnorm(.reps, mean =
                                             as_vector(objective_function(.next_eval_x3_1$x_1,
                                                                          .next_eval_x3_1$x_2,
                                                                          .next_eval_x3_1$x_3)),
                                           sd = .std_dev),
                                 ei_stop = 0)
    }
    
    .x3_0_rows <- bind_cols(tibble(iteration = rep(i, .reps)),
                            .next_eval_x3_0 %>% uncount(.reps),
                            .next_y_obs_x3_0)
    .x3_1_rows <- bind_cols(tibble(iteration = rep(i, .reps)),
                            .next_eval_x3_1 %>% uncount(.reps),
                            .next_y_obs_x3_1)
    .new_eval_rows <- bind_rows(.x3_0_rows,
                                .x3_1_rows)
    .evals_temp <- bind_rows(.evals, .new_eval_rows)
    
    
    .continue_x3_0 <- .x3_0_rows  %>% select(ei_stop) %>%
      distinct() %>% pull() == 0
    .continue_x3_1 <- .x3_1_rows  %>% select(ei_stop) %>%
      distinct() %>% pull() == 0
    
    
    ## performs model fitting only if exploration should continue
    # otherwise skips to end and records the stop
    
    if(.continue_x3_0 + .continue_x3_1 >= 1) {
      
      ## Get fits from single GP to borrow strength across strata
      .yhat_x_3_0 <- get_fits_pers(.evals_temp %>% filter(!is.na(y)) %>% select(x_1, x_2, x_3) %>% as.matrix(),
                                   .evals_temp %>% filter(!is.na(y)) %>% select(y) %>% as.matrix(),
                                   xnew = .new_x %>% filter(x_3 == 0))
      .yhat_x_3_1 <- get_fits_pers(.evals_temp %>% filter(!is.na(y)) %>% select(x_1, x_2, x_3) %>% as.matrix(),
                                   .evals_temp %>% filter(!is.na(y)) %>% select(y) %>% as.matrix(),
                                   xnew = .new_x %>% filter(x_3 == 1))
      .yhat <- .yhat_x_3_0 %>%
        bind_rows(.yhat_x_3_1)
      
      if(.continue_x3_0 == TRUE & .continue_x3_1 == TRUE){
        
        .new_eval_rows <- .yhat %>%
          group_by(x_3) %>%
          filter(mu == min(mu)) %>%
          distinct() %>%
          slice_sample(n = 1) %>%
          mutate(lb_95_rec = mu - qstudent_t(0.975, df = df)*sigma_hat,
                 ub_95_rec = mu + qstudent_t(0.975, df = df)*sigma_hat) %>% 
          rename(x_1_rec = x_1,
                 x_2_rec = x_2,
                 mu_rec = mu,
                 sigma_hat_rec = sigma_hat) %>%
          ungroup() %>%
          left_join(.new_eval_rows, by = "x_3")
        
      } else if (.continue_x3_1 == TRUE){
        .new_eval_rows <- .yhat %>%
          filter(x_3 == 1) %>%
          filter(mu == min(mu)) %>%
          distinct() %>%
          slice_sample(n = 1) %>%
          mutate(lb_95_rec = mu - qstudent_t(0.975, df = df)*sigma_hat, 
                 ub_95_rec = mu + qstudent_t(0.975, df = df)*sigma_hat) %>% 
          rename(x_1_rec = x_1,
                 x_2_rec = x_2,
                 mu_rec = mu,
                 sigma_hat_rec = sigma_hat) %>%
          ungroup() %>%
          add_row(x_3 = 0) %>%
          left_join(.new_eval_rows, by = "x_3")
      } else if (.continue_x3_0 == TRUE){
        .new_eval_rows <- .yhat %>%
          filter(x_3 == 0) %>%
          filter(mu == min(mu)) %>%
          distinct() %>%
          slice_sample(n = 1) %>%
          mutate(lb_95_rec = mu - qstudent_t(0.975, df = df)*sigma_hat, 
                 ub_95_rec = mu + qstudent_t(0.975, df = df)*sigma_hat) %>% 
          rename(x_1_rec = x_1,
                 x_2_rec = x_2,
                 mu_rec = mu,
                 sigma_hat_rec = sigma_hat) %>%
          ungroup() %>%
          add_row(x_3 = 1) %>%
          left_join(.new_eval_rows, by = "x_3")
      }
      
      .post_samps_f_x_3_0 <- rstudent_t(3000,
                                   df = .new_eval_rows %>% filter(x_3 == 0) %>% pull(df) %>% unique(),
                                   mu =  .new_eval_rows %>% filter(x_3 == 0) %>% pull(mu_rec) %>% unique(),
                                   sigma = .new_eval_rows  %>% filter(x_3 == 0) %>% pull(sigma_hat_rec) %>% unique())
      .post_samps_f_x_3_1 <-  rstudent_t(3000,
                                    df = .new_eval_rows %>% filter(x_3 == 1) %>% pull(df) %>% unique(),
                                    mu =  .new_eval_rows %>% filter(x_3 == 1) %>% pull(mu_rec) %>% unique(),
                                    sigma = .new_eval_rows  %>% filter(x_3 == 1) %>% pull(sigma_hat_rec) %>% unique())
      
      # calculate posterior squared error loss 
      .rmse_x_3_0 <- sqrt(mean((.post_samps_f_x_3_0 - optimal_value_x_3_0)^2))
      .rmse_x_3_1 <- sqrt(mean((.post_samps_f_x_3_1 - optimal_value_x_3_1)^2))
      
      # calculate euclidean distance between recommended dose and optimal here
      .opt_dose_x_3_0 <- .new_eval_rows %>%
        filter(x_3 == 0) %>% select(x_1_rec, x_2_rec) %>% distinct()
      .opt_dose_x_3_1 <- .new_eval_rows %>%
        filter(x_3 == 1) %>% select(x_1_rec, x_2_rec) %>% distinct()
      .euc_dist_dose_x_3_0 <- sqrt((.opt_dose_x_3_0$x_1_rec - optimal_dose_combo_x_3_0[1])^2 +
                                     (.opt_dose_x_3_0$x_2_rec - optimal_dose_combo_x_3_0[2])^2)
      .euc_dist_dose_x_3_1 <- sqrt((.opt_dose_x_3_1$x_1_rec - optimal_dose_combo_x_3_1[1])^2 +
                                     (.opt_dose_x_3_1$x_2_rec - optimal_dose_combo_x_3_1[2])^2)
      
      
      .performance_data <- bind_rows(tibble(x_3 = 0,
                                            rmse = .rmse_x_3_0,
                                            euc_dist_dose = .euc_dist_dose_x_3_0),
                                            #post_samps_f = list(.post_samps_f_x_3_0)),
                                     tibble(x_3 = 1,
                                            rmse = .rmse_x_3_1,
                                            euc_dist_dose = .euc_dist_dose_x_3_1))
                                            #post_samps_f = list(.post_samps_f_x_3_1)))
      
      .new_eval_rows <- .new_eval_rows %>%
        left_join(.performance_data, by = "x_3")
      
    }
    
    .evals <- .evals %>%
      bind_rows(.new_eval_rows)
    
    if(.continue_x3_0 == FALSE & .continue_x3_1 == FALSE){
      break
    }
    
  }

  bind_cols(tibble(simulation = rep(simulation, nrow(.evals)),
                   delta = rep(delta, nrow(.evals)),
                   reps = rep(reps, nrow(.evals)),
                   max_iterations = rep(iterations, nrow(.evals)),
                   min_stop_iteration = rep(min_stop_iteration, nrow(.evals)),
                   optimal_value_x_3_0 = rep(optimal_value_x_3_0, nrow(.evals)),
                   optimal_value_x_3_1 = rep(optimal_value_x_3_1, nrow(.evals)),
                   seed = rep(seed, nrow(.evals))),
            .evals)
}
