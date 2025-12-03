

fit_model = function(df, model, samples){

  
  datastan = list(N = nrow(df),
                  binom_y = df$Correct,
                  RT = df$RT,
                  Conf = df$Confidence,
                  X = df$X,
                  minRT = min(df$RT),
                  ACC = df$Correct)
  
  
  fit = model$sample(data = datastan,
                    iter_warmup = samples,
                    iter_sampling = samples,
                    parallel_chains = 4,
                    adapt_delta = 0.90)
}


get_model_predictions_model_1 = function(fit){
  
  
  n_draws = 2000
  workers = 2
  memory = 1000 * 1024^2
  
  # Group-level parameters (from gm)
  parameters = c("alpha", "beta", "lapse",
                 "rt_int", "rt_slope", "rt_prec",
                 "conf_prec", "meta_un", "meta_bias")
  
  # Extract group means (gm)
  df_param = as_draws_df(fit$draws("gm")) %>%
    select(-contains(".")) %>%
    rename_with(~parameters) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw, names_to = "variable")
  
  # Get average values for subject-specific parameters
  constants = as_draws_df(fit$draws(c("rt_ndt", "c0", "c11",
                                      "rho_p_rt", "rho_p_conf", "rho_rt_conf"))) %>%
    select(-contains(".")) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw) %>%
    extract(name, into = c("variable"),
            regex = "([a-zA-Z0-9_]+)", convert = TRUE) %>%
    group_by(variable) %>%
    summarize(mean = mean(value)) %>%
    pivot_wider(names_from = variable, values_from = mean)
  
  
  library(future.apply)
  # Set up parallel processing
  plan(multisession, workers = 2)
  options(future.globals.maxSize = 5000 * 1024^2)
  
  draws <- 1:n_draws
  
  # Only use the number of draws that the user wants
  dfq = df_param %>% filter(draw %in% draws)
  
  
  # Function to get the draws
  pred_list <- future_lapply(draws, function(d) {
    
    # Extract parameter vectors for this draw
    params <- dfq %>%
      filter(draw == d) %>%
      select(variable, value) %>%
      pivot_wider(names_from = "variable", values_from = "value")
    
    # Generate X values
    x = seq(-1, 1, by = 0.01)
    n_trials = length(x)
    
    # Get probability correct for each trial (using same approach as subject-level)
    prob_faster = psycho(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse) / 2)
    prob_cor = get_prob_cor(prob_faster, x)
    
    # Entropy for RT model
    entropy_t = entropy(prob_faster)
    
    # Theta for confidence (with meta_un)
    prob_faster_conf = psycho(x, params$alpha, exp(params$beta + params$meta_un), brms::inv_logit_scaled(params$lapse) / 2)
    
    # Build correlation matrix from averaged copula parameters
    R = matrix(c(1, constants$rho_p_rt, constants$rho_p_conf,
                 constants$rho_p_rt, 1, constants$rho_rt_conf,
                 constants$rho_p_conf, constants$rho_rt_conf, 1),
               nrow = 3, byrow = TRUE)
    
    # Sample multivariate normals (same as subject-level code)
    z_samples = MASS::mvrnorm(n = n_trials, mu = rep(0, 3), Sigma = R)
    
    # Transform to uniform [0,1] via standard normal CDF
    u_samples = pnorm(z_samples)
    
    # Transform uniforms to marginal distributions
    ACC_pred = rbinom(length(prob_cor), 1, prob_cor)
    
    # 2. RT (lognormal)
    rt_mu = params$rt_int + params$rt_slope * entropy_t
    RT_pred = qlnorm(u_samples[, 2], meanlog = rt_mu, sdlog = exp(params$rt_prec)) + constants$rt_ndt
    
    # Calculate expected RT mean
    rt_mu_expected = exp(rt_mu + exp(params$rt_prec)^2 / 2) + constants$rt_ndt
    
    # 3 Confidence mean (probability of getting it correct from confidence)
    prob_cor_conf = get_conf(ACC_pred, prob_faster_conf, x, params$alpha)
    
    # Apply meta_bias in logit space
    conf_mu_correct = brms::inv_logit_scaled(qlogis(prob_cor_conf) + params$meta_bias)
    
    # Sample confidence values
    conf_pred_correct = qordbeta(u_samples[, 3],
                                 mu = conf_mu_correct,
                                 phi = exp(params$conf_prec),
                                 cutzero = constants$c0,
                                 cutone = exp(constants$c0) + constants$c11)
    
    predictions = data.frame(
      X = x,
      Correct = ACC_pred,  # 1 if correct, 0 if incorrect
      prob = prob_cor,
      RT_pred = RT_pred,
      conf_mu_actual = conf_mu_correct,
      rt_mu = rt_mu_expected,
      Confidence = conf_pred_correct,
      prob_faster = prob_faster,
      draw = d
    )
    
    return(predictions)
  }, future.seed = TRUE)
  
  
  # Flatten nested list and create a tidy long dataframe
  predictions = map_dfr(pred_list, bind_rows)
  
}
