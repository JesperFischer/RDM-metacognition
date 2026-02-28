simulate_data_mcmc = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- (rnorm(1,-1,1))           # noise on evidence (needs to be fixed for now)
  bias <- rnorm(1,0,0.3)           # noise on evidence (needs to be fixed for now)
  
  sigma_m = (rnorm(1,-1,1))
  conf_prec = (rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  sigma_k <- (rnorm(1,-2,1))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.3,0.2)
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  # e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  
  # Compute posterior P(D=1|e)
  
  
  e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  c <- rnorm(n_trials, 0, exp(sigma_k))
  a_gauss <- as.numeric(e > c)  # Direct comparison
  ehat <- rnorm(n_trials, e, exp(sigma_m))
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat*X / (exp(sigma_e)^2 + exp(sigma_m)^2)),
                1-brms::inv_logit_scaled(2*ehat*X / (exp(sigma_e)^2 + exp(sigma_m)^2)))
  
  c_mu = ifelse(c_mu > 0.99, 0.99,ifelse(c_mu < 0.01, 0.01,c_mu))
  
  C = ordbetareg::rordbeta(n_trials,c_mu,conf_prec,c(c0,c1))
  
  # Dataframe
  df = data.frame(
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    bias = bias,
    slope = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2),
    sigma2 = (exp(sigma_e)^2 + exp(sigma_m)^2),
    # e = e,
    ehat = ehat,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    C = C,
    a_gauss = a_gauss
  )
}


simulate_data_mcmc_x1 = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- (rnorm(1,-1,1))           # noise on evidence (needs to be fixed for now)
  bias <- rnorm(1,0,0.3)           # noise on evidence (needs to be fixed for now)
  
  sigma_m = (rnorm(1,-1,1))
  conf_prec = (rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  sigma_k <- (rnorm(1,-2,1))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.3,0.2)
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  # e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  
  # Compute posterior P(D=1|e)
  
  
  e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  c <- rnorm(n_trials, 0, exp(sigma_k))
  a_gauss <- as.numeric(e > c)  # Direct comparison
  ehat <- rnorm(n_trials, e, exp(sigma_m))
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat*1 / (exp(sigma_e)^2 + exp(sigma_m)^2)),
                1-brms::inv_logit_scaled(2*ehat*1 / (exp(sigma_e)^2 + exp(sigma_m)^2)))
  
  c_mu = ifelse(c_mu > 0.99, 0.99,ifelse(c_mu < 0.01, 0.01,c_mu))
  
  C = ordbetareg::rordbeta(n_trials,c_mu,conf_prec,c(c0,c1))
  
  # Dataframe
  df = data.frame(
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    bias = bias,
    slope = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2),
    sigma2 = (exp(sigma_e)^2 + exp(sigma_m)^2),
    # e = e,
    ehat = ehat,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    C = C,
    a_gauss = a_gauss
  )
}
