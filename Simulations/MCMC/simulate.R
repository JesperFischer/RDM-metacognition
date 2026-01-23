pacman::p_load(tidyverse,posterior,bayesplot,cmdstanr,patchwork,ordbetareg,shiny,statConfR)

simulate_data = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- abs(rnorm(1,1,2))           # noise on evidence (needs to be fixed for now)
  sigma_m = abs(rnorm(1,1,2))
  sigma_k <- abs(rnorm(1,1,2))                 # criterion sd on the evidence scale (gaussian distribution)
  
  conf_prec = abs(rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  mean_choice <- rnorm(1,0,0.05)
  mu_k = 0
  
  # Coherence level (stimulus intensity)
  X <- abs(rnorm(n_trials,0,1))
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  e <- rnorm(n_trials, mean = mean_choice + D * X, sd = sigma_e)
  
  # Compute posterior P(D=1|e)
  
  ehat = rnorm(n_trials,e,sigma_m)
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  p_gauss = pnorm((e-mu_k)/sigma_k)                  #  probability of response on trial t given the citerion
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat*X / (sigma_e^2 + sigma_m^2)),
                1-brms::inv_logit_scaled(2*ehat*X / (sigma_e^2 + sigma_m^2)))
  
  
  
  C = ordbetareg::rordbeta(n_trials,c_mu,conf_prec,c(c0,c1))
  
  # Dataframe
  df = data.frame(
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    mean_choice = mean_choice,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    e = e,
    ehat = ehat,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    C = C,
    a_gauss = a_gauss
  )
}  



simulate_data_small_space = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- abs(rnorm(1,0.5,0.5))           # noise on evidence (needs to be fixed for now)
  sigma_m = abs(rnorm(1,0.5,0.5))
  sigma_k <- abs(rnorm(1,0.5,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  
  conf_prec = abs(rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  mean_choice <- rnorm(1,0,0.05)
  mu_k = 0
  
  # Coherence level (stimulus intensity)
  X <- abs(rnorm(n_trials,0,1))
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  e <- rnorm(n_trials, mean = mean_choice + D * X, sd = sigma_e)
  
  # Compute posterior P(D=1|e)
  
  ehat = rnorm(n_trials,e,sigma_m)
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  p_gauss = pnorm((e-mu_k)/sigma_k)                  #  probability of response on trial t given the citerion
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat*X / (sigma_e^2 + sigma_m^2)),
                1-brms::inv_logit_scaled(2*ehat*X / (sigma_e^2 + sigma_m^2)))
  
  
  
  C = ordbetareg::rordbeta(n_trials,c_mu,conf_prec,c(c0,c1))
  
  # Dataframe
  df = data.frame(
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    mean_choice = mean_choice,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    e = e,
    ehat = ehat,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    C = C,
    a_gauss = a_gauss
  )
}  




simulate_data_smallest_space = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- abs(rnorm(1,2,1))           # noise on evidence (needs to be fixed for now)
  sigma_m = abs(rnorm(1,2,1))
  sigma_k <- abs(rnorm(1,0,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  
  conf_prec = abs(rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  mean_choice <- rnorm(1,0,0.05)
  mu_k = 0
  
  # Coherence level (stimulus intensity)
  X <- abs(rnorm(n_trials,0,1))
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  e <- rnorm(n_trials, mean = mean_choice + D * X, sd = sigma_e)
  
  # Compute posterior P(D=1|e)
  
  ehat = rnorm(n_trials,e,sigma_m)
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  p_gauss = pnorm((e-mu_k)/sigma_k)                  #  probability of response on trial t given the citerion
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat*X / (sigma_e^2 + sigma_m^2)),
                1-brms::inv_logit_scaled(2*ehat*X / (sigma_e^2 + sigma_m^2)))
  
  
  
  C = ordbetareg::rordbeta(n_trials,c_mu,conf_prec,c(c0,c1))
  
  # Dataframe
  df = data.frame(
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    mean_choice = mean_choice,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    e = e,
    ehat = ehat,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    C = C,
    a_gauss = a_gauss
  )
}  