
simulate_data_p = function(n = 1000, sigma_e,sigma_m,sigma_k,bias, seed = 1997){
  
  set.seed(seed)
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  # sigma_e <- (rnorm(1,-1,0.5))           # noise on evidence (needs to be fixed for now)
  # bias <- rnorm(1,0,0.2)           # noise on evidence (needs to be fixed for now)
  
  # sigma_m = (rnorm(1,-1,0.5))
  conf_prec = (rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  # sigma_k <- (rnorm(1,-2,0.25))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.5,1)
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  
  
  ehat = rnorm(n_trials,e,(exp(sigma_m)))
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  # k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  
  # p_gauss = pnorm((D * X - bias) / sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))                  #  probability of response on trial t given the citerion
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  # k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  
  p_gauss = pnorm(e / exp(sigma_k))                  #  probability of response on trial t given the citerion
  
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat*(X) / (exp(sigma_e)^2 + exp(sigma_m)^2)),
                1-brms::inv_logit_scaled(2*ehat*(X) / (exp(sigma_e)^2 + exp(sigma_m)^2)))
  
  # Dataframe
  df = data.frame(
    p_gauss = p_gauss,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    a_gauss = a_gauss
  ) %>% mutate(   c0 = c0,
                   c1 = c1,
                   sigma_k = sigma_k,
                   sigma_e = sigma_e,
                   conf_prec = conf_prec,
                   sigma_m = sigma_m,
                   bias = bias)
  
  
  
}

simulate_data_margin = function(n = 1000, sigma_e,sigma_m,sigma_k,bias, seed = 1997){
  
  set.seed(seed)
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  # sigma_e <- (rnorm(1,-1,0.5))           # noise on evidence (needs to be fixed for now)
  # bias <- rnorm(1,0,0.2)           # noise on evidence (needs to be fixed for now)
  
  # sigma_m = (rnorm(1,-1,0.5))
  conf_prec = (rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  # sigma_k <- (rnorm(1,-2,0.25))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.5,1)
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  # e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  
  # Compute posterior P(D=1|e)
  
  # ehat = rnorm(n_trials,D * X - bias,sqrt(exp(sigma_m)^2 + exp(sigma_e)^2))
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  # k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  
  p_gauss = pnorm((D * X - bias) / sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))                  #  probability of response on trial t given the citerion
  
  
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  

  # After observing the choice
  a_gauss = rbinom(n_trials, 1, p_gauss)
  
  # Compute evidence conditional on choice (Mills ratio)
  mu_e = D * X - bias
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  z = mu_e / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  action_correction = ifelse(a_gauss == 1,
                  correction_factor * dnorm(z) / pnorm(z),
                  correction_factor * dnorm(z) / pnorm(-z))  # Note: pnorm(-z) here!
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  sigma2 = exp(sigma_m)^2 + exp(sigma_e)^2
  
  marginalization_correction = (1/(sqrt(1 + (Cc^2*(X)^2) / sigma2)))
  
  
  c_mu = ifelse(a_gauss == 1,
                pnorm((((mu_e + action_correction) * (X) * Cc) / sigma2) * marginalization_correction),
                1-pnorm((((mu_e - action_correction) * (X) * Cc) / sigma2) * marginalization_correction))
  
  
  
  
  
  # Dataframe
  df = data.frame(
    p_gauss = p_gauss,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    # C = C,
    a_gauss = a_gauss) %>% 
    mutate(    
      c0 = c0,
      c1 = c1,
      sigma_k = sigma_k,
      sigma_e = sigma_e,
      conf_prec = conf_prec,
      sigma_m = sigma_m,
      bias = bias,
    )
}









simulate_data_p_bias = function(n = 1000, sigma_e,sigma_m,sigma_k,bias, seed = 1997){
  
  set.seed(seed)
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  # sigma_e <- (rnorm(1,-1,0.5))           # noise on evidence (needs to be fixed for now)
  # bias <- rnorm(1,0,0.2)           # noise on evidence (needs to be fixed for now)
  
  # sigma_m = (rnorm(1,-1,0.5))
  conf_prec = (rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  # sigma_k <- (rnorm(1,-2,0.25))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.5,1)
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  e <- rnorm(n_trials, mean = D * X- bias, sd = exp(sigma_e))
  
  
  ehat = rnorm(n_trials,e,(exp(sigma_m)))
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  # k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  
  # p_gauss = pnorm((D * X - bias) / sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))                  #  probability of response on trial t given the citerion
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  # k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  
  p_gauss = pnorm((e) / exp(sigma_k))                  #  probability of response on trial t given the citerion
  
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*ehat / (exp(sigma_e)^2 + exp(sigma_m)^2)),
                1-brms::inv_logit_scaled(2*ehat / (exp(sigma_e)^2 + exp(sigma_m)^2)))
  
  # Dataframe
  df = data.frame(
    p_gauss = p_gauss,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    a_gauss = a_gauss
  ) %>% mutate(   c0 = c0,
                  c1 = c1,
                  sigma_k = sigma_k,
                  sigma_e = sigma_e,
                  conf_prec = conf_prec,
                  sigma_m = sigma_m,
                  bias = bias)
  
  
  
}

simulate_data_margin_bias = function(n = 1000, sigma_e,sigma_m,sigma_k,bias, seed = 1997){
  
  set.seed(seed)
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  # sigma_e <- (rnorm(1,-1,0.5))           # noise on evidence (needs to be fixed for now)
  # bias <- rnorm(1,0,0.2)           # noise on evidence (needs to be fixed for now)
  
  # sigma_m = (rnorm(1,-1,0.5))
  conf_prec = (rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  # sigma_k <- (rnorm(1,-2,0.25))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.5,1)
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  # Simulate evidence
  # e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  
  # Compute posterior P(D=1|e)
  
  # ehat = rnorm(n_trials,D * X - bias,sqrt(exp(sigma_m)^2 + exp(sigma_e)^2))
  
  # Evidence-based stochastic policy (guassian distribution for citera)
  # k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # citerion on trial t
  
  p_gauss = pnorm((D * X - bias) / sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))                  #  probability of response on trial t given the citerion
  
  
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  
  # After observing the choice
  a_gauss = rbinom(n_trials, 1, p_gauss)
  
  # Compute evidence conditional on choice (Mills ratio)
  mu_e = D * X - bias
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  z = mu_e / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  action_correction = ifelse(a_gauss == 1,
                             correction_factor * dnorm(z) / pnorm(z),
                             correction_factor * dnorm(z) / pnorm(-z))  # Note: pnorm(-z) here!
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  sigma2 = exp(sigma_m)^2 + exp(sigma_e)^2
  
  marginalization_correction = (1/(sqrt(1 + (Cc^2) / sigma2)))
  
  
  c_mu = ifelse(a_gauss == 1,
                pnorm(((((mu_e + action_correction)) * Cc) / sigma2) * marginalization_correction),
                1-pnorm(((((mu_e - action_correction)) * Cc) / sigma2) * marginalization_correction))
  
  
  
  
  
  
  # Dataframe
  df = data.frame(
    p_gauss = p_gauss,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    # C = C,
    a_gauss = a_gauss) %>% 
    mutate(    
      c0 = c0,
      c1 = c1,
      sigma_k = sigma_k,
      sigma_e = sigma_e,
      conf_prec = conf_prec,
      sigma_m = sigma_m,
      bias = bias,
    )
}

sim_trials = function(n = 300, seed){
  
  set.seed(seed)
  
  n_trials <- n        # number of trials
  sigma_e <- (rnorm(1,-2,1))           # noise on evidence (needs to be fixed for now)
  beta <- rnorm(1,0,0.5)           # noise on evidence (needs to be fixed for now)
  
  sigma_m = (rnorm(1,-1,1))
  confprec = exp(rnorm(1,2,2))
  meta_bias = rnorm(1,0,0.3)
  lapse = rnorm(1,-4,2)
  
  sigma_k <- (rnorm(1,-2,1))                 # criterion sd on the evidence scale (gaussian distribution)
  
  
  sigmam_beta = rnorm(1,0,0.2)
  meta_bias_beta = rnorm(1,0,0.2)
  
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  c0 = -abs(rnorm(1,5,1))            # Cut points for ordered beta regression, lower bound
  c1 = abs(rnorm(1,5,1))             # Cut points for ordered beta regression, upper bound
  
  # -----------------------------
  # Coherence level (stimulus intensity)
  # -----------------------------
  offsets  = seq(-0.2,0.2,length.out = 10)
  w = c(1, 2, 3, 2, 1, 1, 2, 3, 2, 1)
  
  
  stim_values_per_block =
    rep((brms::inv_logit_scaled(
      (brms::inv_logit_scaled(beta) - 0.5) * 2 +
        rep(offsets / 1/sqrt(exp(sigma_e)^2 + exp(sigma_k)^2), times = w)
    ) - 0.5) * 2,3)
  
  
  X = rep(stim_values_per_block,6)
  
  interval = runif(length(X), 1,5)
  
  # -----------------------------
  # Decision
  # -----------------------------
  
  
  
  sigma1 = exp(sigma_e)^2 + exp(sigma_k)^2
  
  
  mu = (X - (brms::inv_logit_scaled(beta)-0.5)*2)
  
  theta = brms::inv_logit_scaled(lapse)/2 + (1-2*brms::inv_logit_scaled(lapse)/2) * pnorm( mu / sqrt(sigma1))
  
  resp = rbinom(length(X), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(X < 0 & resp == 0, 1,
               ifelse(X > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  
  conf_mu = numeric(length(X))
  
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  
  z = mu / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  mills_pos  = exp(dnorm(z, log=TRUE) - pnorm(z, log.p=TRUE))
  mills_neg  = exp(dnorm(z, log=TRUE) - pnorm(-z, log.p=TRUE))
  
  e_cond = ifelse(resp == 1,
                  mu + correction_factor * mills_pos,
                  mu - correction_factor * mills_neg)
  
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  
  sigma2 = exp(sigma_e)^2 + exp(sigma_m +  sigmam_beta * interval)^2
  
  conf_mu = ifelse(resp == 1,
                   pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2),
                   1-pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2))
  
  conf_mu = ifelse(conf_mu > 0.999, 0.999, ifelse(conf_mu < 0.001,0.001,conf_mu))
  
  conf = numeric(length(conf_mu))
  
  for(i in 1:length(X)){
    conf[i] = rordbeta(n = 1, mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu[i]) + meta_bias + meta_bias_beta * interval[i]), phi = exp(confprec), cutpoints = c(c0, c1))  # Simulated confidence values
  }
  
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = beta,
    sigma_e = sigma_e,
    sigma_k = sigma_k, 
    sigma_m = sigma_m,
    meta_bias = meta_bias,
    confprec = confprec,
    lapse = lapse,
    sigmam_beta = sigmam_beta,
    meta_bias_beta = meta_bias_beta,
    X = X,
    interval = interval,
    ACC = ACC,
    resp = resp,
    c_mu = conf_mu,
    C = conf,
    theta = theta
  ) %>% mutate(trial = 1:n())
  
  
  dataplot = df %>% pivot_longer(cols = c("resp","C")) %>% 
    mutate(ACC = as.factor(ifelse(name == "resp",NA,ACC))) %>% 
    group_by(ACC,name,X) %>% 
    summarize(mean = mean(value),
              se = sd(value) / sqrt(n())) %>% 
    ggplot(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = ACC))+
    geom_pointrange()+
    geom_smooth()+
    facet_wrap(~name,ncol = 1)+
    theme(legend.position = "top")
  
  effectplot = df %>% 
    filter(ACC == 1) %>% 
    mutate(cut_interval = as.factor(cut(interval,3))) %>% 
    group_by(X,cut_interval) %>% 
    summarize(mean = mean(c_mu),
              se = sd(c_mu) / sqrt(n())) %>% 
    ggplot()+
    geom_pointrange(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = (cut_interval)))+
    geom_text(aes(x = 0, y = 0.8, label = paste0("bias = ",round(unique(meta_bias_beta),2), " \n meta_un = ",round(unique(sigmam_beta),2))))+
    geom_smooth(aes(x = X, y = mean, col = (cut_interval)),se = F)+
    theme(legend.position = "top")
  
  plot = dataplot | effectplot 
  
  plot
  
  df$plot <- c(list(plot), rep(list(NA), nrow(df) - 1))
  
  df$id = rnorm(1,0,10000)
  
  return(df)
}


sim_trials_fix = function(n = 300, seed, sigmam_beta,meta_bias_beta){
  
  set.seed(seed)
  
  n_trials <- n        # number of trials
  sigma_e <- -2
  beta <- 0
  
  sigma_m = -1
  confprec = exp(2)
  meta_bias = 0
  lapse = -4  
  sigma_k <- -2
  
  
  c0 = -5
  c1 = 5
  
  
  # -----------------------------
  # Coherence level (stimulus intensity)
  # -----------------------------
  offsets  = seq(-0.2,0.2,length.out = 10)
  w = c(1, 2, 3, 2, 1, 1, 2, 3, 2, 1)
  
  
  stim_values_per_block =
    rep((brms::inv_logit_scaled(
      (brms::inv_logit_scaled(beta) - 0.5) * 2 +
        rep(offsets / 1/sqrt(exp(sigma_e)^2 + exp(sigma_k)^2), times = w)
    ) - 0.5) * 2,3)
  
  
  X = rep(stim_values_per_block,6)
  
  interval = runif(length(X), 1,5)
  
  # -----------------------------
  # Decision
  # -----------------------------
  
  
  
  sigma1 = exp(sigma_e)^2 + exp(sigma_k)^2
  
  
  mu = (X - (brms::inv_logit_scaled(beta)-0.5)*2)
  
  theta = brms::inv_logit_scaled(lapse)/2 + (1-2*brms::inv_logit_scaled(lapse)/2) * pnorm( mu / sqrt(sigma1))
  
  resp = rbinom(length(X), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(X < 0 & resp == 0, 1,
               ifelse(X > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  
  conf_mu = numeric(length(X))
  
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  
  z = mu / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  mills_pos  = exp(dnorm(z, log=TRUE) - pnorm(z, log.p=TRUE))
  mills_neg  = exp(dnorm(z, log=TRUE) - pnorm(-z, log.p=TRUE))
  
  e_cond = ifelse(resp == 1,
                  mu + correction_factor * mills_pos,
                  mu - correction_factor * mills_neg)
  
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  
  sigma2 = exp(sigma_e)^2 + exp(sigma_m +  sigmam_beta * interval)^2
  
  conf_mu = ifelse(resp == 1,
                   pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2),
                   1-pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2))
  
  conf_mu = ifelse(conf_mu > 0.999, 0.999, ifelse(conf_mu < 0.001,0.001,conf_mu))
  
  conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + meta_bias_beta * interval[i])
  
  conf = numeric(length(conf_mu))
  
  # for(i in 1:length(X)){
  # conf[i] = rordbeta(n = 1, mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu[i]) + meta_bias + meta_bias_beta * interval[i]), phi = exp(confprec), cutpoints = c(c0, c1))  # Simulated confidence values
  # }
  
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = beta,
    sigma_e = sigma_e,
    sigma_k = sigma_k, 
    sigma_m = sigma_m,
    meta_bias = meta_bias,
    confprec = confprec,
    lapse = lapse,
    sigmam_beta = sigmam_beta,
    meta_bias_beta = meta_bias_beta,
    X = X,
    interval = interval,
    ACC = ACC,
    resp = resp,
    c_mu = conf_mu,
    C = conf,
    theta = theta
  ) %>% mutate(trial = 1:n())
  
  
  dataplot = df %>% pivot_longer(cols = c("resp","C")) %>% 
    mutate(ACC = as.factor(ifelse(name == "resp",NA,ACC))) %>% 
    group_by(ACC,name,X) %>% 
    summarize(mean = mean(value),
              se = sd(value) / sqrt(n())) %>% 
    ggplot(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = ACC))+
    geom_pointrange()+
    geom_smooth()+
    facet_wrap(~name,ncol = 1)+
    theme(legend.position = "top")
  
  effectplot = df %>% 
    filter(ACC == 1) %>% 
    mutate(cut_interval = as.factor(cut(interval,3))) %>% 
    group_by(X,cut_interval) %>% 
    summarize(mean = mean(c_mu),
              se = sd(c_mu) / sqrt(n())) %>% 
    ggplot()+
    geom_pointrange(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = (cut_interval)))+
    geom_text(aes(x = 0, y = 0.8, label = paste0("bias = ",round(unique(meta_bias_beta),2), " \n meta_un = ",round(unique(sigmam_beta),2))))+
    geom_smooth(aes(x = X, y = mean, col = (cut_interval)),se = F)+
    theme(legend.position = "top")
  
  plot = dataplot | effectplot 
  
  plot
  
  df$plot <- c(list(plot), rep(list(NA), nrow(df) - 1))
  
  df$id = rnorm(1,0,10000)
  
  return(df)
}





library(tidyverse)
library(cowplot)

Cc <- 2 / 1.7   # approximation constant linking inverse logit to probit

lambda <- function(z) dnorm(z) / pnorm(z)   # inverse Mills ratio

rho_fn <- function(sigma_e, sigma_c, sigma_m) {
  sigma_e^2 / (sqrt(sigma_e^2 + sigma_c^2) * sqrt(sigma_e^2 + sigma_m^2))
}

K_fn <- function(sigma_e, sigma_m) {
  1 / sqrt((1 + Cc^2) / (sigma_e^2 + sigma_m^2))
}

I_fn <- function(XD, beta, sigma_e, sigma_m) {
  Cc * (XD - beta) / (sigma_e^2 + sigma_m^2)
}

# S_i for a given action
S_fn <- function(XD, beta, sigma_e, sigma_c, sigma_m, action) {
  rho  <- rho_fn(sigma_e, sigma_c, sigma_m)
  z    <- (XD - beta) / (sqrt(sigma_e^2 + sigma_c^2))
  lam  <- ifelse(action == 1, lambda(z), lambda(-z))
  ifelse(action == 1, (rho * lam * Cc) / sqrt(sigma_e^2 + sigma_m^2), (-rho * lam * Cc) / sqrt(sigma_e^2 + sigma_m^2))
}

# Full confidence mu_c
mu_c_fn <- function(XD, beta, sigma_e, sigma_c, sigma_m, action) {
  K   <- K_fn(sigma_e, sigma_m)
  I_i <- I_fn(XD, beta, sigma_e, sigma_m)
  S_i <- S_fn(XD, beta, sigma_e, sigma_c, sigma_m, action)
  if_else(action == 1,
          pnorm((I_i + (S_i)) * K),   # = pnorm((I_i + S_fn_raw) * K) — see derivation
          1 - pnorm((I_i + (S_i)) * K))
}
