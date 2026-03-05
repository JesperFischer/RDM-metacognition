simulate_data_mcmc = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- (rnorm(1,-2,1))           # noise on evidence (needs to be fixed for now)
  bias <- rnorm(1,0,0.3)           # noise on evidence (needs to be fixed for now)
  
  sigma_m = (rnorm(1,-1,0.5))
  conf_prec = abs(rnorm(1,50,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  sigma_k <- (rnorm(1,-1,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.3,0.3)
  
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
  sigma_e <- (rnorm(1,-2,1))           # noise on evidence (needs to be fixed for now)
  bias <- rnorm(1,0,0.3)           # noise on evidence (needs to be fixed for now)
  
  sigma_m = (rnorm(1,-1,0.5))
  conf_prec = abs(rnorm(1,50,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  sigma_k <- (rnorm(1,-1,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.3,0.3)
  
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

library(future)
library(future.apply)
library(ordbetareg)

pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes,
               brms, patchwork, cowplot,ggpubr,flextable)


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


