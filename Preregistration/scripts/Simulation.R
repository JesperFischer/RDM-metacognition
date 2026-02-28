
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
  
  e = qnorm(p_gauss) * (exp(sigma_k))
  

  # After observing the choice
  a_gauss = rbinom(n_trials, 1, p_gauss)
  
  # Compute evidence conditional on choice (Mills ratio)
  mu_e = D * X - bias
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  z = mu_e / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  e_cond = ifelse(a_gauss == 1,
                  mu_e + correction_factor * dnorm(z) / pnorm(z),
                  mu_e - correction_factor * dnorm(z) / pnorm(-z))  # Note: pnorm(-z) here!
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  sigma2 = exp(sigma_m)^2 + exp(sigma_e)^2
  
  c_mu = ifelse(a_gauss == 1,
                pnorm((1/(sqrt(1 + (Cc^2*X^2) / sigma2))) * (e_cond * X * Cc) / sigma2),
                1-pnorm((1/(sqrt(1 + (Cc^2*X^2) / sigma2))) * (e_cond * X * Cc) / sigma2))

  
  
  
  
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

