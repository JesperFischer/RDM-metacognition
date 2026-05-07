
simulate_data_heu_x1 = function(n = 1000){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  sigma_e <- (rnorm(1,-1.5,0.5))           # noise on evidence (needs to be fixed for now)
  bias <- rnorm(1,0,0.3)           # noise on evidence (needs to be fixed for now)
  bias = 0
  sigma_m = (rnorm(1,-1,0.5))
  conf_prec = rnorm(1,3,1)
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
  sigma_k <- (rnorm(1,-1,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  # sigma_k = -10
  # Coherence level (stimulus intensity)
  X <- rtruncnorm(n_trials,0,1,0.5,0.3)
  # X <- 0.3
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  XD = X*D
  # Simulate evidence
  # e <- rnorm(n_trials, mean = D * X - bias, sd = exp(sigma_e))
  
  # Compute posterior P(D=1|e)
  
  
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  
  mu =  D * X - bias
  
  theta = pnorm( mu / sigma_total)
  
  resp = rbinom(length(D), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(XD < 0 & resp == 0, 1,
               ifelse(XD > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  
  conf_mu = numeric(length(X))
  
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
  
  sigma2 = exp(sigma_e)^2 + exp(sigma_m)^2
  
  conf_mu = ifelse(resp == 1,
                   pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2),
                   1-pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2))
  
  conf_mu = ifelse(conf_mu > 0.999, 0.999, ifelse(conf_mu < 0.001,0.001,conf_mu))
  
  C = numeric(length(conf_mu))
  
  for(i in 1:length(D)){
    C[i] = rordbeta(n = 1, mu = conf_mu[i], phi = exp(conf_prec), cutpoints = c(c0, c1))  # Simulated confidence values
  }
  
  
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
    D = D,
    X =X,
    ACC = ACC,
    c_mu = conf_mu,
    C = C,
    a_gauss = resp
  )
}