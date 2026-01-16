pacman::p_load(tidyverse,posterior,bayesplot,cmdstanr,patchwork,ordbetareg,shiny,statConfR)

simulate_data = function(n = 1000){
  
  
  
  # -----------------------------
  # Set parameters
  # -----------------------------
  n_trials <- n        # number of trials
  v <- abs(rnorm(1,0,0.5))           # noise on evidence (needs to be fixed for now)
  v_beta <- abs(rnorm(1,2,0.5))           # noise on evidence (needs to be fixed for now)
  
  k = abs(rnorm(1,1,0.5))
  a <- abs(rnorm(1,4,2))                 # criterion sd on the evidence scale (gaussian distribution)
  
  sigma_e <- abs(rnorm(1,1,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  sigma_e = 1
  
  sigma_m <- abs(rnorm(1,1,0.5))                 # criterion sd on the evidence scale (gaussian distribution)
  
  
  conf_prec = abs(rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  # mean_choice <- rnorm(1,0,0.05)
  mu_k = 0
  
  # Coherence level (stimulus intensity)
  X <- abs(rnorm(n_trials,0,1))
  
  # stimuli (left or right)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus
  
  vs <- v + v_beta * X * D   # length n
  
  res <- lapply(
    vs,
    function(v_i) simulate_ddm_collapsing(
      v = v_i,
      a = a,
      k = k,
      s = sigma_e
    )
  )
  df <- do.call(rbind, lapply(res, as.data.frame))
  
  mu_k = 0
  sigma_k <- 0   # motor / readout noise
  
  df$ehat <- df$e + rnorm(n_trials, 0, sigma_k)
  

  p_gauss = pnorm((df$e-mu_k)/sigma_k)                  #  probability of response on trial t given the citerion
  a_gauss = rbinom(n_trials,1,p_gauss)               # binary response
  
  c_mu = ifelse(a_gauss == 1, brms::inv_logit_scaled(2*df$ehat*X / (sigma_e^2 + sigma_m^2)),
                1-brms::inv_logit_scaled(2*df$ehat*X / (sigma_e^2 + sigma_m^2)))
  
  
  
  C = ordbetareg::rordbeta(n_trials,c_mu,conf_prec,c(c0,c1))
  
  # Dataframe
  df = data.frame(
    v = v,
    v_beta = v_beta,
    k = k,
    a = a,
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    # mean_choice = mean_choice,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    e = df$e,
    ehat = df$ehat,
    RT = df$rt,
    D = D,
    X =X,
    ACC = ifelse(a_gauss == 0 & D == -1, 1, ifelse(a_gauss == 1 & D == 1, 1,0)),
    c_mu = c_mu,
    C = C,
    a_gauss = a_gauss
  )
  
  return(df)
}  



simulate_ddm_collapsing <- function(
    v = 0.3,
    a = 1,
    k = 0.5,
    z = 0,
    s = 1,
    dt = 0.001,
    t_max = 5
) {
  
  t <- seq(0, t_max, by = dt)
  x <- z
  
  boundary <- function(t) return(a * exp(-k * t))
  
  for (i in 2:length(t)) {
    x <- x + v * dt + s * sqrt(dt) * rnorm(1)
    b <- boundary(t[i])
    
    if (x >= b) {
      return(list(
        rt = t[i],
        choice = 1,
        x_hit = x,
        e = b
      ))
    }
    
    if (x <= -b) {
      return(list(
        rt = t[i],
        choice = 0,
        x_hit = x,
        e = -b
      ))
    }
  }
  
  list(rt = NA, choice = NA, x_hit = NA, e = NA)
}
