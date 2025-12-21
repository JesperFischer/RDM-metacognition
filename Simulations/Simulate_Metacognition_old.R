simulate_trial <- function(D, X, alpha, beta, beta_m, delta,conf_prec = 10,cut0 = -4,cut1 = 4,
                           model_type = "independent_biased",
                           mixture_weight = 0.5) {
  # Standardized actor's evidence (unit variance)
  e <- rnorm(1, mean = 0, sd = 1)
  e_tilde <- beta * (D * X - alpha) + e
  
  
  rt = -abs(e_tilde)
  #rlnorm(1,-abs(e_tilde),0.2)+0.1
  
  # Actor's decision
  a_i <- ifelse(e_tilde > 0, 1, -1)
  
  # Correct or incorrect?
  correct <- (a_i == D)
  
  # Asymmetric priors for confidence bias (common to all models)
  p_correct <- 0.5 + delta
  p_incorrect <- 0.5 - delta
  
  # Compute confidence based on model type
  if (model_type == "shared") {
    # Rater only has access to actor's evidence e_i
    mu_e_correct <- beta * (a_i * X - alpha)
    like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
    prob_correct <- p_correct * like_e_correct
    
    mu_e_incorrect <- beta * (-a_i * X - alpha)
    like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
    prob_incorrect <- p_incorrect * like_e_incorrect
    
    
    
  } 

  else if (model_type == "real_independent_biased") {
    # Rater has independent evidence m_i WITH bias
    m = rnorm(1,mean=0, sd=1)
    m_tilde <- beta_m * (D * X - alpha) + m
    
    mu_m_correct <- beta_m * (a_i * X - alpha)
    like_m_correct <- dnorm(m_tilde, mean = mu_m_correct, sd = 1)
    prob_correct <- p_correct * like_m_correct
    
    mu_m_incorrect <- beta_m * (-a_i * X - alpha)
    like_m_incorrect <- dnorm(m_tilde, mean = mu_m_incorrect, sd = 1)
    prob_incorrect <- p_incorrect  * like_m_incorrect
    
  }   
    
  else if (model_type == "independent_unbiased") {
    # Rater has independent evidence m_i WITHOUT bias
    m_tilde <- rnorm(1, mean = beta_m * D * X, sd = 1)
    
    mu_e_correct <- beta * (a_i * X - alpha)
    mu_m_correct <- beta_m * a_i * X
    like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
    like_m_correct <- dnorm(m_tilde, mean = mu_m_correct, sd = 1)
    prob_correct <- p_correct * like_e_correct * like_m_correct
    
    mu_e_incorrect <- beta * (-a_i * X - alpha)
    mu_m_incorrect <- beta_m * (-a_i * X)
    like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
    like_m_incorrect <- dnorm(m_tilde, mean = mu_m_incorrect, sd = 1)
    prob_incorrect <- p_incorrect * like_e_incorrect * like_m_incorrect
    
  } 
  
  else if (model_type == "independent_biased") {
    # Rater has independent evidence m_i WITH bias
    m = rnorm(1,mean=0, sd=1)
    m_tilde <- beta_m * (D * X - alpha) + m
    
    mu_e_correct <- beta * (a_i * X - alpha)
    mu_m_correct <- beta_m * (a_i * X - alpha)
    like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
    like_m_correct <- dnorm(m_tilde, mean = mu_m_correct, sd = 1000)
    prob_correct <- p_correct * like_e_correct * like_m_correct
    
    mu_e_incorrect <- beta * (-a_i * X - alpha)
    mu_m_incorrect <- beta_m * (-a_i * X - alpha)
    like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
    like_m_incorrect <- dnorm(m_tilde, mean = mu_m_incorrect, sd = 1000)
    prob_incorrect <- p_incorrect * like_e_incorrect * like_m_incorrect
    
  } 
  
  else if (model_type == "gaussian_metanoise") {
    # Gaussian meta noise: r_conf ~ N(e_tilde, beta_m^2)
    # beta_m represents metacognitive noise (sigma_meta)
    r_conf <- rnorm(1, mean = e_tilde, sd = 1/beta_m)
    
    # Check for cross-over: if r_conf and e_tilde on opposite sides of decision criterion (0)
    crossover <- sign(e_tilde) != sign(r_conf)
    
    if (crossover) {
      # Cross-over detected: set confidence to maximum certainty
      confidence <- 1.0
      return(data.frame(
        correct = correct,
        confidence = confidence,
        action = a_i,
        model_type = model_type
      ))
    } else {
      # No cross-over: use r_conf for confidence computation
      mu_e_correct <- beta * (a_i * X - alpha)
      like_e_correct <- dnorm(r_conf, mean = mu_e_correct, sd = 1)
      prob_correct <- p_correct * like_e_correct
      
      mu_e_incorrect <- beta * (-a_i * X - alpha)
      like_e_incorrect <- dnorm(r_conf, mean = mu_e_incorrect, sd = 1)
      prob_incorrect <- p_incorrect * like_e_incorrect
    }
    
  } 
  
  else if (model_type == "decay") {
    # Decay model: signal decay + Gaussian metacognitive noise
    # Evidence decays by factor mixture_weight (p_decay), then Gaussian noise added
    # r_conf ~ N(e_tilde * p_decay, beta_m^2)
    decayed_evidence <- e_tilde * mixture_weight
    r_conf <- rnorm(1, mean = decayed_evidence, sd = 1/beta_m)
    
    # Check for cross-over: if r_conf and e_tilde on opposite sides of decision criterion (0)
    crossover <- sign(e_tilde) != sign(r_conf)
    
    if (crossover) {
      # Cross-over detected: set confidence to maximum certainty
      confidence <- 1.0
      return(data.frame(
        correct = correct,
        confidence = confidence,
        action = a_i,
        model_type = model_type
      ))
    } else {
      # No cross-over: use r_conf for confidence computation
      mu_e_correct <- beta * (a_i * X - alpha)
      like_e_correct <- dnorm(r_conf, mean = mu_e_correct, sd = 1)
      prob_correct <- p_correct * like_e_correct
      
      mu_e_incorrect <- beta * (-a_i * X - alpha)
      like_e_incorrect <- dnorm(r_conf, mean = mu_e_incorrect, sd = 1)
      prob_incorrect <- p_incorrect * like_e_incorrect
    }
    
  } 
  
  else if (model_type == "mixture_shared_biased") {
    # Mixture: with prob mixture_weight use shared, otherwise independent_biased
    if (runif(1) < mixture_weight) {
      # Use shared model
      mu_e_correct <- beta * (a_i * X - alpha)
      like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
      prob_correct <- p_correct * like_e_correct
      
      mu_e_incorrect <- beta * (-a_i * X - alpha)
      like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
      prob_incorrect <- p_incorrect * like_e_incorrect
    } else {
      # Use independent_biased model
      m_tilde <- rnorm(1, mean = beta_m * (D * X - alpha), sd = 1)
      
      mu_e_correct <- beta * (a_i * X - alpha)
      mu_m_correct <- beta_m * (a_i * X - alpha)
      like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
      like_m_correct <- dnorm(m_tilde, mean = mu_m_correct, sd = 1)
      prob_correct <- p_correct * like_e_correct * like_m_correct
      
      mu_e_incorrect <- beta * (-a_i * X - alpha)
      mu_m_incorrect <- beta_m * (-a_i * X - alpha)
      like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
      like_m_incorrect <- dnorm(m_tilde, mean = mu_m_incorrect, sd = 1)
      prob_incorrect <- p_incorrect * like_e_incorrect * like_m_incorrect
    }
    
  } 
  
  else if (model_type == "mixture_shared_unbiased") {
    # Mixture: with prob mixture_weight use shared, otherwise independent_unbiased
    if (runif(1) < mixture_weight) {
      # Use shared model
      mu_e_correct <- beta * (a_i * X - alpha)
      like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
      prob_correct <- p_correct * like_e_correct
      
      mu_e_incorrect <- beta * (-a_i * X - alpha)
      like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
      prob_incorrect <- p_incorrect * like_e_incorrect
    } else {
      # Use independent_unbiased model
      m_tilde <- rnorm(1, mean = beta_m * D * X, sd = 1)
      
      mu_e_correct <- beta * (a_i * X - alpha)
      mu_m_correct <- beta_m * a_i * X
      like_e_correct <- dnorm(e_tilde, mean = mu_e_correct, sd = 1)
      like_m_correct <- dnorm(m_tilde, mean = mu_m_correct, sd = 1)
      prob_correct <- p_correct * like_e_correct * like_m_correct
      
      mu_e_incorrect <- beta * (-a_i * X - alpha)
      mu_m_incorrect <- beta_m * (-a_i * X)
      like_e_incorrect <- dnorm(e_tilde, mean = mu_e_incorrect, sd = 1)
      like_m_incorrect <- dnorm(m_tilde, mean = mu_m_incorrect, sd = 1)
      prob_incorrect <- p_incorrect * like_e_incorrect * like_m_incorrect
    }
    
  } 
  
  else {
    stop("Unknown model_type. Use: 'shared', 'independent_unbiased', 'independent_biased', 'gaussian_metanoise', 'decay', 'mixture_shared_biased', or 'mixture_shared_unbiased'")
  }
  
  
  
  confidence <- prob_correct / (prob_correct + prob_incorrect)
  
  # confidence = ordbetareg::rordbeta(n = 1, mu = confidence,conf_prec,cutpoints = c(cut0,cut1))
  # confidence = rnorm(1,confidence,0.05)
  
  
  return(data.frame(
    correct = correct,
    confidence = confidence,
    action = a_i,
    e = e,
    rt = rt,
    # m = m,
    model_type = model_type
  ))
}
