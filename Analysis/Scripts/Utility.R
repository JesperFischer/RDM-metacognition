#utilt

# Helper functions
psycho = function(x, alpha, beta, lapse) {
  lapse + (1 - 2 * lapse) * (brms::inv_logit_scaled(beta * (x - alpha)))
}

entropy = function(p) {
  -p * log(p) - (1-p) * log(1-p)
}

get_conf = function(ACC, theta, x, alpha) {
  # ACC=1 means correct response
  # If x > alpha and ACC=1 (correct), return theta
  # If x > alpha and ACC=0 (incorrect), return 1-theta
  # If x < alpha and ACC=1 (correct), return 1-theta
  # If x < alpha and ACC=0 (incorrect), return theta
  ifelse(ACC == 1 & x > alpha, theta,
         ifelse(ACC == 1 & x < alpha, 1 - theta,
                ifelse(ACC == 0 & x > alpha, 1 - theta,
                       ifelse(ACC == 0 & x < alpha, theta, 0.5))))
}

get_prob_cor = function(theta, x) {
  # If x > 0, P(correct) = theta
  # If x < 0, P(correct) = 1-theta
  ifelse(x > 0, theta,
         ifelse(x < 0, 1 - theta, 0.5))
}


qordbeta <- function(p, mu, phi, cutzero, cutone) {
  
  # ensure p and mu are same length
  n <- max(length(p), length(mu))
  p  <- rep(p,  length.out = n)
  mu <- rep(mu, length.out = n)
  
  # Beta parameters
  alpha <- mu * phi
  beta  <- (1 - mu) * phi
  
  # mixture weights (logistic cutpoints)
  p0    <- plogis(cutzero)            # mass at 0
  p1    <- 1 - plogis(cutone)         # mass at 1
  p_mid <- 1 - p0 - p1                # mass in (0,1)
  
  # initialize output
  y <- numeric(n)
  
  # regions
  y[p < p0] <- 0
  y[p > (1 - p1)] <- 1
  
  # continuous region: rescale p to (0,1)
  idx <- p >= p0 & p <= (1 - p1)
  p_rescaled <- (p[idx] - p0) / p_mid
  y[idx] <- qbeta(p_rescaled, alpha[idx], beta[idx])
  
  return(y)
}
