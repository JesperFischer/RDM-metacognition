# Simulation functions for heuristic model
library(tidyverse)
library(brms)
library(ordbetareg)
library(truncnorm)

#' Simulate data from the heuristic model with X=1
#' 
#' This function simulates data from the static heuristic model of perceptual
#' decision making and confidence ratings, with stimulus intensity fixed at X=1.
#' 
#' Parameters:
#' - sigma_e: evidence noise (on log scale)
#' - sigma_k: criterion variability (on log scale) 
#' - sigma_m: metacognitive noise (on log scale)
#' - bias: bias in evidence (in units of stimulus)
#' - conf_prec: precision parameter for confidence ratings (phi in ordered beta)
#'
simulate_data_mcmc_x1_custom <- function(
  n = 500,
  sigma_e = -2,
  sigma_k = -1,
  sigma_m = -1,
  bias = 0,
  metabias = 0,
  conf_prec = 50,
  n_coherence_levels = 7,
  seed = NULL
) {
  
  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Convert log-scale parameters to natural scale
  sigma_e_real <- exp(sigma_e)
  sigma_k_real <- exp(sigma_k)
  sigma_m_real <- exp(sigma_m)
  
  # Set random cut points for ordered beta regression
  c0 <- -abs(rnorm(1, 5, 1))
  c1 <- abs(rnorm(1, 5, 1))
  
  # Generate stimulus intensities
  # Create different coherence levels evenly distributed
  X_values <- seq(0.05, 0.95, length.out = n_coherence_levels)
  
  # Repeat each coherence level equally across trials
  trials_per_level <- n %/% n_coherence_levels
  remainder <- n %% n_coherence_levels
  
  X <- c(
    rep(X_values, times = trials_per_level),
    X_values[1:remainder]
  )
  
  # Shuffle the order
  order_idx <- sample(1:length(X))
  X <- X[order_idx]
  
  # Generate true stimulus (left = -1, right = 1)
  D <- sample(c(-1, 1), n, replace = TRUE)
  
  # Generate evidence
  e <- rnorm(n, mean = D * X - bias, sd = sigma_e_real)
  
  # Generate variable criterion
  c <- rnorm(n, 0, sigma_k_real)
  
  # Binary choice based on evidence vs criterion
  a_gauss <- as.numeric(e > c)
  
  # Generate noisy readout of evidence (for metacognitive judgments)
  ehat <- rnorm(n, e, sigma_m_real)
  
  # Compute probability of being correct using logistic function
  # This assumes X=1, and uniform prior on D
  c_mu <- ifelse(
    a_gauss == 1,
    brms::inv_logit_scaled(2 * ehat * 1 / (sigma_e_real^2 + sigma_m_real^2)),
    1 - brms::inv_logit_scaled(2 * ehat * 1 / (sigma_e_real^2 + sigma_m_real^2))
  )
  
  c_mu = brms::inv_logit_scaled(brms::logit_scaled(c_mu) + metabias)
  
  # Clamp confidence to valid range
  # c_mu <- pmax(0.01, pmin(0.99, c_mu))
  

  # Compute accuracy
  ACC <- ifelse(
    a_gauss == 0 & D == -1, 1,
    ifelse(a_gauss == 1 & D == 1, 1, 0)
  )
  
  # Create output dataframe
  df <- data.frame(
    c0 = c0,
    c1 = c1,
    sigma_k = sigma_k,
    sigma_e = sigma_e,
    sigma_m = sigma_m,
    bias = bias,
    sigma_e_real = sigma_e_real,
    sigma_k_real = sigma_k_real,
    sigma_m_real = sigma_m_real,
    slope = sqrt(sigma_e_real^2 + sigma_k_real^2),
    sigma2 = (sigma_e_real^2 + sigma_m_real^2),
    ehat = ehat,
    D = D,
    X = X,
    ACC = ACC,
    c_mu = c_mu,
    a_gauss = a_gauss
  )
  
  return(df)
}

