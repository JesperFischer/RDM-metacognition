# Simulation function for the TRUE model (Equation 2) with omega weighting parameter
# Uses the inverse logit function directly (NOT the Phi approximation)
# This simulates confidence judgements with a weight on the inverse Mills ratio adjustment

#' Simulate confidence data with omega parameter (TRUE MODEL)
#' 
#' @param n_trials Number of trials to simulate
#' @param sigma_e Standard deviation of evidence encoding noise
#' @param sigma_c Standard deviation of criterion noise (choice variability)
#' @param sigma_m Standard deviation of metacognitive noise
#' @param omega Weight on the inverse Mills ratio adjustment (1 = normative, >1 = overweight choice, <1 = underweight choice, 0 = ignore choice)
#' @param beta Evidence bias parameter (centered at 0.5 on probability scale)
#' @param X_values Vector of stimulus strengths (unsigned)
#' @param lapse Lapse rate for choices (probability of random choice)
#' @param conf_prec Precision parameter for beta regression of confidence ratings


#' @return A data frame with columns: D (true state), X (stimulus strength), a (choice), 
#'         e (true evidence), r (metacognitive evidence), C (confidence rating), 
#'         C_mu (expected confidence), ACC (accuracy)

simulate_with_omega <- function(n_trials = 1000,
                                 sigma_e = 1.0,
                                 sigma_c = 0.5,
                                 sigma_m = 1.0,
                                 omega = 1.0,
                                 beta = 0.0,
                                 X_values = c(0.01, 0.5, 1, 2),
                                 lapse = 0.01,
                                 conf_prec = 100) {
  
  # Helper function for inverse Mills ratio
  lambda <- function(z) {
    # For numerical stability
    if (z > -5) {
      return(dnorm(z) / pnorm(z))
    } else {
      # Asymptotic expansion for z << 0
      return(-1 / z)
    }
  }
  
  # Vectorized lambda
  lambda_vec <- Vectorize(lambda)
  
  # Setup
  X <- sample(X_values, n_trials, replace = TRUE)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)
  
  # Generate evidence
  e <- rnorm(n_trials, mean = X * D - beta, sd = sigma_e)
  
  # Generate criterion noise
  c <- rnorm(n_trials, mean = 0, sd = sigma_c)
  
  # Generate choices (with lapse)
  p_choice <- lapse/2 + (1 - 2*lapse/2) * pnorm((e - c) / sigma_c)
  a <- rbinom(n_trials, 1, pnorm(e / sigma_c))
  
  # Calculate accuracy
  ACC <- ifelse(a == 1 & D == 1, 1, 
                ifelse(a == 0 & D == -1, 1, 0))
  
  C = 2/1.7
  
  # Compute expected metacognitive evidence with omega weighting
  # Using the TRUE model with inverse logit (not Phi approximation)
  mu <- X * D - beta
  sigma1 <- sigma_e^2 + sigma_c^2
  sigma2 <- sigma_e^2 + sigma_m^2
  
  # Conditional expectation with omega weight on inverse Mills ratio
  z_z <- mu / sqrt(sigma1)
  
  mills_adjustment <- omega * (sigma_e^2 / sqrt(sigma1)) * lambda_vec(z_z)
  
  mills_adjustment_neg <- omega * (sigma_e^2 / sqrt(sigma1)) * lambda_vec(-z_z)
  
  # Expected metacognitive evidence conditional on choice
  E_r_given_a1 <- mu + mills_adjustment
  E_r_given_a0 <- mu - mills_adjustment_neg
  
  correction = 1/(sqrt(1 + ((C^2 * X^2) / (sigma2))))
  
  # Expected confidence using TRUE inverse logit function (Equation 2)
  # P(correct | r, a, X) = g(2 * r * X / (sigma_e^2 + sigma_m^2))
  # where g(x) = 1/(1+exp(-x)) is the inverse logit
  C_mu <- numeric(n_trials)
  for (i in 1:n_trials) {
    if (a[i] == 1) {
      # Use expected r given a=1
      arg <-  ((E_r_given_a1[i] * X[i] * C) / sigma2)* correction
      C_mu[i] <- plogis(arg)  # plogis is R's inverse logit function
    } else {
      # Use expected r given a=0 (which is a=-1 in the model)
      arg <- ((E_r_given_a0[i] * X[i] * C) / sigma2) * correction
      C_mu[i] <- 1 - plogis(arg)
    }
  }
  
  # Ensure C_mu is in valid range
  C_mu <- pmax(0.001, pmin(0.999, C_mu))
  
  # Generate confidence ratings using ordered beta regression
  # (simplified version - just using beta regression on (0,1))
  C <- rbeta(n_trials, 
             shape1 = C_mu * conf_prec,
             shape2 = (1 - C_mu) * conf_prec)
  
  # Create output data frame
  df <- data.frame(
    trial = 1:n_trials,
    D = D,
    X = X,
    XD = X * D,
    e = e,
    a = a,
    ACC = ACC,
    C_mu = C_mu,
    C = C,
    omega = omega,
    sigma_e = sigma_e,
    sigma_c = sigma_c,
    sigma_m = sigma_m
  )
  
  return(df)
}


plot_confidence_curves <- function(df, title = NULL) {
  
  if (!require("tidyverse")) {
    stop("tidyverse package required")
  }
  
  if (is.null(title)) {
    title <- sprintf("ω = %.2f, σ_e = %.2f, σ_c = %.2f, σ_m = %.2f",
                    unique(df$omega), unique(df$sigma_e), 
                    unique(df$sigma_c), unique(df$sigma_m))
  }
  
  p <- df %>%
    mutate(SS = X * D,
           Accuracy = factor(ACC, levels = c(0, 1), labels = c("Error", "Correct"))) %>%
    group_by(SS, Accuracy) %>%
    summarize(
      mean_conf = mean(C_mu),
      se_conf = sd(C_mu) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = SS, y = mean_conf, color = Accuracy, group = Accuracy)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_conf - 2*se_conf, 
                      ymax = mean_conf + 2*se_conf), 
                  width = 0.1, alpha = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "Signed Stimulus Strength (X·D)",
      y = "Expected Confidence",
      color = "Accuracy"
    ) +
    scale_color_manual(values = c("Error" = "#D55E00", "Correct" = "#009E73")) +
    ylim(0, 1) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 12)
    )
  
  return(p)
}


plot_confidence_curves_action <- function(df, title = NULL) {
  
  if (!require("tidyverse")) {
    stop("tidyverse package required")
  }
  
  if (is.null(title)) {
    title <- sprintf("ω = %.2f, σ_e = %.2f, σ_c = %.2f, σ_m = %.2f",
                     unique(df$omega), unique(df$sigma_e), 
                     unique(df$sigma_c), unique(df$sigma_m))
  }
  
  p <- df %>%
    mutate(SS = X * D,
           Accuracy = factor(ACC, levels = c(0, 1), labels = c("Error", "Correct"))) %>%
    mutate(a = as.factor(a)) %>% 
    group_by(SS, a) %>%
    summarize(
      mean_conf = mean(C_mu),
      se_conf = sd(C_mu) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = SS, y = mean_conf, color = a, group = a)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_conf - 2*se_conf, 
                      ymax = mean_conf + 2*se_conf), 
                  width = 0.1, alpha = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "Signed Stimulus Strength (X·D)",
      y = "Expected Confidence",
      color = "a"
    ) +
    scale_color_manual(values = c("#D55E00", "#009E73")) +
    ylim(0, 1) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 12)
    )
  
  return(p)
}



plot_confidence_curves(simulate_with_omega(sigma_e = exp(-0.5), sigma_c = 0, sigma_m = 0, omega = 1))

plot_confidence_curves_action(simulate_with_omega(sigma_e = exp(-0.5), sigma_c = 0, sigma_m = 0, omega = 1.5))
