# Simulate Agent with Dynamic Prior and Serial Dependence
# Shows how learning rate (eta) affects serial dependence magnitude

library(tidyverse)
library(patchwork)

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# Simulation Parameters
# ============================================================================

simulate_agent <- function(n_trials = 500,
                          eta = 0.1,
                          sigma_e = 1,      # evidence encoding noise
                          sigma_c = 1,      # choice criterion noise
                          sigma_m = 0.5,    # metacognitive noise
                          D_true = 0.5,     # true stimulus distribution: P(D=1)
                          alpha_0 = 1,      # initial pseudocount for D=1
                          beta_0 = 1) {     # initial pseudocount for D=-1
  
  # Storage
  results <- tibble(
    trial = 1:n_trials,
    D_i = NA_real_,           # stimulus
    evidence = NA_real_,      # perceptual evidence
    choice = NA_real_,        # choice: 1 or -1
    confidence = NA_real_,    # confidence rating
    P_D = NA_real_,           # agent's subjective P(D=1)
    alpha = NA_real_,         # pseudocount for D=1
    beta = NA_real_           # pseudocount for D=-1
  )
  
  # Initialize
  alpha <- alpha_0
  beta <- beta_0
  
  for (t in 1:n_trials) {
    # Generate stimulus
    results$D_i[t] <- ifelse(runif(1) < D_true, 1, -1)
    
    # Compute subjective prior
    P_D <- alpha / (alpha + beta)
    results$P_D[t] <- P_D
    results$alpha[t] <- alpha
    results$beta[t] <- beta
    
    # Generate perceptual evidence
    X_i <- results$D_i[t]
    results$evidence[t] <- rnorm(1, 
                                 mean = X_i * results$D_i[t], 
                                 sd = sigma_e)
    
    # Generate choice using dynamic prior
    # P(a=1 | D, X, P(D)) = Phi((D*X - (sigma_e^2/2*X)*log(P(D)/(1-P(D)))) / sqrt(sigma_e^2 + sigma_c^2))
    
    log_odds <- log(P_D / (1 - P_D))
    bias_term <- (sigma_e^2 / (2 * X_i)) * log_odds
    threshold <- -bias_term  # agent chooses 1 if evidence > threshold
    
    results$choice[t] <- ifelse(results$evidence[t] > threshold, 1, -1)
    
    # Generate confidence from readout of evidence using inverse logit
    # Matches part2-current-model.qmd metacognition model
    # r_i ~ N(e_i, sigma_m): noisy readout of evidence
    r_i <- rnorm(1, mean = results$evidence[t], sd = sigma_m)
    
    # Confidence = P(correct | r_i, a, X_i)
    # For a=1: g(2*r_i*X_i / (sigma_e^2 + sigma_m^2))
    # For a=-1: 1 - g(2*r_i*X_i / (sigma_e^2 + sigma_m^2))
    # where g(x) = 1/(1+exp(-x)) is inverse logit
    
    logit_term <- 2 * r_i * X_i / (sigma_e^2 + sigma_m^2)
    
    if (results$choice[t] == 1) {
      results$confidence[t] <- plogis(logit_term)  # inverse logit
    } else {
      results$confidence[t] <- 1 - plogis(logit_term)
    }
    
    # Update beliefs for next trial
    # (If this is not the last trial)
    if (t < n_trials) {
      c_t <- results$confidence[t]  # confidence just generated
      a_t <- results$choice[t]       # choice just made
      
      if (a_t == 1) {
        alpha <- alpha + eta * c_t
        beta <- beta + eta * (1 - c_t)
      } else {
        alpha <- alpha + eta * (1 - c_t)
        beta <- beta + eta * c_t
      }
    }
  }
  
  return(results)
}

# ============================================================================
# Run simulations at different learning rates
# ============================================================================

eta_values <- c(0.01, 0.05, 0.1, 0.2, 0.5)

all_simulations <- map_df(eta_values, ~{
  simulate_agent(n_trials = 500, eta = .x) %>%
    mutate(eta = .x)
})

# ============================================================================
# Compute Serial Dependence
# ============================================================================

serial_dep <- all_simulations %>%
  group_by(eta) %>%
  mutate(
    prev_choice = lag(choice),
    trial_num = row_number()
  ) %>%
  filter(!is.na(prev_choice), trial_num > 2) %>%  # Skip first trial
  summarise(
    # P(a_i = 1 | a_{i-1} = 1)
    p_repeat = mean(choice[prev_choice == 1] == 1),
    # P(a_i = 1 | a_{i-1} = -1)
    p_switch = mean(choice[prev_choice == -1] == 1),
    # Serial dependence = P(repeat) - P(switch)
    serial_dep_magnitude = mean(choice[prev_choice == 1] == 1) - 
                           mean(choice[prev_choice == -1] == 1),
    .groups = 'drop'
  )

print("Serial Dependence by Learning Rate:")
print(serial_dep)

# ============================================================================
# Visualization 1: Serial Dependence vs Learning Rate
# ============================================================================

p1 <- serial_dep %>%
  ggplot(aes(x = eta, y = serial_dep_magnitude)) +
  geom_point(size = 4, color = "steelblue") +
  geom_line(color = "steelblue", linetype = "dashed") +
  labs(
    x = "Learning Rate (η)",
    y = "Serial Dependence\nP(repeat) - P(switch)",
    title = "Serial Dependence Increases with Learning Rate"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# ============================================================================
# Visualization 2: Trial-by-trial Prior Evolution
# ============================================================================

p2 <- all_simulations %>%
  filter(eta %in% c(0.01, 0.1, 0.5)) %>%
  mutate(eta = factor(eta, levels = c(0.01, 0.1, 0.5),
                      labels = c("η=0.01 (slow learner)",
                                "η=0.1 (moderate)",
                                "η=0.5 (fast learner)"))) %>%
  ggplot(aes(x = trial, y = P_D, color = eta)) +
  geom_line(alpha = 0.7, size = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
  facet_wrap(~eta) +
  labs(
    x = "Trial",
    y = "P(D = 1)",
    title = "Evolution of Agent's Belief About Stimulus Distribution",
    color = "Learning Rate"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# ============================================================================
# Visualization 3: Choice Autocorrelation
# ============================================================================

compute_autocorr <- function(data, max_lag = 20) {
  lags <- 1:max_lag
  acf_vals <- map_dbl(lags, ~{
    cor(data$choice[1:(nrow(data) - .)], 
        data$choice[(. + 1):nrow(data)],
        use = "complete.obs")
  })
  
  tibble(lag = lags, acf = acf_vals)
}

p3 <- all_simulations %>%
  group_by(eta) %>%
  nest() %>%
  mutate(acf_data = map(data, compute_autocorr)) %>%
  unnest(acf_data) %>%
  select(eta, lag, acf) %>%
  ggplot(aes(x = lag, y = acf, color = factor(eta))) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
  labs(
    x = "Lag (trials)",
    y = "Choice Autocorrelation",
    title = "Choice Autocorrelation Increases with Learning Rate",
    color = "η"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# ============================================================================
# Combine and save
# ============================================================================

combined <- p1 / p2 / p3 + 
  plot_layout(heights = c(1, 2, 1.5))

ggsave("serial_dependence_simulation.png", combined, width = 12, height = 10, dpi = 300)
print(combined)

# ============================================================================
# Additional Analysis: Posterior Stability
# ============================================================================

posterior_stats <- all_simulations %>%
  group_by(eta) %>%
  summarise(
    P_D_mean = mean(P_D, na.rm = TRUE),
    P_D_sd = sd(P_D, na.rm = TRUE),
    P_D_min = min(P_D, na.rm = TRUE),
    P_D_max = max(P_D, na.rm = TRUE),
    .groups = 'drop'
  )

print("\nPosterior Belief Statistics:")
print(posterior_stats)

cat("\nInterpretation:\n")
cat("- Low η: Agent is conservative, priors stay close to 0.5 (weakly influenced by trials)\n")
cat("- High η: Agent is reactive, priors fluctuate more (quickly adapts to recent evidence)\n")
cat("- Higher η → stronger serial dependence because beliefs drift with recent choices\n")
