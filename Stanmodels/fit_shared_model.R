# Fit "shared" metacognitive model using Stan
# Estimates alpha (bias), beta (sensitivity), delta (confidence bias)

library(tidyverse)
library(rstan)
library(bayesplot)

# =============================================================================
# GENERATE OR LOAD DATA
# =============================================================================

# Option 1: Generate simulated data
source("~/RDM-metacognition/Simulations/Simulate_Metacognition.R")

# True parameters
true_alpha <- 0.05
true_beta <- 5
true_delta <- 0

# Generate data
n_trials <- 300
X_values <- runif(n_trials, 0, 0.3)
D_values <- sample(c(-1, 1), n_trials, replace = TRUE)

sim_data <- map2_dfr(D_values, X_values, ~{
  simulate_trial(D = .x, X = .y, 
                alpha = true_alpha, 
                beta = true_beta, 
                beta_m = 0,  # Not used in shared model
                conf_prec = 10,
                delta = true_delta,
                model_type = "shared")
}) %>%
  mutate(D = D_values, X = X_values)

sim_data %>% mutate(SS = D*X) %>% ggplot(aes(x = SS, y = action))+geom_point()
sim_data %>% mutate(SS = D*X) %>% ggplot(aes(x = SS, y = confidence, col = as.factor(correct)))+geom_point()+
  scale_y_continuous(limits = c(0,1))

# Prepare data for Stan
stan_data <- list(
  N = nrow(sim_data),
  D = sim_data$D,
  X = sim_data$X,
  a = sim_data$action,
  C = sim_data$confidence
)

# =============================================================================
# FIT MODEL
# =============================================================================
mod = cmdstanr::cmdstan_model(here::here("Stanmodels","shared_model.stan"))

# Compile and fit
fit <- mod$sample(
  data = stan_data,
  chains = 4,
  refresh = 0,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
)

fit
