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
true_alpha <- 0.5
true_beta <- 5
true_delta <- -0.1
true_beta_m = 2
true_mix = 0
# Generate data
n_trials <- 2000
X_values <- rnorm(n_trials, true_alpha, 0.5)
D_values <- sample(c(-1, 1), n_trials, replace = TRUE)

sim_data <- map2_dfr(D_values, X_values, ~{
  simulate_trial(D = .x, X = .y, 
                 alpha = true_alpha, 
                 beta = true_beta, 
                 beta_m = true_beta_m,  # Not used in shared model
                 conf_prec = 20,
                 cut0 = -3,
                 cut1 = 3,
                 # mixture_weight = true_mix,
                 delta = true_delta,
                 model_type = "real_independent_biased")
}) %>%
  mutate(D = D_values, X = X_values)

sim_data %>% mutate(SS = D*X) %>% ggplot(aes(x = SS, y = action))+geom_point()
sim_data %>% mutate(SS = D*X) %>% ggplot(aes(x = SS, y = confidence, col = as.factor(correct)))+geom_point()+
  scale_y_continuous(limits = c(0,1))
sim_data %>% mutate(SS = D*X) %>% ggplot(aes(x = SS, y = rt))+geom_point()+geom_smooth()



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
mod = cmdstanr::cmdstan_model(here::here("Stanmodels","independent_biased.stan"))
# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","centered_independent_biased.stan"))

# Compile and fit
fit <- mod$sample(
  data = stan_data,
  chains = 4,
  refresh = 100,
  parallel_chains = 4,
  iter_warmup = 1000,
  adapt_delta = 0.90,
  iter_sampling = 1000
)

fit

mcmc_pairs(fit$draws(c("beta","beta_diff","conf_prec","delta","alpha")),np = nuts_params(fit))



mod = cmdstanr::cmdstan_model(here::here("Stanmodels","full_independent.stan"))
# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","centered_independent_biased.stan"))

# Compile and fit
fit_ind <- mod$sample(
  data = stan_data,
  chains = 4,
  refresh = 100,
  parallel_chains = 4,
  iter_warmup = 1000,
  adapt_delta = 0.90,
  iter_sampling = 1000
)

fit_ind

# mcmc_pairs(fit_ind$draws(c("beta","beta_diff","conf_prec","delta","alpha")),np = nuts_params(fit))


model1 = cmdstanr::cmdstan_model(here::here("Stanmodels","simp_model.stan"))


datastan = list(  N = nrow(sim_data),
                  binom_y = sim_data$D == sim_data$action,
                  Conf = sim_data$confidence,
                  X = sim_data$X * sim_data$D,
                  ACC = sim_data$D == sim_data$action)


fit2 = model1$sample(data = datastan,
                     iter_warmup = 500,
                     iter_sampling = 500,
                     parallel_chains = 4,
                     refresh = 100,
                     adapt_delta = 0.90)

fit2$summary("gm")

# alpha
# beta
# lapse
# confprec
# conf_un
# conf_bias



