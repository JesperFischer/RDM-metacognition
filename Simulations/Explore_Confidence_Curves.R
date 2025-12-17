# Explore Confidence Curves: Understanding Metacognitive Parameters
# This script generates confidence curves as a function of stimulus strength
# to help understand how alpha, beta, delta, and beta_m affect metacognition

library(tidyverse)

# =============================================================================
# PARAMETER EXPLORATION SETUP
# =============================================================================

# Fixed actor parameters for all scenarios
alpha <- 0.0       # Actor's threshold (keep at 0 for symmetry)
# beta <- -1.0        # Actor's sensitivity


scenarios = expand.grid(
  delta = seq(-0.2, 0, by = 0.2),
  alpha = seq(-0.2, 0, by = 0.2),
  beta_m = seq(0, 5, by = 2.5),
  beta = seq(0, 5, by = 2.5),
  model_type = c("shared", "independent_unbiased", "independent_biased"),
  mixture_weight = 0.5,  # Only used for mixture models
  stringsAsFactors = FALSE
)

# Simulation settings
n_trials_per_point <- 30  # More trials for smooth curves
stimulus_strengths <- seq(0, 0.3, by = 0.1)  # Fine grid
true_states <- c(-1, 1)


source("~/RDM-metacognition/Simulations/Simulate_Metacognition.R")
# =============================================================================
# SIMULATION FUNCTION
# =============================================================================
# MODEL TYPES:
# - "shared": Rater only observes actor's evidence e_i (no independent m_i)
# - "independent_unbiased": Rater observes e_i + independent m_i without bias
# - "independent_biased": Rater observes e_i + independent m_i with bias alpha
# - "gaussian_metanoise": Evidence corrupted by metacognitive noise r_conf ~ N(e_i, beta_m^2)
# - "decay": Signal decay + Gaussian noise r_conf ~ N(e_i * mixture_weight, beta_m^2)
# - "mixture_shared_biased": Mixture of shared and independent_biased models
# - "mixture_shared_unbiased": Mixture of shared and independent_unbiased models
# =============================================================================

# =============================================================================
# RUN SIMULATIONS FOR ALL SCENARIOS
# =============================================================================

set.seed(123)

all_results <- scenarios %>%
  rowwise() %>%
  mutate(
    data = list({
      # Simulate trials for this scenario
      trials <- expand_grid(
        trial = 1:n_trials_per_point,
        X = stimulus_strengths,
        D = sample(true_states, n_trials_per_point, replace = TRUE)
      )
      
      results <- trials %>%
        rowwise() %>%
        mutate(
          sim = list(simulate_trial(D, X, alpha, beta, beta_m, delta, model_type, mixture_weight))
        ) %>%
        unnest(sim) %>%
        ungroup()
    }))
      
# Aggregate by X and accuracy
all_results2 = all_results %>% unnest(data) %>% 
  mutate(X_signed = X*D) %>% 
        group_by(X_signed, correct, beta_m, delta, alpha, beta, model_type) %>%
        summarise(
          mean_confidence = mean(confidence),
          se_confidence = sd(confidence) / sqrt(n()),
          n = n(),
          .groups = "drop"
        ) %>% 
  ungroup()


all_results2 %>% filter(alpha == 0, model_type == "shared") %>% 
  ggplot(aes(x = X_signed,y = mean_confidence, ymin = mean_confidence-2*se_confidence, ymax = mean_confidence+ 2*se_confidence, col = interaction(as.factor(beta),correct)))+
  geom_pointrange()+
  scale_color_manual(values = c("yellow","orange","red","lightgreen","green","darkgreen"))+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")

all_results2 %>% filter(alpha == 0) %>% 
  mutate(mean_correct = correct*n) %>%
  group_by(X_signed, alpha,beta) %>% summarize(mean_correct = sum(mean_correct) / sum(n)) %>% 
  ggplot(aes(x = X_signed,y = mean_correct, col = as.factor(beta)))+
  geom_point()+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "top")


all_results2 %>% filter(alpha == -0.2) %>% 
  ggplot(aes(x = X_signed,y = mean_confidence, ymin = mean_confidence-2*se_confidence, ymax = mean_confidence+ 2*se_confidence, col = interaction(as.factor(beta),correct)))+
  geom_pointrange()+
  scale_color_manual(values = c("yellow","orange","red","lightgreen","green","darkgreen"))+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")

all_results2 %>% filter(alpha == -0.2) %>% 
  mutate(mean_correct = correct*n) %>%
  group_by(X_signed, alpha,beta) %>% summarize(mean_correct = sum(mean_correct) / sum(n)) %>% 
  ggplot(aes(x = X_signed,y = mean_correct, col = as.factor(beta)))+
  geom_point()+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "top")





all_results2 %>% filter(alpha == 0.2) %>% 
  ggplot(aes(x = X_signed,y = mean_confidence, ymin = mean_confidence-2*se_confidence, ymax = mean_confidence+ 2*se_confidence, col = interaction(as.factor(beta),correct)))+
  geom_pointrange()+
  scale_color_manual(values = c("yellow","orange","red","lightgreen","green","darkgreen"))+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")

all_results2 %>% filter(alpha == 0.2) %>% 
  mutate(mean_correct = correct*n) %>%
  group_by(X_signed, alpha,beta) %>% summarize(mean_correct = sum(mean_correct) / sum(n)) %>% 
  ggplot(aes(x = X_signed,y = mean_correct, col = as.factor(beta)))+
  geom_point()+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "top")







all_results %>% filter(beta == 2.5) %>% 
  ggplot(aes(x = X,y = mean_confidence, ymin = mean_confidence-2*se_confidence, ymax = mean_confidence+ 2*se_confidence, col = correct))+
  geom_pointrange()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")

all_results %>% filter(beta == 2.5) %>% mutate(mean_correct = correct*n) %>%
  group_by(X,beta_m,delta) %>% summarize(mean_correct = sum(mean_correct) / sum(n)) %>% 
  ggplot(aes(x = X,y = mean_correct))+
  geom_point()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")



all_results %>% filter(beta == 5) %>% 
  ggplot(aes(x = X,y = mean_confidence, ymin = mean_confidence-2*se_confidence, ymax = mean_confidence+ 2*se_confidence, col = correct))+
  geom_pointrange()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")


all_results %>% filter(beta == 5) %>% mutate(mean_correct = correct*n) %>%
  group_by(X,beta_m,delta) %>% summarize(mean_correct = sum(mean_correct) / sum(n)) %>% 
  ggplot(aes(x = X,y = mean_correct))+
  geom_point()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(beta_m~delta, labeller = label_both)+
  theme(legend.position = "top")

# =============================================================================
# COMPARE MODEL TYPES
# =============================================================================
# Simple comparison showing how different model types affect confidence

all_results2 %>% 
  filter(alpha == 0, beta == 5, beta_m == 5, delta == 0) %>% 
  ggplot(aes(x = X_signed, y = mean_confidence, 
             ymin = mean_confidence - 2*se_confidence, 
             ymax = mean_confidence + 2*se_confidence, 
             col = correct)) +
  geom_pointrange() +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~model_type, labeller = label_both) +
  labs(title = "Model Type Comparison",
       subtitle = "alpha=0, beta=5, beta_m=5, delta=0",
       x = "Signed Stimulus Strength",
       y = "Mean Confidence") +
  theme_bw() +
  theme(legend.position = "top")

