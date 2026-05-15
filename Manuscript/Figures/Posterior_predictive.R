pacman::p_load(cmdstanr, posterior, here,tidyverse,ggplot2)

###################################
# Read RDS
###################################

fit = readRDS(here("fits", "Hierarchical_nolapse_50_clamped_standard.rds"))

fit$summary("gm[7]")

###################################
# Marginal posterior beta
###################################
sigma_m_beta = as_draws_df(fit$draws("gm[7]")) %>% select(`gm[7]`) %>% mutate(prior = rnorm(4000,0,0.5)) %>% pivot_longer()

ggplot(sigma_m_beta)+  geom_histogram(aes(x = `gm[7]`), bins = 70, color = "black",fill = "steelblue")+
  geom_histogram(aes(x = prior), bins = 70, color = "black",fill = "lightblue", alpha =0.5)+
  geom_vline(xintercept = 0, color = 'red', linewidth = 1.2, linetype = "dashed")+
  scale_x_continuous(limits = c(-0.25,0.25))+
  theme_classic()

###################################
# Posterior predictive 
###################################
source(here("Sequential sampling", "utility.R"))