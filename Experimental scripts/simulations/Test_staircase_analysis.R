pacman::p_load("tidyverse", "here", "cmdstanr", "ggplot2", 'jsonlite', "purrr", "LaplacesDemon")

## Functions 
weibull = function(x, alpha, beta, gamma, lapse){
  gamma + (1-lapse-gamma)*(1-exp(-(x/alpha)^beta))
}

# load staircase
sim = fromJSON(here("Experimental scripts", "simulations", "sim_data.json"))
sim_combined = map_dfr(sim, ~.x, .id = "sim") %>% group_by(sim) %>% mutate(idx = row_number())

## Plot staircase
ggplot(sim_combined)+geom_point(aes(x = idx, y = C, color = as.factor(out))) + geom_hline(aes(yintercept = alpha))+facet_wrap(~sim, scales = "free", nrow = 2)+theme_classic()

## Fit 
mod = cmdstan_model(here("Experimental scripts", "simulations", "fit_psychometric_weibul.stan"))
results_list = list()

for (i in 1:8){
  current_df = sim_combined %>% filter(sim == i-1)
  data_stan = list(gamma = 0.5, N = nrow(current_df), x = current_df$C, y = current_df$out) 
  
  fit = mod$sample(data = data_stan, chains = 4, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 4)
  
  current_df$lapse_stan = as.numeric(fit$summary("lapse")$mean)
  current_df$thresh_stan = as.numeric(fit$summary("alpha")$mean)
  current_df$slope_stan = as.numeric(fit$summary("beta")$mean)
  
  results_list[[i]] <- current_df
}

sim_combined <- bind_rows(results_list)

psychometric = sim_combined %>% select(alpha, beta, lapse, slope_est, thresh_est, slope_stan, thresh_stan, lapse_stan, sim) %>% 
  distinct() %>% slice(rep(1:n(), each = 100)) %>% group_by(sim) %>% mutate(dotlife = seq(from = 1, to = 30, length.out = 100)) %>% rowwise() %>% 
  mutate(real = weibull(dotlife, alpha, beta, 0.5, lapse)) %>% mutate(SC = weibull(dotlife, thresh_est, slope_est, 0.5,0.03)) %>% 
  mutate(stan = weibull(dotlife, thresh_stan, slope_stan, 0.5, lapse_stan))

## plot fits 

ggplot(psychometric)+geom_line(aes(x= dotlife, y = real), color = "red")+geom_line(aes(x= dotlife, y = SC), linetype = "dashed")+
  geom_line(aes(x= dotlife, y = stan))+facet_wrap(~sim, scales = "free")+theme_classic()

psychometric %>% select(-c(dotlife, stan, SC, real)) %>% distinct() %>% ggplot()+geom_point(aes(x=sim, y = alpha), color = "red", size =3)+ geom_point(aes(x=sim, y = thresh_est), shape = 15)+
  geom_point(aes(x=sim, y = thresh_stan))+labs(title = "bias")+theme_classic()



psychometric %>% select(-c(dotlife, stan, SC, real)) %>% distinct() %>% ggplot()+geom_point(aes(x=sim, y = beta), color = "red", size =3)+ geom_point(aes(x=sim, y = slope_est), shape = 15)+
  geom_point(aes(x=sim, y = slope_stan))+labs(title = "slope")+theme_classic()


