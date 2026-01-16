
pacman::p_load(tidyverse,posterior,bayesplot,cmdstanr,patchwork,ordbetareg,shiny,statConfR)
source(here::here("Simulations","MCMC","simulate.R"))
fit = function(nn){
  
  qq = rnorm(1,1,1)
  
  df = simulate_data(n = nn)
  
  mod = cmdstanr::cmdstan_model(here::here("Simulations","MCMC","stanmodels","MCMC_cent_sim.stan"))
  
  datastan = list(a = df$a_gauss,
                  D = df$D,
                  C = df$C,
                  N = nrow(df),
                  X = df$X)
  
  
  fit = mod$sample(data = datastan,
                   iter_warmup = 500,
                   refresh = 0,
                   iter_sampling = 500,
                   parallel_chains = 4)
  
  div = data.frame(fit$diagnostic_summary())%>% mutate(sim = qq, trial = nn,
                                                       time = fit$time()$total)
  
  
  esti_real = fit$summary(c("lp__","sigma_choice","sigma_m","sigma_e","mean_choice","prec_conf","c0","c1")) %>% 
    mutate(simulated = c(NA,unique(df$sigma_k),unique(df$sigma_m),unique(df$sigma_e),
                         unique(df$mean_choice),unique(df$conf_prec),unique(df$c0),unique(df$c1))) %>% 
    mutate(div = max(div$num_divergent))
  
  
  data = df %>% mutate(conf_bin = cut(C,4),
                       stimulus  = as.factor(D),
                       participant = 1,
                       rating = as.numeric(as.factor(conf_bin)),
                       correct = ACC) %>% drop_na()
  
  
  
  mratio = fitMetaDprime(data)
  
  
  result = list(esti_real%>% mutate(sim = qq, trial = nn),
                df%>% mutate(sim = qq, trial = nn),
                div = div,
                mratio = mratio %>% mutate(sim = qq, trial = nn))
  
  return(result)
}


library(furrr)
library(purrr)
library(tidyverse)

plan(multisession, workers = 10)  # Windows-friendly

safe_function <- possibly(fit, otherwise = "Error")

# Run 20 times in parallel
results_list <- future_map(1:2000, ~ safe_function(100), .progress = T)

saveRDS(results_list,here::here("Simulations","MCMC","Saved models","2000sim_1000trials_centered_MCMC.RData"))

results_list = results_list[-which(results_list == "Error")]
