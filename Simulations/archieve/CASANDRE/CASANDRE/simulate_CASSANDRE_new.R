library(tidyverse)
library(here)
library(cmdstanr)
library(ggplot2)
library(ordbetareg)
library(posterior)
library(furrr)
library(purrr)
library(bayesplot)
library(brms)
library(ggtext)
library(patchwork)

###############################################################################
######## Parameters (model) ###################################################
###############################################################################

N = 150    # Simulations
ns = 1        # number of subjects
ncond = 1
ncont = 1 
sampn = 30
sampx = seq(1/(2*sampn), 1 - 1/(2*sampn), length.out = sampn)

sens <- 1.5
crit <- 0.1
meta <- 0.2

nconf <- 0.5
pconf <- 0.5

guess = 0

oris <- runif(N, -1, 1)
DV <- sens * (oris - crit)

meta_un = rlnorm(N, mean = log(1/sens), sd = meta)

resps <- matrix(0, nrow = N, ncol = 4)

for (t in 1:N) {
  
  # CDFs for confidence thresholds
  raw1 <- pnorm(-nconf, mean = DV[t], sd = meta_un[t])
  raw2 <- pnorm(0, mean = DV[t], sd = meta_un[t])
  raw3 <- pnorm(pconf, mean = DV[t], sd = meta_un[t])
  
  # probabilities of each response bin
  p <- numeric(4)
  p[1] <- guess/4 + (1 - guess) * raw1         # high conf, option 1
  p[2] <- guess/4 + (1 - guess) * (raw2 - raw1) # low conf, option 1
  p[3] <- guess/4 + (1 - guess) * (raw3 - raw2) # low conf, option 2
  p[4] <- guess/4 + (1 - guess) * (1 - raw3)    # high conf, option 2
  
  # sample response
  resps[t,] <- rmultinom(1, size = 1, prob = p)
}

# responses array: ncond x ncont x ns x nt x 4
resps_array <- array(0, dim = c(ncond, ncont, ns, N, 4))
resps_array[1,1,1,,] <- resps

# orientations array: ncond x ncont x nt x ns
oris_array <- array(0, dim = c(ncond, ncont, N, ns))
oris_array[1,1,,1] <- oris

stan_data = list(ns=ns,nt=N,ncond=ncond, resps=resps_array, oris=oris_array, sampn=sampn, sampx=sampx, ncont=ncont)


mod = cmdstan_model("C:/Users/User/Downloads/CASANDRE-Perceptual-BHM.stan")


fit = mod$sample(data = stan_data,
                 iter_warmup = 500,
                 iter_sampling = 500,
                 parallel_chains = 4)



###############################################################################
######## Parameters (model 2)  ################################################
###############################################################################

N = 10000    # Simulations
nB = 10 

sens <- 6
sigma_e = 5
crit <- 47
meta <- 0.1
bias = 1.2
conf_prec = 50

c0 = -5
c1 = 5

increments = seq(-3,3,by=0.5)
increments = rep(increments, times = 10)

N = length(increments)

temps = crit + increments*sigma_e

theta = pnorm((temps-crit)/sigma_e)

resp = rbinom(N,1,theta)

meta_un = exp((log(sigma_e)+bias)+meta*rnorm(N))

conf_mu = pnorm(abs(temps-crit)/(meta_un))

conf = rordbeta(n = N, mu = conf_mu, phi = conf_prec, cutpoints = c(c0, c1))  # Simulated confidence values

#conf = rnorm(n = N, mean = conf_mu, sd = conf_prec)  # Simulated confidence values


ggplot()+geom_point(aes(x=temps,y=conf))+geom_line(aes(x=temps,y=theta))


mod = cmdstan_model(here("Stanmodels", "CANSANDRE", "Emperical dataanalysis", "Single Subject", "CASANDRE.stan"))

stan_data = list(N=N, a=resp, X=temps, C=conf)

fit = mod$sample(data = stan_data,
                 iter_warmup = 500,
                 init = 5,
                 iter_sampling = 500,
                 parallel_chains = 4)

variables = c('crit', 'sigma_e',"prec_conf", "mu_sigma_log", "bias")
fit$summary(variables)

mcmc_trace(fit$draws(variables), np = nuts_params(fit))
mcmc_pairs(fit$draws(variables), np = nuts_params(fit))

fit$summary("mu_conf") %>% mutate(sim=conf_mu) %>% ggplot()+geom_point(aes(x=sim,y=mean))+geom_abline()
fit$summary("sigma_dhat") %>% mutate(sim=meta_un) %>% ggplot()+geom_point(aes(x=sim,y=exp(mean)))+
  geom_abline()+scale_x_continuous((limits=c(10,30)))+scale_y_continuous(limits=c(10,30))

meta_un_est = exp(log(sigma_e_est)+sigma_m_est * rnorm(1000))

df <- data.frame(
  x = temps,
  conf = conf
)

ggplot(df, aes(x = x)) +
  geom_point(aes(y = conf), color = "black", alpha = 0.7) 


###############################################################################
######## simulations   ########################################################
###############################################################################

simulate_data = function(){
  # -----------------------------
  # Set parameters
  # -----------------------------
  sigma_e <- exp(rnorm(1,1,0.3)) 
  sigma_m = 0.1
  conf_prec = abs(rnorm(1,100,50))
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  crit = rnorm(1,47,3)
  bias = runif(1, min = 0, max = 1.5)

  increments = seq(-2,2,by=0.5)
  increments = rep(increments, times = 10)
  
  N = length(increments)
  
  temps = crit + increments*sigma_e

  
  theta = pnorm((temps-crit)/sigma_e)
  
  resp = rbinom(N,1,theta)
  
  sigma_e_hat = exp((log(sigma_e)+bias)+sigma_m*rnorm(N))
  
  conf_mu = pnorm(abs(temps-crit)/(sigma_e_hat))
  
  conf = rordbeta(n = N, mu = conf_mu, phi = conf_prec, cutpoints = c(c0, c1))  
  
  # Dataframe
  return(df = data.frame(
    c0 = c0,
    c1 = c1,
    bias = bias,
    sigma_e = sigma_e,
    conf_prec = conf_prec,
    sigma_m = sigma_m,
    sigma_e_hat = sigma_e_hat,
    crit = crit,
    temps = temps,
    c_mu = conf_mu,
    C = conf,
    a = resp)
  )
}  

###############################################################################
######## Model fit   ##########################################################
###############################################################################

mod = cmdstan_model(here("Stanmodels", "CANSANDRE", "Emperical dataanalysis", "Single Subject", "CASANDRE.stan"))


## function for 1 fit
fit = function(mod,N){
  qq = rnorm(1, 1, 1)
  
  df = simulate_data()
  
  datastan = list(a = df$a,
                  C = df$C,
                  N = nrow(df),
                  X = df$temps)
  
  
  fit = mod$sample(data = datastan,
                   iter_warmup = 500,
                   init = 5,
                   iter_sampling = 500,
                   parallel_chains = 4)
  
  div = data.frame(fit$diagnostic_summary())%>% mutate(sim = qq, trial = N,
                                                       time = fit$time()$total)
  
  
  esti_real = fit$summary(c("crit","sigma_e","bias","mu_sigma_log", "prec_conf","c0","c1")) %>% 
    mutate(simulated = c(unique(df$crit),unique(df$sigma_e),unique(df$bias),unique(log(df$sigma_e)+df$bias),
                         unique(df$conf_prec), unique(df$c0),unique(df$c1))) %>% 
    mutate(div = max(div$num_divergent))
  
  draws = fit$draws(variables = c("crit", "sigma_e", "bias"), format = "df")
  
  result = list(esti_real%>% mutate(sim = qq, trial = N),
                df%>% mutate(sim = qq, trial = N),
                div = div %>% mutate(sim = qq, trial = N),
                draws = draws)
  
  return(result)
}


## Model fit 
plan(multisession, workers = 10)  # Make it run on multiple cores (Windows-friendly) 

safe_function <- possibly(fit, otherwise = "Error")   # Prevent crashing when 1 fit is bad

results_list <- future_map(1:100, ~ fit(mod, 90), .progress = T)   # Run 20 times in parallel

saveRDS(results_list,here::here("Simulations","CASANDRE","Parameter_recovery_CASANDRE_PhD.RData"))

sum(results_list == "Error")
results_list = results_list[-which(results_list == "Error")]


map_dfr(results_list,1) %>% 
  mutate(div = as.factor(ifelse(div > 0,1,0))) %>% select(div,sim) %>% distinct() %>% group_by(div) %>% summarise(n())

map_dfr(results_list,1) %>% 
  mutate(div2 = as.factor(ifelse(div > 0,1,0))) %>% filter(div2 == 1) %>% select(div,sim) %>% distinct() %>% ggplot()+ geom_histogram(aes(x=div))

data = 

  
## Plots
PR = map_dfr(results_list,1) %>% filter(div==0 & rhat<1.01) %>% filter(variable %in% c("sigma_e","bias","crit")) %>% 
  filter(!(variable == "sigma_e" & mean>10)) %>% filter(!(variable == "crit" & mean<45 & simulated>47)) %>%
  ggplot(aes(x = simulated, y = mean, ymin=q5, ymax = q95))+
  geom_pointrange(color = "black", alpha = 0.3, size = 0.5)+theme_classic(base_size = 30)+geom_point(size = 3.2)+geom_abline(alpha = 0.5, linetype = "dashed", size = 2, color = "red")+
  facet_wrap(~variable, nrow = 1, scales = "free", labeller = labeller( variable = c( bias = "Bias[hat(sigma)]", crit = "C[d]", sigma_e = "sigma[e]"), ,.default = label_parsed))+ 
  labs(x = "Simulated", y="Estimated", title = "Parameter recovery") + theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
                                                                             axis.title = element_text(size = rel(1.2)),
                                                                             axis.text = element_text(size = rel(1),
                                                                             strip.text = element_text(size = rel(1.2))))

PR1 = map_dfr(results_list,1) %>% filter(div==0 & rhat<1.01) %>% filter(variable == "sigma_e")%>% filter(mean<10) %>% 
      ggplot(aes(x = simulated, y = mean, ymin=q5, ymax = q95))+geom_abline(alpha = 1, linetype = "dashed", size = 2, color = "grey20")+
      geom_pointrange(color = "black", alpha = 0.4, size = 1)+theme_classic(base_size = 30)+ geom_point(size =5)+
        theme(axis.title.x = element_text(angle = 0, hjust = 0.5, size = 70,), 
              plot.title = element_textbox( size = 70, width = unit(1, "npc"),padding = margin(5, 10, 5, 10), margin = margin(10, 0, 10, 0), 
                                            fill = "grey",  linewidth = 1, r = unit(5, "pt"), halign = 0.5))+
      labs(x = NULL, y=NULL, title = " σₑ")+scale_x_continuous(limits = c(1, 6), expand = c(0, 0), breaks = c(1,3.5,6)) +
        scale_y_continuous(limits = c(1,12), expand = c(0, 0), breaks = c(6,12))


PR2 = map_dfr(results_list,1) %>% filter(div==0 & rhat<1.01) %>% filter(variable == "bias") %>% 
      ggplot(aes(x = simulated, y = mean, ymin=q5, ymax = q95))+geom_abline(alpha = 1, linetype = "dashed", size = 2,color = "grey20")+
      geom_pointrange(color = "black", alpha = 0.4, size = 1)+theme_classic(base_size = 30)+ geom_point(size =5)+
      theme(axis.title.x = element_text( angle = 0,hjust = 0.5, size = 70),
            plot.title = element_textbox( size = 60, width = unit(1, "npc"),padding = margin(5,10, 5, 10), margin = margin(10, 0, 10, 0), 
                                          fill = "grey",  linewidth = 1, r = unit(5, "pt"), halign = 0.5))+
      labs(x = NULL, y=NULL, title = "Bias(σ̂)")+scale_x_continuous(limits = c(0, 1.6), expand = c(0, 0), breaks = c(0,0.8,1.6)) +
      scale_y_continuous(limits = c(0,2), expand = c(0, 0), breaks = c(1,2))


PR3 = map_dfr(results_list,1) %>%filter(div==0 & rhat<1.01) %>%  filter(variable == "crit") %>% filter(!(mean<45&simulated>47)) %>% 
      ggplot(aes(x = simulated, y = mean, ymin=q5, ymax = q95))+geom_abline(alpha = 1, linetype = "dashed", size = 2, color = "grey20")+
      geom_pointrange(color = "black", alpha = 0.4, size = 1)+theme_classic(base_size = 30)+ geom_point(size =5)+ 
      theme(axis.title.x = element_text( angle = 0,hjust = 0.5, size = 50), 
            plot.title = element_textbox( size = 60, width = unit(1, "npc"),padding = margin(5,10, 5, 10), margin = margin(10, 0, 10, 0), 
                                          fill = "grey",  linewidth = 1, r = unit(5, "pt"), halign = 0.5))+
      labs(x = NULL, y=NULL, title = "Criterion")+
      scale_x_continuous(limits = c(36, 56), expand = c(0, 0), breaks = c(36,46,56)) +
      scale_y_continuous(limits = c(36,56), expand = c(0, 0), breaks = c(46,56))

sim = 6
results_list[[sim]][[2]] %>% select(temps, a) %>% group_by(temps) %>% summarise(mean = mean(a)) %>% 
  mutate(crit = results_list[[sim]][[1]] %>% filter(variable == "crit") %>% select(mean) %>% as.numeric()) %>% 
  mutate(sigma = results_list[[sim]][[1]] %>% filter(variable == "sigma_e") %>% select(mean) %>% as.numeric()) %>% 
  mutate(est = pnorm((temps-crit)/sigma)) %>% 
  ggplot()+ geom_point(aes(x=temps, y = mean))+geom_line(aes(x=temps,y=est))+theme_classic()

for (i in 1:90){
  sim = i
  p = results_list[[sim]][[2]] %>% select(temps, a) %>% group_by(temps) %>% summarise(mean = mean(a)) %>% 
    mutate(crit = results_list[[sim]][[1]] %>% filter(variable == "crit") %>% select(mean) %>% as.numeric()) %>% 
    mutate(sigma = results_list[[sim]][[1]] %>% filter(variable == "sigma_e") %>% select(mean) %>% as.numeric()) %>% 
    mutate(est = pnorm((temps-crit)/sigma)) %>% ggplot()+ geom_point(aes(x=temps, y = mean))+geom_line(aes(x=temps,y=est))+labs(title = sim)+theme_classic()
  
  print(p)
}


sim_1 = 60


param_df = results_list[[sim_1]][[1]] %>% 
            filter(variable %in% c("crit", "sigma_e", "bias")) %>%
            select(variable, mean)
          
 

posterior_draw = results_list[[sim_1]][[4]]%>% crossing(temps=seq(
  min(unique(results_list[[sim_1]][[2]]$temps)), 
  max(unique(results_list[[sim_1]][[2]]$temps)),
  length.out = 200)) %>% 
  mutate(y_pred = pnorm((temps-crit)/sigma_e)) %>%mutate(c_pred = pnorm(abs(temps-crit)/(bias+sigma_e))) %>%  
  group_by(temps) %>% summarise(ymin = min(y_pred), ymax = max(y_pred), cmin =min(c_pred),cmax = max(c_pred)) %>% 
  mutate(sim = sim_1)



Type_1 = results_list[[sim_1]][[2]] %>% 
          select(temps, a) %>% group_by(temps) %>% summarise(mean = mean(a), sum = sum(a))  %>% mutate(aP = sum+1, bP = 10-sum+1) %>% 
          mutate(lower = qbeta(0.025, aP,bP),upper = qbeta(0.975, aP, bP)) %>% 
          bind_cols(param_df %>% pivot_wider(names_from = variable, values_from = mean)) %>% 
          mutate(est = pnorm((temps-crit)/sigma_e)) %>% 
          ggplot()+ geom_pointrange(aes(x=temps, y = mean, ymin = lower, ymax=upper), color ="#1f77b4", size = 1)+
          geom_smooth(aes(x=temps,y=est), color = "#1f77b4", size = 1.5, alpha = 0.8, se =FALSE, span = 0.5)+
          geom_ribbon(data = posterior_draw, aes(x=temps,ymin=ymin,ymax=ymax), fill = "#1f77b4", alpha=0.1)+
          labs(y = "P(pain)", x= "Temperature (°C)", title = "Binary Response fit")+
          theme_classic(base_size =20)+ 
          theme(plot.title = element_text(margin = margin(b = 30), family = "bold",hjust = 0.5))



df_conf <- bind_rows(results_list[[sim_1]][[2]] %>% mutate(sim = sim_1), results_list[[sim_2]][[2]] %>% mutate(sim = sim_2) ) %>% 
  left_join(param_df, by = "sim") %>% mutate(est = pnorm(abs(temps - crit.y) / (sigma_e.y + bias.y)))


df_mean <- df_conf %>% group_by(sim, temps) %>% summarise(mean_C = mean(C), .groups = "drop")

conf = df_conf %>% 
          ggplot() + 
          geom_point(data = df_mean,aes(x=temps, y= mean_C, color = as.factor(sim)), size = 2.5, alpha = 0.8)+
          geom_smooth(aes(x=temps,y=est, color =as.factor(sim)), size = 0.8, alpha = 1, se = FALSE,span = 0.5)+
          geom_point(aes(x=temps, y = C, color = as.factor(sim)), size =1.5, alpha = 0.5)+
          geom_ribbon(data = posterior, aes(x=temps,ymin=cmin,ymax=cmax, fill = as.factor(sim)), alpha=0.1)+
          scale_color_manual(values = setNames(c("#1f77b4", "#ff7f0e"),c(as.character(sim_1), as.character(sim_2))))+
          scale_fill_manual( values = setNames(c("#1f77b4", "#ff7f0e"), c(as.character(sim_1), as.character(sim_2))), guide = "none")+
          scale_y_continuous(limits=c(0,1))+geom_hline(aes(yintercept=0.5), linetype = "dashed")+ 
          labs(y = "Confidence", x= "Temperature (°C)", title = "Confidence fit")+
          facet_wrap(~ sim, ncol = 1, labeller = labeller( sim = c( "60" = "Participant 60", "93" = "Participant 93")))+
          theme_classic(base_size =20)+ 
          theme(legend.position = "none",
                plot.title = element_text(margin = margin(b = 30), family = "bold",hjust = 0.5))



(Type_1|conf)







