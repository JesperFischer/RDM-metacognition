library(tidyverse)
library(here)
library(cmdstanr)
library(ggplot2)
library(ordbetareg)
library(posterior)
library(furrr)
library(purrr)

###############################################################################
######## Parameters (model) ###################################################
###############################################################################

n = 100  # N simulated values 
x_seq = seq(from = -1, to =1, by = 0.01)    # Sequence of coherence values

beta =  5  # Slope (expressed in precision on normal)
meta_un = -2  # meta uncertainty for correct trials (expressed in precision on normal)
meta_un_IC = 3 # meta uncertainty for incorrect trials (expressed in precision on normal)
alpha = -0.1    # threshold
conf_prec = 10   # confidence precision (precision on beta distribution)

c0 = -5  # Cut points for ordered beta regression
c1 = 5

###############################################################################
######## Simulated values (1 run) #############################################
###############################################################################

D = sample(c(-1,1),n, replace = TRUE)     # True state of the world 
coh = runif(n, min = 0, max = 1)        # coherence
X = D * coh                         # Stimulus 

theta = pnorm(beta * (X-alpha))     # Internal probability of being correct
out = rbinom(n,1,theta)         # simulated response

ACC = ifelse(D == -1 & out == 0,1,ifelse(D == 1 & out == 1,1,0)) # Accuracy based on simulated response

## confidence means for ordered beta reg
conf_mu = numeric(n)
for(i in 1:n){
  if(ACC[i] == 1){
    conf_mu[i] <- pnorm((beta+meta_un) * (abs(X[i] - alpha)))
  } else if(ACC[i] == 0){
    conf_mu[i] <- pnorm(meta_un_IC * (abs(X[i] - alpha)))
  }
}

conf = rordbeta(n=n, mu=conf_mu, phi = conf_prec, cutpoints = c(c0,c1))   # Simulated confidence values

## Data frames
sim_data = data.frame(X, theta, out, ACC, conf_mu, conf)
sim_data$bins = sim_data %>% pull(X) %>% cut(breaks = 20, include.lowest = TRUE)

binned_out = sim_data %>% group_by((bins)) %>%
    summarise(
      bin = mean(X),
      resp = mean(out)
    )

## Plot to check validity of raw simulations

ggplot()+geom_line(data = sim_data, aes(X,theta))+ geom_point(data= binned_out, aes(bin,resp, col ="red"))+theme_minimal()   # Plot to check validity of raw simulations
ggplot(sim_data)+
  geom_line(aes(X,conf_mu, color = factor(ACC, levels = c(0, 1), labels = c("Incorrect", "Correct"))))+
  geom_point(aes(X,conf,color = factor(ACC, levels = c(0, 1), labels = c("Incorrect", "Correct"))))+
  scale_color_manual(values = c("Incorrect" = "salmon", "Correct" = "lightgreen"))+theme_minimal()+theme(legend.position = "none")


#data = data %>% filter(Trialtype=="Main") %>% filter(scale=="conf") %>% filter(resp!="None") %>%filter(RTrating!="None") %>%
  #mutate(up = ifelse(resp == "up", 1, 0)) %>% mutate(e = ifelse(cor_resp == "up", coherence, -coherence)) %>% mutate(SR_conf = as.numeric(as.character(SR_conf)))


###############################################################################
######## Model fit (1 run) ####################################################
###############################################################################

sim_data_list = list(N=n,X = X, a=out, C = conf, ACC = ACC)     # Data list for model fit

model = cmdstan_model(here("Stanmodels", "Heurestic models", "heurestic_no_metabias_Siebe.stan")) # Heurestic model (stan file)

## Model fit 
fit <- model$sample(
  data = sim_data_list,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4,
)


## Check parameter recovery 
fit$summary("alpha")  # threshold
fit$summary("beta")   # Slope (expressed in precision on normal)

fit$summary("meta_un_cor1")   # meta uncertainty for correct trials (expressed in precision on normal)
fit$summary("meta_un_inc1")   # meta uncertainty for incorrect trials (expressed in precision on normal)
fit$summary("meta_un")
fit$summary("meta_un_inc")

fit$summary("conf_prec")    # confidence precision (precision on beta distribution)
fit$summary("c0")         # Cut points for ordered beta regression
fit$summary("c1")

## Extract parameters from fit
alpha_fit = as.numeric(fit$summary("alpha")[,"mean"]) # Extract the threshold fit
beta_fit = as.numeric(fit$summary("beta")[,"mean"])   # Extract the slope fit

meta_un_fit = as.numeric(fit$summary("meta_un_cor1")[,"mean"]) # Extract meta uncertainty for correct trials 
meta_un_IC_fit = as.numeric(fit$summary("meta_un_inc1")[,"mean"]) # Extract meta uncertainty for incorrect trials 

theta_fit = pnorm((x_seq-alpha_fit)*beta_fit) # Estimated best fit 

theta_C_cor = pnorm(abs(x_seq-alpha_fit)*meta_un_fit) # Estimated best fit confidence correct 
theta_C_IC = pnorm(abs(x_seq-alpha_fit)*meta_un_IC_fit) # Estimated best fit confidence incorrect  

## Plot the outcome of the psychometric 
ggplot(data.frame(x_seq,theta_fit))+geom_line(aes(x_seq,theta_fit))+geom_point(data= binned_out, aes(bin,resp, col ="red"))+theme_minimal()
ggplot(data.frame(x_seq,theta_C_cor,theta_C_IC))+geom_line(aes(x_seq,theta_C_cor), col = "green")+geom_line(aes(x_seq,theta_C_IC), col = "red")+
  geom_point(data = sim_data, aes(X,conf,color = factor(ACC, levels = c(0, 1), labels = c("Incorrect", "Correct"))))+
  scale_color_manual(values = c("Incorrect" = "salmon", "Correct" = "lightgreen"))+theme_minimal()+theme(legend.position = "none")


###############################################################################
######## simulate data (extensive) ############################################
###############################################################################

simulate_data <- function(n){
  # -----------------------------
  # Set parameters
  # -----------------------------
  n = n        # number of trials
  
  beta = rnorm(1,1,0.5)       # Slope (expressed in precision on normal, moderate values)
  meta_un_cor = exp(rnorm(1,0.25,0.5))    # Meta uncertainty for correct trials (positive, avoids inversion)
  meta_prec_inc = rnorm(1,0,1)   # Meta uncertainty for incorrect trials (positive, avoids inversion)
  alpha = rnorm(1,0,0.2)        # Threshold
  conf_prec = exp(rnorm(1,2,0.5))    # Confidence precision (precision on beta distribution)
  
  c0 = -abs(rnorm(1,5,1))            # Cut points for ordered beta regression, lower bound
  c1 = abs(rnorm(1,5,1))             # Cut points for ordered beta regression, upper bound
  
  # -----------------------------
  # Coherence level (stimulus intensity)
  # -----------------------------
  coh = runif(n, 0, 1)             # Uniformly sampled stimulus strength
  
  # -----------------------------
  # Stimuli (left or right)
  # -----------------------------
  D = sample(c(-1, 1), n, replace = TRUE)  # True stimulus direction
  
  # -----------------------------
  # Simulate evidence
  # -----------------------------
  X = coh * D                       # Evidence scaled by stimulus
  
  # -----------------------------
  # Decision
  # -----------------------------
  theta = pnorm(exp(beta) * (X - alpha)) # Internal probability of being correct
  resp = rbinom(n, 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(D == -1 & resp == 0, 1,
               ifelse(D == 1 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  conf_mu = numeric(n)
  for(i in 1:n){
    if(ACC[i] == 1){
      conf_mu[i] <- pnorm(exp(beta - meta_un_cor) * abs(X[i] - alpha))       # Confidence for correct trials
    } else {
      conf_mu[i] <- pnorm(meta_prec_inc * abs(X[i] - alpha))    # Confidence for incorrect trials
    }
  }
  
  conf = rordbeta(n = n, mu = conf_mu, phi = conf_prec, cutpoints = c(c0, c1))  # Simulated confidence values
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = exp(beta),
    # exp_beta = exp(beta),
    alpha = alpha, 
    meta_un_cor = meta_un_cor,
    meta_total_cor = exp(beta - meta_un_cor),
    meta_prec_inc = meta_prec_inc,
    conf_prec = conf_prec,
    D = D,
    X = X,
    ACC = ACC,
    resp = resp,
    c_mu = conf_mu,
    C = conf
  )
  
  return(df)
}



# sim_data = simulate_data(100)

#sim_data = df

#ggplot(sim_data)+
  #geom_point(aes(X,resp))
  
#ggplot(sim_data)+
  #geom_point(aes(X,C,color = factor(ACC, levels = c(0, 1), labels = c("Incorrect", "Correct"))))+
  #scale_color_manual(values = c("Incorrect" = "salmon", "Correct" = "lightgreen"))+theme_minimal()+theme(legend.position = "none")



###############################################################################
######## Model fit (extensive) ################################################
###############################################################################
mod <- cmdstan_model(here("Stanmodels", "Heurestic models", "heurestic_no_metabias_Siebe.stan"))

## function for 1 fit
fit = function(nn, mod){
  qq = rnorm(1, 1, 1)
  
  df = simulate_data(n = nn)
  
  datastan = list(a = df$resp,
                  ACC = df$ACC,
                  C = df$C,
                  N = nrow(df),
                  X = df$X)
  
  
  fit = mod$sample(data = datastan,
                   iter_warmup = 500,
                   refresh = 0,
                   adapt_delta = 0.99,
                   iter_sampling = 500,
                   parallel_chains = 4)
  
  div = data.frame(fit$diagnostic_summary())%>% mutate(sim = qq, trial = nn,
                                                       time = fit$time()$total)
  
  
  esti_real = fit$summary(c("lp__","alpha","beta","meta_un_cor","meta_total_cor", "meta_prec_inc", "conf_prec","c0","c1")) %>% 
    mutate(simulated = c(NA,unique(df$alpha),unique(df$beta),unique(df$meta_un_cor),unique(df$meta_total_cor),
                         unique(df$meta_prec_inc),unique(df$conf_prec), unique(df$c0),unique(df$c1))) %>% 
    mutate(div = max(div$num_divergent))
  
  
  result = list(esti_real%>% mutate(sim = qq, trial = nn),
                df%>% mutate(sim = qq, trial = nn),
                div = div %>% mutate(sim = qq, trial = nn))
  
  return(result)
  }


## Model fit 
plan(multisession, workers = 5)  # Make it run on multiple cores (Windows-friendly) 

safe_function <- possibly(fit, otherwise = "Error")   # Prevent crashing when 1 fit is bad

# qq = safe_function(100,mod)

results_list <- future_map(1:200, ~ safe_function(200, mod), .progress = T)   # Run 20 times in parallel

saveRDS(results_list,here::here("Simulations","Heurestic","Parameter_recovery_Hsim_Hfit.RData"))

sum(results_list == "Error")
results_list = results_list[-which(results_list == "Error")]


###############################################################################
######## Check model fit (extensive) ##########################################
###############################################################################

## Check number of runs with divergences
map_dfr(results_list,1) %>% 
  mutate(div = as.factor(ifelse(div > 0,1,0))) %>% select(div,sim) %>% distinct() %>% group_by(div) %>% summarise(n())

## Check parameter recovery 
map_dfr(results_list,1) %>% 
  mutate(div = as.factor(ifelse(div > 0,1,0))) %>% 
  ggplot(aes(x = simulated, y = mean, col = div))+geom_point()+
  geom_abline()+facet_wrap(~variable, scales = "free")

map_dfr(results_list,1) %>% 
  mutate(div = as.factor(ifelse(div > 0,1,0))) %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5, ymax = q95, col = div))+geom_pointrange()+
  geom_abline()+facet_wrap(~variable, scales = "free")

## Check fit time
map_dfr(results_list,3) %>% ggplot(aes(x=time))+geom_histogram(color = "black")+theme_minimal()















