
pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes,
               brms, patchwork, cowplot,ggpubr,flextable)


sim_trials = function(n = 300, seed){
  
  set.seed(seed)
  
  n_trials <- n        # number of trials
  sigma_e <- (rnorm(1,-2,1))           # noise on evidence (needs to be fixed for now)
  beta <- rnorm(1,0,0.5)           # noise on evidence (needs to be fixed for now)
  
  sigma_m = (rnorm(1,-1,1))
  confprec = exp(rnorm(1,2,2))
  meta_bias = rnorm(1,0,0.3)
  lapse = rnorm(1,-4,2)
  
  sigma_k <- (rnorm(1,-2,1))                 # criterion sd on the evidence scale (gaussian distribution)
  
  
  sigmam_beta = rnorm(1,0,0.2)
  meta_bias_beta = rnorm(1,0,0.2)
  
  c0 = -abs(rnorm(1,5,1))
  c1 = abs(rnorm(1,5,1))
  
  
  c0 = -abs(rnorm(1,5,1))            # Cut points for ordered beta regression, lower bound
  c1 = abs(rnorm(1,5,1))             # Cut points for ordered beta regression, upper bound
  
  # -----------------------------
  # Coherence level (stimulus intensity)
  # -----------------------------
  offsets  = seq(-0.2,0.2,length.out = 10)
  w = c(1, 2, 3, 2, 1, 1, 2, 3, 2, 1)
  
  
  stim_values_per_block =
    rep((brms::inv_logit_scaled(
      (brms::inv_logit_scaled(beta) - 0.5) * 2 +
        rep(offsets / 1/sqrt(exp(sigma_e)^2 + exp(sigma_k)^2), times = w)
    ) - 0.5) * 2,3)
  
  
  X = rep(stim_values_per_block,6)
  
  interval = runif(length(X), 1,5)
  
  # -----------------------------
  # Decision
  # -----------------------------
  
  
  
  sigma1 = exp(sigma_e)^2 + exp(sigma_k)^2
  
  
  mu = (X - (brms::inv_logit_scaled(beta)-0.5)*2)
  
  theta = brms::inv_logit_scaled(lapse)/2 + (1-2*brms::inv_logit_scaled(lapse)/2) * pnorm( mu / sqrt(sigma1))
  
  resp = rbinom(length(X), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(X < 0 & resp == 0, 1,
               ifelse(X > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  
  conf_mu = numeric(length(X))
  
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  
  z = mu / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  mills_pos  = exp(dnorm(z, log=TRUE) - pnorm(z, log.p=TRUE))
  mills_neg  = exp(dnorm(z, log=TRUE) - pnorm(-z, log.p=TRUE))
  
  e_cond = ifelse(resp == 1,
                  mu + correction_factor * mills_pos,
                  mu - correction_factor * mills_neg)
  
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  
  sigma2 = exp(sigma_e)^2 + exp(sigma_m +  sigmam_beta * interval)^2
  
  conf_mu = ifelse(resp == 1,
                   pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2),
                   1-pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2))
  
  conf_mu = ifelse(conf_mu > 0.999, 0.999, ifelse(conf_mu < 0.001,0.001,conf_mu))
  
  conf = numeric(length(conf_mu))
  
  for(i in 1:length(X)){
    conf[i] = rordbeta(n = 1, mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu[i]) + meta_bias + meta_bias_beta * interval[i]), phi = exp(confprec), cutpoints = c(c0, c1))  # Simulated confidence values
  }
  
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = beta,
    sigma_e = sigma_e,
    sigma_k = sigma_k, 
    sigma_m = sigma_m,
    meta_bias = meta_bias,
    confprec = confprec,
    lapse = lapse,
    sigmam_beta = sigmam_beta,
    meta_bias_beta = meta_bias_beta,
    X = X,
    interval = interval,
    ACC = ACC,
    resp = resp,
    c_mu = conf_mu,
    C = conf,
    theta = theta
  ) %>% mutate(trial = 1:n())
  
  
  dataplot = df %>% pivot_longer(cols = c("resp","C")) %>% 
    mutate(ACC = as.factor(ifelse(name == "resp",NA,ACC))) %>% 
    group_by(ACC,name,X) %>% 
    summarize(mean = mean(value),
              se = sd(value) / sqrt(n())) %>% 
    ggplot(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = ACC))+
    geom_pointrange()+
    geom_smooth()+
    facet_wrap(~name,ncol = 1)+
    theme(legend.position = "top")
  
  effectplot = df %>% 
    filter(ACC == 1) %>% 
    mutate(cut_interval = as.factor(cut(interval,3))) %>% 
    group_by(X,cut_interval) %>% 
    summarize(mean = mean(c_mu),
              se = sd(c_mu) / sqrt(n())) %>% 
    ggplot()+
    geom_pointrange(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = (cut_interval)))+
    geom_text(aes(x = 0, y = 0.8, label = paste0("bias = ",round(unique(meta_bias_beta),2), " \n meta_un = ",round(unique(sigmam_beta),2))))+
    geom_smooth(aes(x = X, y = mean, col = (cut_interval)),se = F)+
    theme(legend.position = "top")
  
  plot = dataplot | effectplot 
  
  plot
  
  df$plot <- c(list(plot), rep(list(NA), nrow(df) - 1))
  
  df$id = rnorm(1,0,10000)
  
  return(df)
}



sim_trials_fix = function(n = 300, seed, sigmam_beta,meta_bias_beta){
  
  set.seed(seed)
  
  n_trials <- n        # number of trials
  sigma_e <- -2
  beta <- 0
  
  sigma_m = -1
  confprec = exp(2)
  meta_bias = 0
  lapse = -4  
  sigma_k <- -2
  
  
  c0 = -5
  c1 = 5
  
  
  # -----------------------------
  # Coherence level (stimulus intensity)
  # -----------------------------
  offsets  = seq(-0.2,0.2,length.out = 10)
  w = c(1, 2, 3, 2, 1, 1, 2, 3, 2, 1)
  
  
  stim_values_per_block =
    rep((brms::inv_logit_scaled(
      (brms::inv_logit_scaled(beta) - 0.5) * 2 +
        rep(offsets / 1/sqrt(exp(sigma_e)^2 + exp(sigma_k)^2), times = w)
    ) - 0.5) * 2,3)
  
  
  X = rep(stim_values_per_block,6)
  
  interval = runif(length(X), 1,5)
  
  # -----------------------------
  # Decision
  # -----------------------------
  
  
  
  sigma1 = exp(sigma_e)^2 + exp(sigma_k)^2
  
  
  mu = (X - (brms::inv_logit_scaled(beta)-0.5)*2)
  
  theta = brms::inv_logit_scaled(lapse)/2 + (1-2*brms::inv_logit_scaled(lapse)/2) * pnorm( mu / sqrt(sigma1))
  
  resp = rbinom(length(X), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(X < 0 & resp == 0, 1,
               ifelse(X > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  
  conf_mu = numeric(length(X))
  
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  
  z = mu / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  mills_pos  = exp(dnorm(z, log=TRUE) - pnorm(z, log.p=TRUE))
  mills_neg  = exp(dnorm(z, log=TRUE) - pnorm(-z, log.p=TRUE))
  
  e_cond = ifelse(resp == 1,
                  mu + correction_factor * mills_pos,
                  mu - correction_factor * mills_neg)
  
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  
  sigma2 = exp(sigma_e)^2 + exp(sigma_m +  sigmam_beta * interval)^2
  
  conf_mu = ifelse(resp == 1,
                   pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2),
                   1-pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2))
  
  conf_mu = ifelse(conf_mu > 0.999, 0.999, ifelse(conf_mu < 0.001,0.001,conf_mu))
  
  conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + meta_bias_beta * interval[i])
  
  conf = numeric(length(conf_mu))
  
  # for(i in 1:length(X)){
  # conf[i] = rordbeta(n = 1, mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu[i]) + meta_bias + meta_bias_beta * interval[i]), phi = exp(confprec), cutpoints = c(c0, c1))  # Simulated confidence values
  # }
  
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = beta,
    sigma_e = sigma_e,
    sigma_k = sigma_k, 
    sigma_m = sigma_m,
    meta_bias = meta_bias,
    confprec = confprec,
    lapse = lapse,
    sigmam_beta = sigmam_beta,
    meta_bias_beta = meta_bias_beta,
    X = X,
    interval = interval,
    ACC = ACC,
    resp = resp,
    c_mu = conf_mu,
    C = conf,
    theta = theta
  ) %>% mutate(trial = 1:n())
  
  
  dataplot = df %>% pivot_longer(cols = c("resp","C")) %>% 
    mutate(ACC = as.factor(ifelse(name == "resp",NA,ACC))) %>% 
    group_by(ACC,name,X) %>% 
    summarize(mean = mean(value),
              se = sd(value) / sqrt(n())) %>% 
    ggplot(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = ACC))+
    geom_pointrange()+
    geom_smooth()+
    facet_wrap(~name,ncol = 1)+
    theme(legend.position = "top")
  
  effectplot = df %>% 
    filter(ACC == 1) %>% 
    mutate(cut_interval = as.factor(cut(interval,3))) %>% 
    group_by(X,cut_interval) %>% 
    summarize(mean = mean(c_mu),
              se = sd(c_mu) / sqrt(n())) %>% 
    ggplot()+
    geom_pointrange(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = (cut_interval)))+
    geom_text(aes(x = 0, y = 0.8, label = paste0("bias = ",round(unique(meta_bias_beta),2), " \n meta_un = ",round(unique(sigmam_beta),2))))+
    geom_smooth(aes(x = X, y = mean, col = (cut_interval)),se = F)+
    theme(legend.position = "top")
  
  plot = dataplot | effectplot 
  
  plot
  
  df$plot <- c(list(plot), rep(list(NA), nrow(df) - 1))
  
  df$id = rnorm(1,0,10000)
  
  return(df)
}



########################## 
####### from markdown: ###
##########################


```{r}
map_dfr(results_list,3) %>% filter(div == 0) %>%  
  ggplot(aes(x = simulated, y = mean, ymin = q5, ymax = q95))+geom_pointrange()+geom_abline()+
  ggtitle("Simulated vs Heu")+
  facet_wrap(~variable, scales = "free")



map_dfr(results_list,3) %>% filter(div == 0) %>% filter(variable == "mratio") %>%   
  ggplot(aes(x = simulated, y = mean, ymin = q5, ymax = q95))+geom_pointrange()+geom_abline()+
  ggtitle("Simulated vs Heu (mratio)")+
  facet_wrap(~variable, scales = "free")
```


```{r}
map_dfr(results_list,1) %>% filter(div == 0) %>% ggplot(aes(x = as.factor(N), y = mean))+geom_boxplot()+facet_wrap(~variable, scales = "free")

map_dfr(results_list,2)%>% pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "mean",names_to = "variable")%>% 
  ggplot(aes(x = as.factor(N), y = mean))+geom_boxplot()+facet_wrap(~variable, scales = "free")

inner_join(map_dfr(results_list,3) %>% select(dp,ID),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "mean",names_to = "variable") %>% filter(variable == "dprime")) %>% 
  ggplot(aes(x = dp, y = mean))+geom_point()+xlab("simulated")+geom_abline()


inner_join(map_dfr(results_list,3) %>% select(dp,ID),map_dfr(results_list,1) %>% filter(div == 0) %>% filter(variable == "dp") %>% select(mean,q5,q95,ID,N)) %>% 
  ggplot(aes(x = dp, y = mean,ymin = q5, ymax = q95))+geom_pointrange()+xlab("simulated")+geom_abline()






inner_join(map_dfr(results_list,1) %>% filter(div == 0) %>% filter(variable == "mdp") %>% select(mean,q5,q95,ID,N),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "nonestimean",names_to = "variable") %>% filter(variable == "metaD")) %>% 
  ggplot(aes(x = nonestimean, y = mean,ymin = q5, ymax = q95))+geom_pointrange()+xlab("nonestimean")+geom_abline()

inner_join(map_dfr(results_list,3) %>% select(mdp1,ID),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "mean",names_to = "variable") %>% filter(variable == "metaD")) %>% 
  ggplot(aes(x = mdp1, y = mean))+geom_point()+xlab("simulated")+geom_abline()


inner_join(map_dfr(results_list,3) %>% select(mdp2,ID),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "nonestimean",names_to = "variable") %>% filter(variable == "metaD")) %>% 
  ggplot(aes(x = mdp2, y = nonestimean))+geom_point()+xlab("nonestimean")+geom_abline()

inner_join(map_dfr(results_list,3) %>% select(mdp3,ID),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "nonestimean",names_to = "variable") %>% filter(variable == "metaD")) %>% 
  ggplot(aes(x = mdp3, y = nonestimean))+geom_point()+xlab("nonestimean")+geom_abline()

inner_join(map_dfr(results_list,3) %>% select(mdp4,ID),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "nonestimean",names_to = "variable") %>% filter(variable == "metaD")) %>% 
  ggplot(aes(x = mdp4, y = nonestimean))+geom_point()+xlab("nonestimean")+geom_abline()


inner_join(map_dfr(results_list,3) %>% select(mratio,ID),map_dfr(results_list,2)%>% 
             pivot_longer(cols = c("dprime","c","metaD","Ratio"), values_to = "mean",names_to = "variable") %>% filter(variable == "Ratio")) %>% 
  ggplot(aes(x = mratio, y = mean))+geom_point()+xlab("simulated")+geom_abline()



inner_join(map_dfr(results_list,3) %>% select(mratio,ID),
           map_dfr(results_list,1) %>% filter(div == 0) %>% filter(variable == "mratio") %>% select(mean,q5,q95,ID,N)) %>% 
  ggplot(aes(x = mratio, y = mean))+geom_point()+xlab("simulated")+geom_abline()


```


## Plotting

```{r}
# results_list = readRDS("~/RDM-metacognition/Simulations/Heurestic/Saved models/1000sim_300trials_finalHeu_X1.RData")
# results_list = results_list[-which(results_list == "Error")]

inner_join(
  map_dfr(results_list,1) %>% select(sim_n,sigma_k,sigma_e,sigma_m,bias,conf_prec) %>% 
    mutate(beta = brms::logit_scaled(bias/2 + 0.5), bias = NULL,
           confprec = log(conf_prec), conf_prec = NULL) %>% 
    pivot_longer(cols = c(sigma_k,sigma_e,sigma_m,beta,confprec), values_to = "simulated", names_to = "variable"), 
  map_dfr(results_list,2) %>% 
    filter(rhat < 1.03 & ess_bulk > 400, ess_tail > 400 & div == 0) %>% 
    select(mean,q5,q95,variable,sim_n,div) 
) %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5, ymax = q95))+geom_pointrange(alpha = 0.1)+geom_abline()+
  theme_classic()+
  ylab("Estimated [95 % HDI]")+
  xlab("Simulated value")+
  facet_wrap(~variable, scales = "free")

```

```{r}
inner_join(
  map_dfr(results_list,1),
  map_dfr(results_list,5) 
) %>% 
  ggplot(aes(x = 0.6/sqrt(exp(sigma_e) ^2 +exp(sigma_k)^2), y = dprime))+
  geom_point()+
  geom_abline()

inner_join(
  map_dfr(results_list,1),
  map_dfr(results_list,5) 
) %>% 
  ggplot(aes(x = 1/sqrt(exp(sigma_m)^2+exp(sigma_e)^2), y = metaD))+
  geom_point()+
  geom_abline()


inner_join(
  map_dfr(results_list,1),
  map_dfr(results_list,5) 
) %>% 
  ggplot(aes(x = (1/sqrt(exp(sigma_m)^2+exp(sigma_e)^2)) / 
               (0.6/sqrt(exp(sigma_e)^2 + exp(sigma_k)^2)),
             y = Ratio))+
  geom_point()+
  geom_abline()
```



## brms implementation

```{r}
# estimateda

# standard <- fit_metad(N ~ 1,
#   data = databrms,
#   K = 4,
#   prior = prior(normal(0, 1), class = Intercept) +
#     set_prior(
#       "normal(0, 1)",
#       class = c(
#         "dprime", "c",
#         "metac2zero1diff", "metac2zero2diff", "metac2zero3diff",
#         "metac2one1diff", "metac2one2diff", "metac2one3diff"
#       )
#     )
# )

# brms <- brm(bf(N~1),
#   data = aggregate_metad(databrms, K = 6),
#   family = metad(K = 6),
#   stanvars = stanvars_metad(K = 6),
#    prior = prior(normal(0, 1), class = Intercept) +
#            prior(normal(0, 1), class = dprime) +
#            prior(normal(0, 1), class = c) +
#            prior(lognormal(0, 1), class = metac2zero1diff) +
#            prior(lognormal(0, 1), class = metac2zero2diff) +
#            prior(lognormal(0, 1), class = metac2zero3diff) +
#            prior(lognormal(0, 1), class = metac2one1diff),
#            prior(lognormal(0, 1), class = metac2one2diff),
#            prior(lognormal(0, 1), class = metac2one3diff),
#   backend = "cmdstanr"
#   )


# result = list(simulated = simulated %>% mutate(sim_n = qq),
#               estimateda = estimateda%>% mutate(sim_n = qq, model = "action"),
#               data = df%>% mutate(sim_n = qq),
#               diva = diva,
#               mratio = mratio %>% mutate(sim_n = qq))

```


##############################################################