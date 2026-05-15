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

data = load_data()


# Check exclusion
dd = data %>% 
  mutate(subject = as.numeric(as.factor(subject))) %>% 
  filter(Trialtype == "Main", scale == "conf") %>% 
  filter(
    abs(X) < case_when(
      ID == "RDM_reportz_sub001.csv" ~ 0.6,
      ID == "RDM_reportz_sub010.csv" ~ 0.4,
      TRUE ~ Inf
    )) %>% 
  group_by(subject) %>% 
  mutate(
    trial_n = 1:n(), n_trials = n()
  ) %>% 
  ungroup() %>% 
  mutate(
    across(-c("ID","subject"),
           ~ ifelse(
             ID == "RDM_reportz_sub030.csv" & trial_n > (n_trials - 36),
             NA,
             .
           )
    )
  ) 


exclusions_correct = dd %>% group_by(ID) %>% summarize(n = sum(Correct)/n()) %>% filter(n < 0.60)

exclusions_ntrials = dd %>% 
  select(interTrial.interval, interval,coherence,X,Correct,Y,RT,Confidence,ID) %>% 
  mutate(
    no_response = is.na(Correct) |
      is.na(Confidence) |
      is.na(RT) |
      RT < 0.150) %>% 
  group_by(ID) %>%
  summarize(noresp_rate = mean(no_response)) %>%
  filter(noresp_rate > 0.3)

n_excluded_correctness = length(unique(exclusions_correct$ID))
n_excluded_trials = length(unique(exclusions_ntrials$ID))

n_excluded = length(unique(c(exclusions_correct$ID, exclusions_ntrials$ID)))

dd = dd %>% filter(!ID %in% c(exclusions_correct$ID, exclusions_ntrials$ID)) %>% drop_na()


h = pp_hier_nolapse(fit, dd, 10)

n_bins = 10
df = dd


available <- names(as_draws_df(fit$draws()))


parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","sigmam_beta[1]","meta_bias_beta[1]")
if(length(intersect(parameters, available)) > 5){
params <- intersect(parameters, available)}

n_subj = length(unique(df$ID))
subs = rep(unique(df$ID), length(params))

subj_parameters = str_sub(parameters, 1, -4)

subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")


psycho = function(x,beta,sigma_e, sigma_k){
  return(pnorm((x - (brms::inv_logit_scaled(beta)-0.5)*2) / (sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))))
}

confs = function(action,X,beta,sigma_e,sigma_k,sigma_m){
  
  mu_e = X - (brms::inv_logit_scaled(beta)-0.5)*2
  sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  z = mu_e / sigma_total
  
  # Correction factor
  correction_factor = exp(sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  e_cond = ifelse(action == 1,
                  mu_e + correction_factor * dnorm(z) / pnorm(z),
                  mu_e - correction_factor * dnorm(z) / pnorm(-z))  # Note: pnorm(-z) here!
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  sigma2 = exp(sigma_m)^2 + exp(sigma_e)^2
  
  c_mu = ifelse(action == 1,
                pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2),
                1-pnorm((1/(sqrt(1 + (Cc^2*1^2) / sigma2))) * (e_cond * 1 * Cc) / sigma2))
  
  return(c_mu)
  
}

sum = fit$summary(c("gm","tau_u")) %>% dplyr::select(-c(mad,median))
sub = fit$summary(c(subj_parameters)) %>% dplyr::select(-c(mad,median)) %>% mutate(ID = subs)

pred_data_group = as_draws_df(fit$draws("gm")) %>% select(-contains(".")) %>% 
  rename_with(~c(  "beta",
                   "sigma_e",
                   "sigma_k",
                   "sigma_m",
                   "meta_bias",
                   "confprec",
                   "sigmam_beta",
                   "meta_bias_beta")) %>%
  mutate(draw = 1:n()) %>% 
  mutate(x = list(seq(-1,1,by = 0.05))) %>% 
  unnest() %>% 
  mutate(interval = mean(df$interval)) %>% 
  unnest() %>% 
  group_by(draw) %>% 
  mutate(p = psycho(x, beta,sigma_e,sigma_k)) %>% 
  # resp = rbinom(n(),1,p)) %>% 
  # mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
  rowwise() %>% 
  mutate(action = rbinom(1,1,p)) %>% 
  unnest() %>% 
  rowwise() %>% 
  mutate(conf_mu = ifelse(action == 1,
                          confs(1,x,beta, sigma_e, sigma_k,sigma_m + interval * sigmam_beta ),
                          confs(0,x,beta, sigma_e, sigma_k,sigma_m+ interval * sigmam_beta ))) %>% 
  mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + interval * meta_bias_beta)) %>% 
  ungroup()




df1 = bind_rows(
  pred_data_group %>%
    mutate(action = ifelse(action == 1, "up", "down")) %>%
    group_by(x) %>%
    summarize(
      name = "Type-1",
      mean = mean(p, na.rm = T),
      q5 = quantile(p, 0.05),
      q95 = quantile(p, 0.95),
      .groups = "drop"
    ),
  pred_data_group %>%
    mutate(action = ifelse(action == 1, "up", "down")) %>%
    group_by(x, action) %>%
    summarize(name = "Confidence",
              mean = mean(conf_mu, na.rm = T),
              q5 = quantile(conf_mu, 0.05),
              q95 = quantile(conf_mu, 0.95),
              .groups = "drop")
) 



# Prepare observed data
if (!is.null(n_bins)) {
  # Create common bin boundaries based on the range of both datasets
  df_g <- df %>%
    mutate(
      X_bin = cut(
        X,
        breaks = seq(min(X, na.rm = TRUE),
                     max(X, na.rm = TRUE),
                     length.out = n_bins + 1),
        labels = FALSE,
        include.lowest = TRUE
      ),
      # compute bin centers per subject
      X = {
        breaks_i <- seq(min(X, na.rm = TRUE),
                        max(X, na.rm = TRUE),
                        length.out = n_bins + 1)
        centers_i <- (breaks_i[-1] + breaks_i[-length(breaks_i)]) / 2
        centers_i[X_bin]
      }
    ) %>%
    ungroup() %>%
    select(-X_bin)
}else{
  df_g = df
}


behplot = rbind(
  bin = df_g %>%
    # mutate(action = ifelse(Y == 1, "up", "down")) %>%
    mutate(action = NA) %>%
    group_by(X,action) %>%
    summarize(
      name = "Type-1",
      k = sum(Y),
      n = n(),
      mean = k / n,
      
      a_post = k + 1,
      b_post = (n - k) + 1,
      
      q5 = qbeta(0.025, a_post, b_post),
      q95 = qbeta(0.975, a_post, b_post),
      .groups = "drop") %>% mutate(k = NULL, n = NULL, a_post = NULL, b_post = NULL),
  df_g %>%
    mutate(action = ifelse(Y == 1, "up", "down")) %>%
    group_by(X, action) %>%
    summarize(name = "Confidence",
              mean = mean(Confidence, na.rm = T),
              q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
              q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
              .groups = "drop")
) 


group_plot = behplot%>% 
  rename(variable = name) %>% 
  #filter(name != "RT") %>% 
  ggplot() +
  geom_pointrange(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
                  shape = 21, color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
  geom_line(data = df1%>% rename(variable = name), aes(x = x, y = mean, col = action))+
  geom_ribbon(data = df1%>% rename(variable = name), aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)+
  facet_wrap(~variable, ncol = 1, scales = "free")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme_classic(base_size = 14) +
  labs(color = "action", fill = "action",
       y = "Value")

bounds=df %>% group_by(ID) %>% summarize(min = min(X)-0.05, max = max(X)+0.05, by = ((max-min) / 20))



pred_subj = inner_join(as_draws_df(fit$draws(c(  "beta",
                                                 "sigma_e",
                                                 "sigma_k",
                                                 "sigma_m",
                                                 "meta_bias",
                                                 "confprec",
                                                 "sigmam_beta",
                                                 "meta_bias_beta"
))) %>% select(-contains(".")) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(
    sub_idx = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param   = str_remove(name, "\\[\\d+\\]"),
    ID  = subs[sub_idx],
    name = NULL
  ),bounds) %>% 
  pivot_wider(names_from = "param", values_from = "value") %>% 
  group_by(ID) %>% 
  mutate(x = list(seq(unique(min),unique(max),by = unique(by)))) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(interval = mean(df$interval)) %>% 
  unnest() %>% 
  group_by(draw, ID) %>% 
  mutate(p = psycho(x, beta,sigma_e,sigma_k)) %>% 
  rowwise() %>% 
  mutate(action = rbinom(1,1,p)) %>% 
  unnest() %>% 
  rowwise() %>% 
  mutate(conf_mu = ifelse(action == 1,
                          confs(1,x,beta, sigma_e, sigma_k,sigma_m + interval * sigmam_beta ),
                          confs(0,x,beta, sigma_e, sigma_k,sigma_m+ interval * sigmam_beta ))) %>% 
  mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + interval * meta_bias_beta)) %>% 
  ungroup()



df1_sub = bind_rows(
  pred_subj %>%
    mutate(action = ifelse(action == 1, "up", "down")) %>%
    group_by(x,ID) %>%
    summarize(
      name = "Type-1",
      mean = mean(p, na.rm = T),
      q5 = quantile(p, 0.05),
      q95 = quantile(p, 0.95),
      .groups = "drop"
    ),
  pred_subj %>%
    mutate(action = ifelse(action == 1, "up", "down")) %>%
    group_by(x,ID, action) %>%
    summarize(name = "Confidence",
              mean = mean(conf_mu, na.rm = T),
              q5 = quantile(conf_mu, 0.05),
              q95 = quantile(conf_mu, 0.95),
              .groups = "drop")
)  %>% mutate(variable = name, name = NULL)




# Prepare observed data
if (!is.null(n_bins)) {
  # Create common bin boundaries based on the range of both datasets
  df_s <- df %>%
    group_by(ID) %>% 
    mutate(
      X_bin = cut(
        X,
        breaks = seq(min(X, na.rm = TRUE),
                     max(X, na.rm = TRUE),
                     length.out = n_bins + 1),
        labels = FALSE,
        include.lowest = TRUE
      ),
      # compute bin centers per subject
      X = {
        breaks_i <- seq(min(X, na.rm = TRUE),
                        max(X, na.rm = TRUE),
                        length.out = n_bins + 1)
        centers_i <- (breaks_i[-1] + breaks_i[-length(breaks_i)]) / 2
        centers_i[X_bin]
      }
    ) %>%
    ungroup() %>%
    select(-X_bin)
}else{
  df_s = df
}



behplot = rbind(
  bin = df_s %>%
    mutate(action = NA) %>%
    group_by(X,ID,action) %>%
    summarize(
      name = "Type-1",
      k = sum(Y),
      n = n(),
      mean = k / n,
      
      a_post = k + 1,
      b_post = (n - k) + 1,
      
      q5 = qbeta(0.025, a_post, b_post),
      q95 = qbeta(0.975, a_post, b_post),
      .groups = "drop") %>% mutate(k = NULL, n = NULL, a_post = NULL, b_post = NULL),
  df_s %>%
    mutate(action = ifelse(Y == 1, "up", "down")) %>%
    group_by(X,ID, action) %>%
    summarize(name = "Confidence",
              mean = mean(Confidence, na.rm = T),
              q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
              q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
              .groups = "drop")
) %>% mutate(variable = name, name = NULL)




sub_plots = list()
qq = 0
for(name in unique(df1_sub$ID)){
  
  xmin = min(behplot%>% filter(ID == name) %>% .$X) - 0.1
  xmax = max(behplot%>% filter(ID == name) %>% .$X) + 0.1
  
  parameters = pred_subj %>% filter(ID == name) %>% filter(x > xmin & x < xmax) %>% 
    pivot_longer(cols = c("sigma_e","sigma_m","sigma_k")) %>%
    select(ID,name,value) %>% distinct() %>% 
    mutate(name = ifelse(name == "sigma_e","S_e",ifelse(name == "sigma_m","S_m",ifelse(name == "sigma_k","S_c",name)))) %>% 
    group_by(name) %>% summarize(mean = mean(value, na.rm = T),
                                 q5 = quantile(value,0.05, na.rm = T),
                                 q95 = quantile(value,0.95, na.rm = T)) %>% 
    mutate(label = paste0(name," = ",round(mean,2)," [",round(q5,2)," ; ",round(q95,2),"]")) %>% 
    mutate(variable = "Confidence") %>% mutate(param = c(0,0.15,0.3))
  
  qq = qq+1
  
  
  sub_plot = behplot  %>% 
    filter(variable != "RT" & ID == name)  %>% filter(X > xmin & X < xmax) %>% 
    ggplot() +
    geom_pointrange(data = behplot%>% filter(variable != "RT"& ID == name), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
                    shape = 21, color = "black", alpha = 0.5) +
    facet_wrap(~variable, scales = "free_y", ncol = 1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    theme_classic(base_size = 14) +
    labs(color = "action", fill = "action",
         y = "P(a=1) Or P(a=D)", x = "Signed Stimulus (XD)") +
    geom_vline(xintercept = 0, linetype = 2) +
    # labs(subtitle = name)+
    labs(subtitle = paste0("ID =", qq))+
    # scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
    # scale_x_continuous(limits = c(xmin,xmax), breaks = scales::pretty_breaks(n = 5))+
    
    # geom_text(data = parameters, aes(x = 0, y = param,label = label), size = 4)+
    geom_line(data = df1_sub %>% filter(ID == name), aes(x = x, y = mean, col = action))+
    geom_ribbon(data = df1_sub%>% filter(ID == name), aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)
  
  
  sub_plots[[name]] = sub_plot
  
}


n_bins = 10
# Prepare observed data
if (!is.null(n_bins)) {
  # Create common bin boundaries based on the range of both datasets
  df_g <- df %>%
    mutate(
      X_bin = cut(
        X,
        breaks = seq(min(X, na.rm = TRUE),
                     max(X, na.rm = TRUE),
                     length.out = n_bins + 1),
        labels = FALSE,
        include.lowest = TRUE
      ),
      # compute bin centers per subject
      X = {
        breaks_i <- seq(min(X, na.rm = TRUE),
                        max(X, na.rm = TRUE),
                        length.out = n_bins + 1)
        centers_i <- (breaks_i[-1] + breaks_i[-length(breaks_i)]) / 2
        centers_i[X_bin]
      }
    ) %>%
    ungroup() %>%
    select(-X_bin)
}else{
  df_g = df
}


behplot = rbind(
  bin = df_g %>%
    # mutate(action = ifelse(Y == 1, "up", "down")) %>%
    mutate(action = NA) %>%
    group_by(X,action) %>%
    summarize(
      name = "Type-1",
      k = sum(Y),
      n = n(),
      mean = k / n,
      
      a_post = k + 1,
      b_post = (n - k) + 1,
      
      q5 = qbeta(0.025, a_post, b_post),
      q95 = qbeta(0.975, a_post, b_post),
      .groups = "drop") %>% mutate(k = NULL, n = NULL, a_post = NULL, b_post = NULL),
  df_g %>%
    mutate(action = ifelse(Y == 1, "up", "down")) %>%
    group_by(X, action) %>%
    summarize(name = "Confidence",
              mean = mean(Confidence, na.rm = T),
              q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
              q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
              .groups = "drop")
) 




pred_subj = inner_join(as_draws_df(fit$draws(c(  "beta",
                                                 "sigma_e",
                                                 "sigma_k",
                                                 "sigma_m",
                                                 "meta_bias",
                                                 "confprec",
                                                 "sigmam_beta",
                                                 "meta_bias_beta"
))) %>% select(-contains(".")) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(
    sub_idx = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param   = str_remove(name, "\\[\\d+\\]"),
    ID  = subs[sub_idx],
    name = NULL
  ),bounds) %>% 
  pivot_wider(names_from = "param", values_from = "value") %>% 
  group_by(ID) %>% 
  mutate(x = list(seq(-1,1,by = 0.05))) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(interval = mean(df$interval)) %>% 
  unnest() %>% 
  group_by(draw, ID) %>% 
  mutate(p = psycho(x, beta,sigma_e,sigma_k)) %>% 
  rowwise() %>% 
  mutate(action = rbinom(1,1,p)) %>% 
  unnest() %>% 
  rowwise() %>% 
  mutate(conf_mu = ifelse(action == 1,
                          confs(1,x,beta, sigma_e, sigma_k,sigma_m + interval * sigmam_beta ),
                          confs(0,x,beta, sigma_e, sigma_k,sigma_m+ interval * sigmam_beta ))) %>% 
  mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + interval * meta_bias_beta)) %>% 
  ungroup()


df1_sub = bind_rows(
  pred_subj %>%
    mutate(action = ifelse(action == 1, "up", "down")) %>%
    group_by(ID,x) %>%
    summarize(
      name = "Type-1",
      mean = mean(p, na.rm = T),
      q5 = quantile(p, 0.05),
      q95 = quantile(p, 0.95),
      .groups = "drop"
    ),
  pred_subj %>%
    mutate(action = ifelse(action == 1, "up", "down")) %>%
    group_by(ID,x, action) %>%
    summarize(name = "Confidence",
              mean = mean(conf_mu, na.rm = T),
              q5 = quantile(conf_mu, 0.05),
              q95 = quantile(conf_mu, 0.95),
              .groups = "drop")
)  %>% mutate(variable = name, name = NULL)


average_sub_average_group = df1_sub %>% group_by(x,variable,action) %>% summarize(mean = mean(mean), q5 = mean(q5),q95 = mean(q95)) %>% 
  ggplot(aes(x = x, y = mean, col = action))+geom_line()+geom_ribbon(aes(ymin  = q5, ymax = q95, fill = action), alpha = 0.5)+
  geom_pointrange(data = behplot%>% filter(name != "RT") %>% rename(variable = name), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
                  shape = 21, color = "black", alpha = 0.5)+
  theme_classic(base_size = 14) +
  labs(color = "action", fill = "action", x = "Signed Stimulus (XD)") +
  geom_vline(xintercept = 0, linetype = 2) +
  # labs(subtitle = name)+
  facet_wrap(~variable, ncol = 1, scales = "free")



  
