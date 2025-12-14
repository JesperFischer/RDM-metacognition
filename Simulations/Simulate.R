
simulate_data = function(S = 10, N = 50){
  
  fit <- readRDS("~/Multivariate-Confidence-modeling/Saved models/VMP_extero_cohort1.rds")
  # Extract posterior means
  gm_draws <- fit$draws("gm") %>% as_draws_matrix() %>% colMeans()
  tau_draws <- fit$draws("tau_u") %>% as_draws_matrix() %>% colMeans()
  cor_draws <- fit$draws("correlation_matrix") %>% as_draws_matrix() %>% colMeans()
  
  cor_mat = matrix(ncol = 12, nrow = 12)
  k = 0
  for(i in 1:12){
    for(j in 1:12){
      k = k+1
      cor_mat[i,j] = cor_draws[k]
    }
  }
  
  
  # Scale by tau
  Sigma <- diag(tau_draws) %*% cor_mat %*% diag(tau_draws)  # covariance matrix
  
  # Simulate P parameters for S subjects
  indi_params <- MASS::mvrnorm(n = S, mu = gm_draws, Sigma = Sigma)
  
  colnames(indi_params) <- c("alpha","beta","lapse","rt_int","rt_slope","rt_prec","rt_stim","conf_int","conf_ACC","conf_entropy","conf_entropy_ACC", "conf_prec")
  
  rho13 = data.frame(cop_cor = fit$draws(c("rho_p_conf"))%>% as_draws_matrix() %>%
                       colMeans()) %>%
    summarize(mean = mean(cop_cor),
              sd = sd(cop_cor),
              q5 = quantile(cop_cor,0.05),
              q95 = quantile(cop_cor,0.95))
  
  
  rho12 = data.frame(cop_cor = fit$draws(c("rho_p_rt"))%>% as_draws_matrix() %>%
                       colMeans()) %>%
    summarize(mean = mean(cop_cor),
              sd = sd(cop_cor),
              q5 = quantile(cop_cor,0.05),
              q95 = quantile(cop_cor,0.95))
  
  
  rho23 = data.frame(cop_cor = fit$draws(c("rho_rt_conf"))%>% as_draws_matrix() %>%
                       colMeans()) %>%
    summarize(mean = mean(cop_cor),
              sd = sd(cop_cor),
              q5 = quantile(cop_cor,0.05),
              q95 = quantile(cop_cor,0.95))
  
  
  c11 = data.frame(cop_cor = fit$draws(c("c11"))%>% as_draws_matrix() %>%
                     colMeans()) %>%
    summarize(mean = mean(cop_cor),
              sd = sd(cop_cor),
              q5 = quantile(cop_cor,0.05),
              q95 = quantile(cop_cor,0.95))
  
  
  c00 = data.frame(cop_cor = fit$draws(c("c0"))%>% as_draws_matrix() %>%
                     colMeans()) %>%
    summarize(mean = mean(cop_cor),
              sd = sd(cop_cor),
              q5 = quantile(cop_cor,0.05),
              q95 = quantile(cop_cor,0.95))
  
  
  indi_params = as.data.frame(indi_params) %>% mutate(beta = exp(beta),
                                                      lapse = brms::inv_logit_scaled(lapse) / 2,
                                                      rt_prec = exp(rt_prec),
                                                      conf_prec = exp(conf_prec)) %>%
    mutate(rt_ndt = rnorm(n(),0.35,0.05),
           rho12 = rnorm(n(),rho12$mean,rho12$sd),
           rho23 = rnorm(n(),rho23$mean,rho23$sd),
           rho13 = rnorm(n(),rho13$mean,rho13$sd),
           c0 = rnorm(n(),c00$mean,c00$sd),
           c1 = c0 + exp(rnorm(n(),c11$mean,c11$sd)),
    ) %>%
    mutate(subject = 1:n())
  
  
  generate_trials = function(params,N){
    
    x = rnorm(N,0,4)
    
    p = psychometric(x,params)
    rt_mu = RT_mean(p,params)
    
    us = get_copula_vec(params, length(x))
    
    bin = qbinom(us$u_bin,1,p)
    ACC = ifelse(x > 0 & bin == 1, 1, ifelse(x < 0 & bin == 0,1,0))
    conf_mu = conf_mean(p,ACC,params)
    
    rts = qlnorm(us$u_rt,rt_mu, params$rt_prec) + params$rt_ndt
    
    conf = qordbeta(us$u_vas,
                    brms::inv_logit_scaled(conf_mu),
                    params$conf_prec,
                    params$c0,
                    params$c1)
    
    
    predictions = data.frame(bin = bin, rts = rts,x = x, ACC = ACC,
                             prob = p, rt_mu = exp(rt_mu) + params$rt_ndt,
                             conf = conf, mu_conf = brms::inv_logit_scaled(conf_mu))
    
  }
  
  
  df = indi_params %>% rowwise() %>% mutate(resps = list(generate_trials(cur_data(), N = N)))
  
  plot = df %>%
    unnest(resps) %>%
    pivot_longer(cols = c("bin","rts","conf"), names_to = "name", values_to = "value") %>%
    pivot_longer(cols = c("prob","rt_mu","mu_conf"), names_to = "means", values_to = "mean_value") %>%
    mutate(
      means = ifelse(means == "prob","bin",ifelse(means == "rt_mu","rts","conf"))
    ) %>%
    filter(name == means) %>%
    mutate(ACC = as.factor(ACC)) %>%
    ggplot(aes(x = x, y = value)) +
    # conditional coloring: only color if name == "conf"
    geom_point(aes(col = as.factor(ifelse(name == "conf", ACC, NA))), alpha = 0.6, show.legend = TRUE) +
    geom_line(aes(y = mean_value, col = as.factor(ifelse(name == "conf", ACC, NA))), show.legend = FALSE) +
    facet_grid(name ~ subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 8)
    ) +
    labs(x = "Stimulus (x)", y = "Response", col = "ACC")
  
  return(list(df %>%
                unnest(resps), plot))
  
}


# simulate_data()


fit_data = function(){
  
  qqs = simulate_data()
  qqs[[2]]
  
  qq = qqs[[1]]
  
  # qq = qq[[1]]   %>% filter(subject == 1)
  t_p_s = qq %>% group_by(subject) %>% summarize(n = n())
  
  ends <- cumsum(t_p_s$n)
  
  # Calculate the start points
  starts <- c(1, head(ends, -1) + 1)
  
  # mod = cmdstan_model(here::here("Stanmodels","ss","SS_Confidence_Standard.stan"))
  mod = cmdstan_model(here::here("Stanmodels","Confidence_Standard.stan"))
  
  aa <- mod$sample(
    data = list(N = nrow(qq),
                S = length(unique(qq$subject)),
                starts = starts,
                minRT = unique(qq %>% group_by(subject) %>% summarize(minRT = min(rts)) %>% .$minRT),
                ends = ends,
                t_p_s = t_p_s$n,
                X = qq$x,
                ACC = qq$ACC,
                S_id = qq$subject,
                RT = (qq$rts),
                Conf = qq$conf,
                binom_y = qq$bin),
    refresh = 10,
    init = 0,
    iter_sampling = 500,
    iter_warmup = 500,
    adapt_delta = 0.95,
    parallel_chains = 4)
  
  
  
  param_names = c("alpha","beta","lapse","rt_int","rt_slope","rt_prec","rt_stim",
                  "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC", "conf_prec",
                  "rho_p_rt","rho_p_conf","rho_rt_conf","rt_ndt","c0","c1")
  
  
  fitted_group = aa$summary(c("gm","tau_u")) %>% mutate(simulated = c(gm_draws,tau_draws)) %>% mutate(type = c(rep("means",12), rep("taus",12)))
  
  fitted_group %>% ggplot(aes(x = simulated, y = mean, ymin = q5, ymax =q95))+geom_pointrange()+geom_abline()+facet_wrap(~type, scales = "free")
  
  fitted_subj = aa$summary(param_names) %>%
    mutate(
      subject = str_extract(variable, "(?<=\\[)\\d+(?=\\])"),       # extract number inside brackets
      subject = as.integer(subject),                             # convert to integer
      variable = str_remove(variable, "\\[.*\\]")                   # remove [number]
    )
  
  simulated = qq %>% rename(rho_rt_conf = rho23, rho_p_rt = rho12, rho_p_conf = rho13) %>%
    select(param_names,subject) %>% distinct() %>%
    pivot_longer(param_names, names_to = "variable", values_to = "simulated")
  
  
  fitted_subj = inner_join(fitted_subj,simulated)
  
  fitted_subj %>% ggplot(aes(x = simulated, y = mean, ymin = q5, ymax =q95))+geom_pointrange()+geom_abline()+facet_wrap(~variable, scales = "free")
  
}


psychometric = function(x,df){
  df$lapse + (1-2*df$lapse)*brms::inv_logit_scaled(df$beta * (x-df$alpha))
}

entropy = function(p){
  -p * log(p) - (1-p) * log(1-p)
}

RT_mean = function(p,df){
  
  entropy_t = entropy(p)
  df$rt_int + entropy_t * df$rt_slope
  
}

conf_mean = function(p,ACC,df){
  
  entropy_t = entropy(p)
  df$conf_int + entropy_t * df$conf_entropy + df$conf_ACC * ACC + df$conf_entropy_ACC * ACC * entropy_t
  
}


get_copula_vec = function(df,n){
  
  Sigma <- matrix(c(
    1,       df$rho12, df$rho13,
    df$rho12,  1,      df$rho23,
    df$rho13,   df$rho23, 1
  ), ncol = 3)
  
  df = data.frame(mnormt::rmnorm(n = n, mean = c(0,0,0),
                                 varcov = Sigma)) %>% rename(Bin = X1,
                                                             RT = X2,
                                                             VAS = X3)
  
  us_1 = pnorm(df$Bin)
  us_2 = pnorm(df$RT)
  us_3 = pnorm(df$VAS)
  
  return(data.frame(u_bin = us_1,
                    u_rt = us_2,
                    u_vas = us_3))
}


