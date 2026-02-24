library(future)
library(future.apply)
library(ordbetareg)

pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes,
               brms, patchwork, cowplot,ggpubr,flextable)



extract_Lu = function(row, P) {
  row %>%
    select(matches("^L_u\\.")) %>%
    as.numeric() %>%
    matrix(nrow = P, ncol = P)
}
get_bf = function(prior,post){
  library(logspline)
  post_fit  <- logspline(post)
  prior_fit <- logspline(prior)
  
  post_at_0  <- dlogspline(0, post_fit)
  prior_at_0 <- dlogspline(0, prior_fit)
  
  BF_01 <- prior_at_0 / post_at_0
  
  return(BF_01)
}

sim_trials = function(df){
  
  c0 = -abs(rnorm(1,5,1))            # Cut points for ordered beta regression, lower bound
  c1 = abs(rnorm(1,5,1))             # Cut points for ordered beta regression, upper bound
  
  # -----------------------------
  # Coherence level (stimulus intensity)
  # -----------------------------
  offsets  = seq(-0.2,0.2,length.out = 10)
  w = c(1, 2, 3, 2, 1, 1, 2, 3, 2, 1)
  
  
  stim_values_per_block =
    rep((brms::inv_logit_scaled(
      (brms::inv_logit_scaled(df$beta) - 0.5) * 2 +
        rep(offsets / 1/sqrt(exp(df$sigma_e)^2 + exp(df$sigma_k)^2), times = w)
    ) - 0.5) * 2,3)
  
  
  X = rep(stim_values_per_block,6)
  
  interval = runif(length(X), 1,5)
  
  # -----------------------------
  # Decision
  # -----------------------------
  
  
  
  sigma1 = exp(df$sigma_e)^2 + exp(df$sigma_k)^2


  mu = (X - (brms::inv_logit_scaled(df$beta)-0.5)*2)
  
  theta = brms::inv_logit_scaled(df$lapse)/2 + (1-2*brms::inv_logit_scaled(df$lapse)/2) * pnorm( mu / sqrt(sigma1))
  
  resp = rbinom(length(X), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(X < 0 & resp == 0, 1,
               ifelse(X > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------

  conf_mu = numeric(length(X))
  
  sigma_total = sqrt(exp(df$sigma_k)^2 + exp(df$sigma_e)^2)
  
  z = mu / sigma_total
  
  # Correction factor
  correction_factor = exp(df$sigma_e)^2 / sigma_total
  
  # Conditional expectation (different Mills ratios for each choice)
  mills_pos  = exp(dnorm(z, log=TRUE) - pnorm(z, log.p=TRUE))
  mills_neg  = exp(dnorm(z, log=TRUE) - pnorm(-z, log.p=TRUE))
  
  e_cond = ifelse(resp == 1,
                  mu + correction_factor * mills_pos,
                  mu - correction_factor * mills_neg)
  
  
  # Now use e_cond in the confidence formula
  Cc = 2/1.7
  
  sigma2 = exp(df$sigma_e)^2 + exp(df$sigma_m +  df$sigmam_beta * interval)^2
  
  conf_mu = ifelse(resp == 1,
                pnorm((1/(sqrt(1 + (Cc^2*abs(X)^2) / sigma2))) * (e_cond * abs(X) * Cc) / sigma2),
                1-pnorm((1/(sqrt(1 + (Cc^2*abs(X)^2) / sigma2))) * (e_cond * abs(X) * Cc) / sigma2))
  
  # conf_mu = ifelse(conf_mu > 0.999, 0.999, ifelse(conf_mu < 0.001,0.001,conf_mu))
  
  conf = numeric(length(conf_mu))
  
  for(i in 1:length(X)){
    conf[i] = rordbeta(n = 1, mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu[i]) + df$meta_bias + df$meta_bias_beta * interval[i]), phi = exp(df$confprec), cutpoints = c(c0, c1))  # Simulated confidence values
  }
  
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = df$beta,
    sigma_e = df$sigma_e,
    sigma_k = df$sigma_k, 
    sigma_m = df$sigma_m,
    meta_bias = df$meta_bias,
    confprec = df$confprec,
    lapse = df$lapse,
    sigmam_beta = df$sigmam_beta,
    meta_bias_beta = df$meta_bias_beta,
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
    ggplot(aes(x = X, y = mean, ymin = mean-2*se, ymax = mean+2*se, col = (ACC)))+
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
    geom_text(aes(x = 0, y = 0.8, label = paste0("bias = ",round(unique(df$meta_bias_beta),2), " \n meta_un = ",round(unique(df$sigmam_beta),2))))+
    geom_smooth(aes(x = X, y = mean, col = (cut_interval)),se = F)+
    theme(legend.position = "top")
  
  plot = dataplot | effectplot 
  
  df$plot <- c(list(plot), rep(list(NA), nrow(df) - 1))
  
  return(df)
}

simulate_subjects = function(draw_row, S = 20, P = 9) {
  
  gm =
    draw_row %>%
    select(matches("^gm\\.")) %>%
    as.numeric()
  
  tau =
    draw_row %>%
    select(matches("^tau_u\\.")) %>%
    as.numeric()
  
  L_u = extract_Lu(draw_row, P)
  
  map_dfr(seq_len(S), function(i) {
    
    z_new = rnorm(P)
    
    param_new =
      gm + (diag(tau) %*% L_u %*% z_new)[,1]
    
    tibble(
      subj_id        = i,
      beta         = param_new[1],
      sigma_e          = param_new[2],
      sigma_k     = param_new[3],
      sigma_m   = param_new[4],
      meta_bias   = param_new[5],
      lapse      = param_new[6],
      confprec          = param_new[7],
      sigmam_beta   = param_new[8],
      meta_bias_beta = param_new[9]
    )
  })
  
}

fitter = function(new_subjects,subs){
  
  sim_id = rnorm(1,0,100)
  
  dd = new_subjects %>% select(sim_d,subj_id) %>% unnest() %>% filter(subj_id %in% 1:subs)
  
  model = cmdstanr::cmdstan_model(here::here("Stanmodels","Heurestic models",
                                             "Emperical dataanalysis","hierarchical","Full_Heurestic.stan"))
  
  
  dd =  dd %>%
    arrange(subj_id)
  
  starts_ends = dd %>%
    mutate(row = row_number()) %>%
    group_by(subj_id) %>%
    summarise(
      start = first(row),
      end   = last(row),
      .groups = "drop"
    )
  
  
  datastan = list(N = nrow(dd),
                  S = length(unique(dd$subj_id)),
                  starts = starts_ends$start,
                  ends = starts_ends$end,
                  S_id = as.numeric(as.factor(dd$subj_id)),
                  a = dd$resp,
                  C = dd$C,
                  X = abs(dd$X),
                  XD = dd$X,
                  interval = dd$interval,
                  ACC = dd$ACC)
  
  
  pf = model$pathfinder(data = datastan)
  

  fit <-model$sample(
    data = datastan,
    refresh = 10,
    iter_sampling = 500,
    iter_warmup = 500,
    adapt_delta = 0.95,
    max_treedepth = 12,
    init  = pf,
    parallel_chains = 4)
  
  prior_sd = 0.5
  
  divs = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent),
                                                            mean_tree = mean(num_max_treedepth),
                                                            mean_energi = mean(ebfmi))
  
  bf8 = get_bf(post = as.numeric(as_draws_df(fit$draws("gm[8]")) %>% .$`gm[8]`),prior = rnorm(4000,0,prior_sd))
  bf9 = get_bf(post = as.numeric(as_draws_df(fit$draws("gm[9]")) %>% .$`gm[9]`),prior = rnorm(4000,0,prior_sd))
  
  BF = data.frame(bf8 = bf8, bf9 = bf9, n_subs = subs,sim_id = sim_id,
                  div = divs$mean_div, tree = divs$mean_tree, energi = divs$mean_energi, draw_id = unique(new_subjects$draw),prior_sd = prior_sd)
  
  
  sim_subj = dd %>% select(c0,c1,beta,sigma_e,sigma_k,sigma_m,meta_bias,confprec,lapse,sigmam_beta,meta_bias_beta,subj_id) %>% distinct() %>% 
    pivot_longer(-subj_id, names_to = "variable",values_to = "simulated")
  
  esti_subj = fit$summary(c("c0","c1","beta","sigma_e","sigma_k","sigma_m","meta_bias","confprec","lapse","sigmam_beta","meta_bias_beta")) %>% 
    mutate(
      subj_id = str_extract(variable, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
      variable   = str_remove(variable, "\\[\\d+\\]")
    ) 
  
  subj_param = inner_join(sim_subj, esti_subj)%>% mutate(n_subs = subs, sim_id = sim_id, div = divs$mean_div,
                                                         tree = divs$mean_tree, energi = divs$mean_energi, draw_id = unique(new_subjects$draw),prior_sd = prior_sd)
  
  sim_group = new_subjects %>% 
    select(paste0("tau_u.",1:9,"."),paste0("gm.",1:9,".")) %>%
    distinct() %>% 
    pivot_longer(everything(), names_to = "variable",values_to = "simulated") %>% 
    mutate(variable = sub("\\.(\\d+)\\.$", "[\\1]", variable))
  
  esti_group = fit$summary(c("gm","tau_u"))
  
  group_param = inner_join(sim_group,esti_group) %>% mutate(n_subs = subs, sim_id = sim_id, div = divs$mean_div,
                                                            tree = divs$mean_tree, energi = divs$mean_energi, draw_id = unique(new_subjects$draw),prior_sd = prior_sd)
  
  return(list(BF,subj_param,group_param))
  
}


sim_new_subjects = function(max_subjects = 100){
  
  make_draws = function(){
    # groupmean = fit_hier$draws(c("gm")) %>% as_draws_df() %>% select(-contains("."))
    # between_subj = fit_hier$draws(c("tau_u")) %>% as_draws_df()%>% select(-contains("."))
    # Lu = fit_hier$draws(c("L_u")) %>% as_draws_df()%>% select(-contains("."))
    # 
    # post_draws =
    #   groupmean %>%
    #   mutate(draw = row_number()) %>%
    #   left_join(
    #     between_subj %>% mutate(draw = row_number()),
    #     by = "draw"
    #   ) %>%
    #   left_join(
    #     Lu %>% mutate(draw = row_number()),
    #     by = "draw"
    #   )
    # 
    # write.csv(post_draws,here::here("Simulations","Heurestic","Sequential Sampling","post_draws.csv"))
  }

  post_draws = read.csv(here::here("Simulations","Heurestic","Sequential Sampling","post_draws.csv"))
  
  draws_id = sample(nrow(post_draws), 1)

  new_subjects =
    post_draws %>% filter(draw %in% draws_id) %>% 
    mutate(`gm.8.` = 0,
           `gm.9.` = 0) %>%
    rowwise() %>% 
    mutate(sim = list(simulate_subjects(cur_data(), S = max_subjects,P = 9))) %>% 
    unnest(sim) %>% 
    ungroup() %>% 
    rowwise() %>% 
    mutate(sim_d = list(sim_trials(cur_data())))
  
  
  
  # new_subjects %>% select(sim_d,subj_id) %>% unnest(sim_d) %>% 
  #   group_by(X,subj_id) %>% 
  #   summarize(mean = mean(resp)) %>% 
  #   ggplot(aes(x = X, y = mean))+
  #   facet_wrap(~subj_id)+
  #   geom_point()
  
  # new_subjects %>% select(sim_d,subj_id) %>% unnest(sim_d) %>% 
  #       group_by(X,subj_id,ACC) %>% 
  #   summarize(mean = mean(C)) %>% 
  #   ggplot(aes(x = X, y = mean, col = as.factor(ACC)))+
  #   facet_wrap(~subj_id)+
  #   geom_point()
  
  fits = fitter(new_subjects,5)
  
  
  
  
  plan(multisession, workers = 5)  # Windows-friendly
  
  n_subj <- c(5,6,7,8,9, 10,12, 14,20,25,30,40,50)
  library(furrr)
  
  
  for(i in 10:20){
    
    draws_id = sample(nrow(post_draws), 1)
    
    print(draws_id)
    
    new_subjects =
      post_draws %>% filter(draw %in% draws_id) %>%
      mutate(`gm.8.` = 0.05) %>%
      rowwise() %>% 
      mutate(sim = list(simulate_subjects(cur_data(), S = max_subjects,P = 9))) %>% 
      unnest(sim) %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(sim_d = list(sim_trials(cur_data())))
    
    safe_function <- possibly(
      function(n) fitter(new_subjects, n),
      otherwise = "Error"
    )
    
    results_list <- future_map(
      n_subj,
      ~ safe_function(.x),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
    
    saveRDS(results_list,here::here(paste0("test_effect_0.05_both_ss",i,".RData")))
    
    
  }

  
  files = list.files((here::here("Simulations","Heurestic","Sequential Sampling","results")), full.names = T)
  files = files[3:5]
  results_list = list()
  q = 0
  for(i in 1:length(files)){
    q = q+1
    results_list[[q]] = readRDS(files[i])
  }
  
  
  # results_list = results_list[-which(results_list == "Error")]
  
  purrr::map_dfr(
    results_list,
    function(x) {
      valid <- purrr::keep(x, ~ !identical(.x, "Error"))
      dplyr::bind_rows(purrr::map(valid, purrr::pluck, 1))
    }
  ) %>% 
    # filter(div == 0) %>%
    # filter(draw_id == 2100 |draw_id == 1535) %>% 
    # mutate(draw_id = if_else(row_number() <= 6, draw_id + 1, draw_id)) %>% 
    pivot_longer(cols = c("bf8","bf9")) %>% 
    ggplot(aes(x = n_subs, y =value, group = draw_id))+
    geom_line()+
    geom_point()+
    facet_wrap(~name)+
    geom_hline(yintercept = c(1/30), linetype = 2)+
    # scale_y_continuous(limits = c(0,1))+
    theme_classic()
  
  
  purrr::map_dfr(
    results_list[1:length(results_list)],
    ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 2))
  ) %>%   
    # mutate(draw_id = if_else(row_number() <= 6, draw_id + 1, draw_id)) %>% 
    ggplot(aes(x = simulated, y = mean,ymin = q5,ymax = q95, group = draw_id))+
    geom_pointrange()+
    facet_grid(variable~n_subs, scales = "free")+
    theme_classic()+
    geom_abline(col = "red")
  
  
  
  purrr::map_dfr(
    results_list[1:length(results_list)],
    ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 3))
  ) %>%   
    filter(grepl("gm",variable)) %>% 
    filter(variable %in% c("gm[9]","gm[8]")) %>% 
    filter(div == 0 & tree == 0) %>% 
    # mutate(draw_id = if_else(row_number() <= 6, draw_id + 1, draw_id)) %>% 
    ggplot(aes(x = simulated, y = mean,ymin = q5,ymax = q95, group = draw_id))+
    geom_pointrange()+
    facet_grid(variable~n_subs, scales = "free")+
    theme_classic()+
    geom_abline(col = "red")
  

  
  purrr::map_dfr(
    results_list[1:length(results_list)],
    ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 3))
  ) %>%   
    filter(grepl("tau_u",variable)) %>% 
    # mutate(draw_id = if_else(row_number() <= 6, draw_id + 1, draw_id)) %>% 
    ggplot(aes(x = simulated, y = mean,ymin = q5,ymax = q95, group = draw_id))+
    geom_pointrange()+
    facet_grid(variable~n_subs, scales = "free")+
    theme_classic()+
    geom_abline(col = "red")
  
  
  
  
    
}
