library(future)
library(future.apply)
library(ordbetareg)
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
  offsets  = c(-5,-4,-3,-2,-1, 1, 2, 3, 4, 5)
  w = c(1, 2, 3, 2, 1, 1, 2, 3, 2, 1)
  
  
  stim_values_per_block =
    rep((brms::inv_logit_scaled(
      (brms::inv_logit_scaled(df$alpha1) - 0.5) * 2 +
        rep(offsets / exp(df$beta1), times = w)
    ) - 0.5) * 2,3)
  
  X = rep(stim_values_per_block,6)
  
  
  # -----------------------------
  # Decision
  # -----------------------------
  theta = brms::inv_logit_scaled(df$lapse)/2 + (1-2*brms::inv_logit_scaled(df$lapse)/2) * pnorm(exp(df$beta1) * (X - (brms::inv_logit_scaled(df$alpha1) - 0.5) * 2))
  
  resp = rbinom(length(X), 1, theta)        # Simulated response
  
  # Accuracy based on simulated response
  ACC = ifelse(X < 0 & resp == 0, 1,
               ifelse(X > 0 & resp == 1, 1, 0))
  
  # -----------------------------
  # Simulated confidence
  # -----------------------------
  conf_mu = numeric(length(X))
  for(i in 1:length(X)){
    if(ACC[i] == 1){
      conf_mu[i] <- pnorm(exp(df$beta1 - exp(df$meta_un_cor1)) * abs(X[i] - (brms::inv_logit_scaled(df$alpha1) - 0.5) * 2))       # Confidence for correct trials
    } else {
      conf_mu[i] <- pnorm(df$meta_un_inc1 * abs(X[i] - (brms::inv_logit_scaled(df$alpha1) - 0.5) * 2))    # Confidence for incorrect trials
    }
  }
  
  conf = rordbeta(n = length(X), mu = conf_mu, phi = exp(df$conf_prec1), cutpoints = c(c0, c1))  # Simulated confidence values
  
  # -----------------------------
  # Dataframe
  # -----------------------------
  df = data.frame(
    c0 = c0,
    c1 = c1,
    beta = df$beta1,
    lapse = df$lapse,
    alpha = df$alpha1, 
    meta_un_cor = df$meta_un_cor1,
    meta_un_inc1 = df$meta_un_inc1,
    conf_prec1 = df$conf_prec1,
    X = X,
    interval = runif(length(X), 1,5),
    ACC = ACC,
    resp = resp,
    c_mu = conf_mu,
    C = conf,
    theta = theta
  )
  
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
      alpha1         = param_new[1],
      beta1          = param_new[2],
      conf_prec1     = param_new[3],
      meta_un_cor1   = param_new[4],
      meta_un_inc1   = param_new[5],
      meta_bias      = param_new[6],
      lapse          = param_new[7],
      meta_un_beta   = param_new[8],
      meta_bias_beta = param_new[9]
    )
  })
}

fitter = function(new_subjects,subs){
  
  sim_id = rnorm(1,0,100)
  
  dd = new_subjects %>% select(sim_d,subj_id) %>% unnest() %>% filter(subj_id %in% 1:subs)
  
  model = cmdstanr::cmdstan_model(here::here("Stanmodels","Heurestic models",
                                             "Emperical dataanalysis","hierarchical","lapse","full_model3_lapse_optimized.stan"))
  
  
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
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 10,
    iter_sampling = 1000,
    iter_warmup = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    # init  = 0,
    parallel_chains = 4)
  
  prior_sd = 0.5
  
  divs = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent),
                                                            mean_tree = mean(num_max_treedepth),
                                                            mean_energi = mean(ebfmi))
  
  bf8 = get_bf(post = as.numeric(as_draws_df(fit$draws("gm[8]")) %>% .$`gm[8]`),prior = rnorm(4000,0,prior_sd))
  bf9 = get_bf(post = as.numeric(as_draws_df(fit$draws("gm[9]")) %>% .$`gm[9]`),prior = rnorm(4000,0,prior_sd))
  
  BF = data.frame(bf8 = bf8, bf9 = bf9, n_subs = subs,sim_id = sim_id,
                  div = divs$mean_div, tree = divs$mean_tree, energi = divs$mean_energi, draw_id = unique(new_subjects$draw),prior_sd = prior_sd)
  
  sim_subj = dd %>% select(c0,c1,beta,lapse,alpha,meta_un_cor,meta_un_inc1,conf_prec1,subj_id) %>% distinct() %>% 
    pivot_longer(-subj_id, names_to = "variable",values_to = "simulated")
  
  esti_subj = fit$summary(c("c0","c1","beta1","lapse","alpha1","meta_un_cor1","meta_un_inc1","conf_prec1")) %>% 
    mutate(
      subj_id = str_extract(variable, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
      variable   = str_remove(variable, "\\[\\d+\\]")
    ) %>% 
    mutate(variable = ifelse(variable == "alpha1","alpha",ifelse(variable == "meta_un_cor1","meta_un_cor",ifelse(variable == "beta1","beta",variable))))
  
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
  
  
  
  
  plan(multisession, workers = 6)  # Windows-friendly
  safe_function <- possibly(
    function(n) fitter(new_subjects, n),
    otherwise = "Error"
  )
  
  for(i in 1:5){
    n_subj <- c(5, 10, 15,20,30,50)
    
    library(furrr)
    results_list <- future_map(
      n_subj,
      ~ safe_function(.x),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
    
    saveRDS(results_list,here::here(paste0("test_ss",i,".RData")))
    
    
  }

  results_list = results_list[-which(results_list == "Error")]
  
  files = list.files((here::here("Simulations","Heurestic","Sequential Sampling","results")), full.names = T)
  
  results_list = list()
  
  for(i in 1:length(files)){
    results_list[[i]] = readRDS(files[i])
  }
  
  purrr::map_dfr(
    results_list[1:length(results_list)],
    ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 1))
  ) %>%   
    # mutate(draw_id = if_else(row_number() <= 6, draw_id + 1, draw_id)) %>% 
    pivot_longer(cols = c("bf8","bf9")) %>% 
    ggplot(aes(x = n_subs, y =value, group = draw_id))+geom_line()+facet_wrap(~name)+
    geom_hline(yintercept = c(1/30), linetype = 2)+
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
