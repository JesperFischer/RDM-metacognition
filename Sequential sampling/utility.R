
load_data = function(){
  
  
  paths = list.files(here::here("Sequential sampling","Data"), recursive = T, full.names = T)[grepl("\\d\\.csv$",list.files(here::here("Sequential sampling","Data"), recursive = T, full.names = T))]
  
  df = data.frame()
  for(path in paths){
    print(path)
    subjectID = strsplit(path, '/')[[1]][[9]]
    
    
    dd = read.csv(path) %>% mutate(X = NULL) %>% 
      mutate(ID = subjectID) %>% 
      mutate(resp = ifelse(resp == "down",0,1),
             interval = as.numeric(interval),
             D = ifelse(dots.direction == 270,-1,1),
             stim = ifelse(dots.direction == 270, -coherence,coherence)) %>%
      mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
                          ifelse(resp == 0 & stim < 0, 1,0))) %>% filter(Trialtype == "Main") %>% 
      filter(SR_conf != "None") %>% filter(scale == "conf")
    
    if ("sub" %in% names(dd)) {
      names(dd)[names(dd) == "sub"] = "subject"
    }
    
    dd = dd %>% mutate(SR_conf = as.numeric(SR_conf)) %>%
      rename(X = stim, RT = RTdec, Correct = cor, Y = resp, Confidence = SR_conf) %>% mutate(subject = as.numeric(sub(".*sub(\\d{3})\\.csv$", "\\1", path))) %>% 
      select(X,Correct,Y,RT,Confidence,scale,D, coherence,Trialtype,ACC,RTrating,
             interval,interTrial.interval,
             subject,age,gender,ID)
    
    
    df = rbind(df,dd)
  }
  
  return(df)
  
}



plot_beh_data = function(df,n_bins,ACC = F, r_data = F){
  subjects = unique(df$subject)
  IDs = unique(df$ID)
  # If more than one subject → apply function per subject
  if(length(subjects) != 1){
    
    df_split = split(df, df$ID)
    
    out = lapply(df_split, function(d){
      plot_beh_data(d, n_bins = n_bins, ACC = ACC, r_data  = r_data)
    })
    if(is.ggplot(out[[1]])){
      for(i in 1:length(out)){
        out[[i]] = out[[i]] + ggtitle(unique(df_split[[i]]$ID))
      }      
    }
    return(out)
  }
  
  
  if(r_data == T){
    return(df %>% mutate(trial = 1:n()))
  }
  
  
  # Prepare observed data
  if (!is.null(n_bins)) {
    # Create common bin boundaries based on the range of both datasets
    all_X <- c(df$X)
    X_range <- range(all_X, na.rm = TRUE)
    bin_breaks <- seq(X_range[1], X_range[2], length.out = n_bins + 1)
    
    # Calculate bin centers (midpoints)
    bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2
    
    # Bin the data using the common breaks and assign bin centers
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
             X = bin_centers[X_bin]) %>%
      select(-X_bin)%>% mutate(trial = 1:n())
    
  }
  
  if(ACC){
    bin = df %>% mutate(trial = 1:n())%>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X,ID) %>%
      summarize(
        name = "Type-1",
        mean = mean(ACC, na.rm = T),
        q5 = mean(ACC) - 2* (mean(ACC, na.rm = T) * (1-mean(ACC, na.rm = T)) / sqrt(n())),
        q95 = mean(ACC) + 2 * (mean(ACC, na.rm = T) * (1-mean(ACC, na.rm = T)) / sqrt(n())),
        .groups = "drop"
      )
  }else{
    bin = df %>% 
      mutate(trial = 1:n()) %>% 
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X,ID) %>%
      summarize(
        name = "Type-1",
        
        k = sum(Y),
        n = n(),
        mean = k / n,
        
        a_post = k + 1,
        b_post = (n - k) + 1,
        
        q5 = qbeta(0.025, a_post, b_post),
        q95 = qbeta(0.975, a_post, b_post),
        .groups = "drop"
      )
  }
  
  
  
  
  df1 = bind_rows(
    bin,
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT, na.rm = T),
                q5 = mean(RT, na.rm = T) - 2 * (sd(RT, na.rm = T) / sqrt(n())),
                q95 = mean(RT, na.rm = T) + 2 * (sd(RT, na.rm = T) / sqrt(n())),
                .groups = "drop"),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct,ID) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence, na.rm = T),
                q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                .groups = "drop")
  ) 
  
  
  # Plot 1: Expected means (main plot)
  plot_mean = df1 %>% 
    ggplot() +
    geom_pointrange(data = df1, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                    shape = 21, color = "black", alpha = 0.5) +
    (if (!is.null(n_bins)) 
      geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1)
     else NULL) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")
  
  
  
  
  return(plot_mean)
  
  
}

plot_beh_data_tog = function(df,n_bins,ACC = F, r_data = F){
  subjects = unique(df$subject)
  IDs = unique(df$ID)
  # If more than one subject → apply function per subject
  if(length(subjects) != 1){
    
    df_split = split(df, df$subject)
    
    out = lapply(df_split, function(d){
      plot_beh_data_tog(d, n_bins = n_bins, ACC = ACC, r_data  = r_data)
    })
    
    if(is.ggplot(out[[1]][[1]])){
      for(i in 1:length(out)){
        out[[i]][[1]] = out[[i]][[1]] + ggtitle(unique(df_split[[i]]$ID))
        out[[i]][[2]] = out[[i]][[2]] + ggtitle(unique(df_split[[i]]$ID))
      }      
    }
    return(out)
  }
  
  
  if(r_data == T){
    return(df %>% mutate(trial = 1:n()))
  }
  
  
  # Prepare observed data
  if (!is.null(n_bins)) {
    # Create common bin boundaries based on the range of both datasets
    all_X <- c(df$X)
    X_range <- range(all_X, na.rm = TRUE)
    bin_breaks <- seq(X_range[1], X_range[2], length.out = n_bins + 1)
    
    # Calculate bin centers (midpoints)
    bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2
    
    # Bin the data using the common breaks and assign bin centers
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
             X = bin_centers[X_bin]) %>%
      select(-X_bin)%>% mutate(trial = 1:n())
    
  }
  
  if(ACC){
    bin = df %>% mutate(trial = 1:n())%>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X,ID) %>%
      summarize(
        name = "Type-1",
        mean = mean(ACC, na.rm = T),
        q5 = mean(ACC) - 2* (mean(ACC, na.rm = T) * (1-mean(ACC, na.rm = T)) / sqrt(n())),
        q95 = mean(ACC) + 2 * (mean(ACC, na.rm = T) * (1-mean(ACC, na.rm = T)) / sqrt(n())),
        .groups = "drop"
      )
  }else{
    bin = df %>% 
      mutate(trial = 1:n()) %>% 
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X,ID) %>%
      summarize(
        name = "Type-1",
        
        k = sum(Y),
        n = n(),
        mean = k / n,
        
        a_post = k + 1,
        b_post = (n - k) + 1,
        
        q5 = qbeta(0.025, a_post, b_post),
        q95 = qbeta(0.975, a_post, b_post),
        .groups = "drop"
      )
  }
  
  
  
  
  df1 = bind_rows(
    bin,
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT, na.rm = T),
                q5 = mean(RT, na.rm = T) - 2 * (sd(RT, na.rm = T) / sqrt(n())),
                q95 = mean(RT, na.rm = T) + 2 * (sd(RT, na.rm = T) / sqrt(n())),
                .groups = "drop"),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct,ID) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence, na.rm = T),
                q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                .groups = "drop")
  ) 
  
  
  # Plot 1: Expected means (main plot)
  plot_mean_ACC = data.frame() %>% 
    ggplot() +
    geom_pointrange(data = df1 %>% filter(name == "Confidence"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                    shape = 21, color = "black", size = 0.8) +
    geom_pointrange(data = df1 %>% filter(name == "Type-1"), aes(x = X, y = mean, ymin = q5, ymax = q95), size = 0.8,
                    shape = 21, color = "black", fill = "black") +
    
    (if (!is.null(n_bins)) 
      geom_line(data = df1 %>% filter(name == "Confidence"), aes(x = X, y = mean, color = Correct),
                size = 0.8) 
     else NULL) +
    (if (!is.null(n_bins)) 
      geom_line(data = df1 %>% filter(name == "Type-1"), aes(x = X, y = mean), size = 0.8)
     else NULL) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    scale_color_manual(values = c("red","green","black"))+
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")
  
  
  
  df1 = bind_rows(
    bin,
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT, na.rm = T),
                q5 = mean(RT, na.rm = T) - 2 * (sd(RT, na.rm = T) / sqrt(n())),
                q95 = mean(RT, na.rm = T) + 2 * (sd(RT, na.rm = T) / sqrt(n())),
                .groups = "drop"),
    
    df %>%
      mutate(Y = ifelse(Y == 1, "Right", "Left")) %>%
      group_by(X, Y,ID) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence, na.rm = T),
                q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                .groups = "drop")
  ) 
  
  
  plot_mean_dec = data.frame() %>% 
    ggplot() +
    geom_pointrange(data = df1 %>% filter(name == "Confidence"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Y),
                    shape = 21, color = "black", size = 0.8) +
    geom_pointrange(data = df1 %>% filter(name == "Type-1"), aes(x = X, y = mean, ymin = q5, ymax = q95), size = 0.8,
                    shape = 21, color = "black", fill = "black") +
    
    (if (!is.null(n_bins)) 
      geom_line(data = df1 %>% filter(name == "Confidence"), aes(x = X, y = mean, color = Y),
                size = 0.8) 
     else NULL) +
    (if (!is.null(n_bins)) 
      geom_line(data = df1 %>% filter(name == "Type-1"), aes(x = X, y = mean), size = 0.8)
     else NULL) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    scale_color_manual(values = c("red","green","black"))+
    theme_classic(base_size = 14) +
    labs(color = "Y", fill = "Y",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")
  
  
  return(list(plot_mean_ACC,plot_mean_dec))
  
  
}



fit_model_ss = function(df,
                        model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Single Subject.stan")),
                        samples = 1000){
  
  
  datastan = list(N = nrow(df),
                  a = df$Y,
                  RT = df$RT,
                  D = df$D,
                  C = df$Confidence,
                  X = df$X * df$D,
                  XD = df$X,
                  minRT = min(df$RT),
                  ACC = df$Correct,
                  starts = 1,
                  ends = nrow(df),
                  interval = df$interval)
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 100,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.95,
    max_treedepth = 12,
    parallel_chains = 4)
  
  if(grepl("Single Subject",model$stan_file())){
    fit$save_object(here::here("Sequential sampling","Fits","Single Subject",paste0("fit_",unique(df$ID),".rds")))
  }
  
  return(fit)
  
}


fit_model_ss_nolapse = function(df,
                        model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Single Subject nolapse.stan")),
                        samples = 1000){
  
  
  datastan = list(N = nrow(df),
                   a = df$Y,
                  RT = df$RT,
                  D = df$D,
                  C = df$Confidence,
                  X = df$X * df$D,
                  XD = df$X,
                  minRT = min(df$RT),
                  ACC = df$Correct,
                  starts = 1,
                  ends = nrow(df),
                  interval = df$interval)
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 100,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.95,
    max_treedepth = 12,
    parallel_chains = 4)
  
  
  if(grepl("Single Subject",model$stan_file())){
    fit$save_object(here::here("Sequential sampling","Fits","Single Subject nolapse",paste0("fit_",unique(df$ID),".rds")))
  }
  
  return(fit)
  
}




dia = function(fit,df){
  
  if(length(fit) > 1 & is.list(fit)){
    out <- lapply(names(fit), function(nm) {
      dia(fit[[nm]], df %>% filter(ID == nm))
    })
    
    return(out)
    
  }
  
  available <- names(as_draws_df(fit$draws()))

  parameters = c("c0","c11","beta","sigma_e","sigma_k","sigma_m","confprec_p","meta_bias")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
  }
  
  
  sum = fit$summary(params) %>% mutate(ID = unique(df$ID)) %>% dplyr::select(-c(mad,median))
  
  plot = as_ggplot(mcmc_pairs(fit$draws(c(params)),np = nuts_params(fit)))+
    plot_annotation(
      title = paste(unique(df$ID), collapse = ", ")
    )
  
  plot2 = mcmc_trace(fit$draws(c(params)),np = nuts_params(fit))+ ggtitle(unique(df$ID))
  
  return(list(sum,plot,plot2, fit$diagnostic_summary()))
  
}


pp = function(fit,df,n_bins){
  
  if(length(fit) > 1 & is.list(fit)){
    out <- lapply(names(fit), function(nm) {
      pp(fit[[nm]], df %>% filter(ID == nm),n_bins)
    })
    
    return(out)
    
  }
  
  available <- names(as_draws_df(fit$draws()))
  
  parameters = c("c0","c11","beta","sigma_e","sigma_k","sigma_m","confprec","meta_bias","lapse")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    behplot = plot_beh_data(df,n_bins,F, r_data = T)
    # Prepare observed data
    if (!is.null(n_bins)) {
      # Create common bin boundaries based on the range of both datasets
      behplot <- behplot %>%
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
    }
    
    
    df_behplot = bind_rows(
      behplot %>% 
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X,ID) %>%
        summarize(
          name = "Type-1",
          
          k = sum(Y),
          n = n(),
          mean = k / n,
          
          a_post = k + 1,
          b_post = (n - k) + 1,
          
          q5 = qbeta(0.025, a_post, b_post),
          q95 = qbeta(0.975, a_post, b_post),
          .groups = "drop"
        ),
      behplot %>%
        mutate(action = ifelse(Y == 1, "up", "down")) %>%
        group_by(X, action,ID) %>%
        summarize(name = "Confidence",
                  mean = mean(Confidence, na.rm = T),
                  q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                  q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                  .groups = "drop")
    ) 
    
    
    
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
                    pnorm((1/(sqrt(1 + (Cc^2*abs(X)^2) / sigma2))) * (e_cond * abs(X) * Cc) / sigma2),
                    1-pnorm((1/(sqrt(1 + (Cc^2*abs(X)^2) / sigma2))) * (e_cond * abs(X) * Cc) / sigma2))
      
      return(c_mu)
      
    }
    
    pred_data = as_draws_df(fit$draws(parameters)) %>% select(-contains(".")) %>% 
      # rename_with(~c("mu","prec","lapse")) %>% 
      mutate(draw = 1:n()) %>% 
      mutate(x = list(seq(-1,1,by = 0.05))) %>% 
      unnest() %>% 
      mutate(interval = mean(df$interval)) %>% 
      group_by(draw) %>% 
      mutate(p = psycho(x, beta,sigma_e,sigma_k)) %>% 
      rowwise() %>% 
      mutate(action = rbinom(1,1,p)) %>% 
      unnest() %>% 
      rowwise() %>% 
      mutate(conf_mu = ifelse(action == 1,
                              confs(1,x,beta, sigma_e, sigma_k,sigma_m),
                              confs(0,x,beta, sigma_e, sigma_k,sigma_m))) %>% 
      mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias)) %>% 
      ungroup() %>% 
      mutate(ID = unique(behplot$ID)[1])
    
    
    bin = pred_data %>%
      mutate(action = ifelse(action == 1, "up", "down")) %>%
      group_by(ID,x) %>%
      summarize(
        name = "Type-1",
        mean = mean(p, na.rm = T),
        q5 = quantile(p, 0.05),
        q95 = quantile(p, 0.95),
        .groups = "drop"
      )
    
    df1 = bind_rows(
      bin,
      pred_data %>%
        mutate(action = ifelse(action == 1, "up", "down")) %>%
        group_by(x, action,ID) %>%
        summarize(name = "Confidence",
                  mean = mean(conf_mu, na.rm = T),
                  q5 = quantile(conf_mu, 0.05, na.rm = T),
                  q95 = quantile(conf_mu, 0.95, na.rm = T),
                  .groups = "drop")
    ) 
    
    
    
    pp_plot = df_behplot %>% 
      ggplot() +
      geom_pointrange(data = df_behplot, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
                      shape = 21, color = "black", alpha = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
      theme_classic(base_size = 14) +
      labs(color = "action", fill = "action",
           y = "Value") +
      geom_vline(xintercept = 0, linetype = 2) +
      labs(subtitle = unique(df_behplot$ID)[1])+
      scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
      # theme(legend.position = "top")+
      geom_line(data = df1, aes(x = x, y = mean, col = action))+
      geom_ribbon(data = df1, aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)
    
    
    pp_plot
    
    
    
    
    return(pp_plot)
    
  }
  
  
}




fit_model_hier = function(dd,
                        model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Hierarchical.stan")),
                        samples = 1000){
  
  
  
  datastan = list(N = nrow(dd),
                  S = length(unique(dd$ID)),
                  S_id = as.numeric(as.factor(dd$ID)),
                  a = dd$Y,
                  a_vec = dd$Y,
                  RT = dd$RT,
                  D = dd$D,
                  C = dd$Confidence,
                  X = dd$X * dd$D,
                  XD = dd$X,
                  minRT = min(dd$RT),
                  ACC = dd$Correct,
                  interval = dd$interval)
  
  pf = model$pathfinder(data = datastan, psis_resample = F, calculate_lp = F)
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 100,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.95,
    init = pf,
    max_treedepth = 12,
    parallel_chains = 4)
  
  fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical.rds"))
  
  
  if(grepl("Hierarchical",model$stan_file())){
    
    if(dir.exists(here::here("Sequential sampling","Fits","Hierarchical"))){
      fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical.rds"))
    }else{
      dir.create(here::here("Sequential sampling","Fits","Hierarchical"))
      fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical.rds"))
      
    }
    
  }
  
  return(fit)
  
}


fit_model_hier_nolapse = function(dd,
                          model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Hierarchical_nolapse.stan")),
                          samples = 1000){
  model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Hierarchical_nolapse_clamped.stan"))
  
  samples = 1000
  
  datastan = list(N = nrow(dd),
                  S = length(unique(dd$ID)),
                  S_id = as.numeric(as.factor(dd$ID)),
                  a = dd$Y,
                  a_vec = dd$Y,
                  RT = dd$RT,
                  D = dd$D,
                  C = dd$Confidence,
                  X = dd$X * dd$D,
                  XD = dd$X,
                  minRT = min(dd$RT),
                  ACC = dd$Correct,
                  interval = dd$interval)
  
  pf = model$pathfinder(data = datastan, psis_resample = F, calculate_lp = F)
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 100,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.99,
    init = pf,
    max_treedepth = 12,
    parallel_chains = 4)
  
  fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical_nolapse_39_clamped.rds"))
  
  if(grepl("Hierarchical",model$stan_file())){
    
    if(dir.exists(here::here("Sequential sampling","Fits","Hierarchical"))){
      fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical_nolapse.rds"))
    }else{
      dir.create(here::here("Sequential sampling","Fits","Hierarchical"))
      fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical_nolapse.rds"))
      
    }
    
  }
  
  return(fit)
  
}



dia_hier = function(fit,df,lapse){
  
  
  available <- names(as_draws_df(fit$draws()))
 
  if(lapse){
    parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","lapse[1]","confprec[1]","meta_bias[1]","sigmam_beta[1]","meta_bias_beta[1]")
    
  }else{
    parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","sigmam_beta[1]","meta_bias_beta[1]")
  }
  
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    n_subj = length(unique(df$ID))
    subs = rep(unique(df$ID), length(params))
    
    subj_parameters = str_sub(parameters, 1, -4)
    
    subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
  }
  
  
  
  
  sum = fit$summary(c("gm","tau_u")) %>% dplyr::select(-c(mad,median))
  sub = fit$summary(c(subj_parameters)) %>% dplyr::select(-c(mad,median)) %>% mutate(ID = subs)
  
  
  
  plot = as_ggplot(mcmc_pairs(fit$draws(c(c("gm","tau_u"))),np = nuts_params(fit)))+
    plot_annotation(
      title = paste(unique(df$ID), collapse = ", ")
    )
  
  plot2 = mcmc_trace(fit$draws(c(c("gm","tau_u"))),np = nuts_params(fit))+ ggtitle(unique(df$ID))
  
  plot3 = mcmc_trace(fit$draws(subj_parameters),np = nuts_params(fit))
  
  return(list(sum,sub,plot,plot2,plot3, fit$diagnostic_summary()))
  
}


pp_hier_nolapse = function(fit,df,n_bins){
  
  
  available <- names(as_draws_df(fit$draws()))
  
  
  parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","sigmam_beta[1]","meta_bias_beta[1]")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
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
      mutate(x = list(seq(-1,1,by = 0.1))) %>% 
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
      labs(color = "action", fill = "action",
           y = "P(a=1) Or P(a=D)", x = "Signed Stimulus (XD)") +
      geom_vline(xintercept = 0, linetype = 2) +
      # labs(subtitle = name)+
      facet_wrap(~variable, ncol = 1, scales = "free")
    
    
    
    
    return(list(group_plot,sub_plots,average_sub_average_group))
    
    
  }
  
  
}

pp_hier_nolapse_centbeta = function(fit,df,n_bins){
  
  
  # available <- names(as_draws_df(fit$draws()))
  
  
  parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","sigmam_beta[1]","meta_bias_beta[1]")
  # if(length(intersect(parameters, available)) > 5){
  # params <- intersect(parameters, available)
  
  params = parameters
  
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
    rename_with(~c(  "sigma_e",
                     "sigma_k",
                     "sigma_m",
                     "meta_bias",
                     "confprec",
                     "sigmam_beta",
                     "meta_bias_beta")) %>%
    bind_cols(as_draws_df(fit$draws("gm_beta")) %>% select(-contains(".")) %>% rename(beta = `gm_beta`)) %>%
    mutate(draw = 1:n()) %>% 
    mutate(x = list(seq(-1,1,by = 0.1))) %>% 
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
  ) %>% rename(variable = name)
  
  
  group_plot = behplot %>% 
    #filter(name != "RT") %>% 
    ggplot() +
    geom_pointrange(data = behplot, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
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
  
  # df_s = df
  
  
  behplot = rbind(
    bin = df_s %>%
      mutate(action = NA) %>%
      group_by(ID,X,action) %>%
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
      group_by(ID,X, action) %>%
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
    labs(color = "action", fill = "action",
         y = "P(a=1) Or P(a=D)", x = "Signed Stimulus (XD)") +
    geom_vline(xintercept = 0, linetype = 2) +
    # labs(subtitle = name)+
    facet_wrap(~variable, ncol = 1, scales = "free")
  
  
  
  
  return(list(group_plot,sub_plots,average_sub_average_group))
  
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
