fit_model_ss = function(df, model,samples){
  
  

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
                interval = df$interTrial.interval)


  fit <-model$sample(
    data = datastan,
    refresh = 300,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.95,
    max_treedepth = 12,
    init  = 0,
    parallel_chains = 4)
  
  if(grepl("MCMC",model1$stan_file())){
    fit$save_object(here::here("Experimental analysis","MCMC","saved models",paste0("fit_",unique(df$ID),".rds")))
  }
  if(grepl("CANSANDRE",model1$stan_file())){
    fit$save_object(here::here("Experimental analysis","CANSANDRE","saved models",paste0("fit_",unique(df$ID),".rds")))
  }
  if(grepl("Heurestic models",model1$stan_file())){
    fit$save_object(here::here("Experimental analysis","Heurestic models","saved models",paste0("fit_",unique(df$ID),".rds")))
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
  parameters = c("sigma_choice","mean_choice","prec_conf","sigma_m","sigma_e")
  
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
  }
  parameters = c("c0","c11","alpha1","beta1","alpha","beta","conf_prec1","meta_un_cor1","meta_un_inc1")
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
  parameters = c("sigma_choice","mean_choice","prec_conf","sigma_m","sigma_e")
  
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
  }
  parameters = c("c0","c11","alpha1","beta1","conf_prec1","meta_un_cor1","meta_un_inc1")
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
    
    behplot = plot_beh_data(df,n_bins,F, r_data = T)
    
    psycho = function(x,alpha,beta){
      return(pnorm(beta * (x-alpha)))
    }
    
    psycho_ACC = function(x,alpha,beta){
      return(pnorm(beta * abs(x-alpha)))
    }
    
    pred_data = as_draws_df(fit$draws(parameters)) %>% select(-contains(".")) %>% 
      # rename_with(~c("mu","prec","lapse")) %>% 
      mutate(draw = 1:n()) %>% 
      mutate(x = list(seq(min(df$X)-0.2,max(df$X)+0.2,by = 0.05))) %>% unnest() %>% 
      group_by(draw) %>% 
      mutate(p = psycho(x, (brms::inv_logit_scaled(alpha1)-0.5)*2,exp(beta1)),
             resp = rbinom(n(),1,p)) %>% 
      mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
      mutate(conf_mu = ifelse(ACC == 1, psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1) + meta_un_cor1), psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1) + meta_un_inc1))) %>% 
      mutate(ID = unique(behplot$ID)[1])
    
    
    bin = pred_data %>%
      mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
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
      mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
      group_by(x, Correct,ID) %>%
      summarize(name = "Confidence",
                mean = mean(conf_mu, na.rm = T),
                q5 = quantile(conf_mu, 0.05),
                q95 = quantile(conf_mu, 0.95),
                .groups = "drop")
  ) 
    
  
  
  pp_plot = behplot %>% filter(name != "RT") %>% 
    ggplot() +
    geom_pointrange(data = behplot%>% filter(name != "RT"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                    shape = 21, color = "black", alpha = 0.5) +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(subtitle = unique(behplot$ID)[1])+
    scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
    # theme(legend.position = "top")+
    geom_line(data = df1, aes(x = x, y = mean, col = Correct))+
    geom_ribbon(data = df1, aes(x = x, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.5)
  
  
  pp_plot
  
  
  
  
  return(pp_plot)
  
}

}



dia_hier = function(fit,df){
  

  available <- names(as_draws_df(fit$draws()))
  parameters = c("sigma_choice[1]","mean_choice[1]","prec_conf[1]","sigma_m[1]","sigma_e[1]")
  
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
  }
  parameters = c("c0[1]","c11[1]","alpha1[1]","beta1[1]","conf_prec1[1]","meta_un_cor1[1]","meta_un_inc1[1]")
  if(length(intersect(parameters, available)) != 0){
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


pp_hier = function(fit,df,n_bins){
  

  
  available <- names(as_draws_df(fit$draws()))
  parameters = c("sigma_choice[1]","mean_choice[1]","prec_conf[1]","sigma_m[1]","sigma_e[1]")
  
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
  }
  
  parameters = c("c0[1]","c11[1]","alpha1[1]","beta1[1]","conf_prec1[1]","meta_un_cor1[1]","meta_un_inc1[1]")
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
    
    n_subj = length(unique(df$ID))
    subs = rep(unique(df$ID), length(params))
    
    subj_parameters = str_sub(parameters, 1, -4)
    
    subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
  
    
    psycho = function(x,alpha,beta){
      return(pnorm(beta * (x-alpha)))
    }
    
    psycho_ACC = function(x,alpha,beta){
      return(pnorm(beta * abs(x-alpha)))
    }
    
    sum = fit$summary(c("gm","tau_u")) %>% dplyr::select(-c(mad,median))
    sub = fit$summary(c(subj_parameters)) %>% dplyr::select(-c(mad,median)) %>% mutate(ID = subs)
    
    pred_data_group = as_draws_df(fit$draws("gm")) %>% select(-contains(".")) %>% 
      rename_with(~c(  "alpha1",
                       "beta1",
                       "conf_prec1",
                       "meta_un_cor1",
                       "meta_un_inc1")) %>%
      mutate(draw = 1:n()) %>% 
      mutate(x = list(seq(min(df$X)-0.2,max(df$X)+0.2,by = 0.05))) %>% unnest() %>% 
      group_by(draw) %>% 
      mutate(p = psycho(x, (brms::inv_logit_scaled(alpha1)-0.5)*2,exp(beta1)),
             resp = rbinom(n(),1,p)) %>% 
      mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
      mutate(conf_mu = ifelse(ACC == 1, psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1) + meta_un_cor1), psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1) + meta_un_inc1)))
    
    
    df1 = bind_rows(
      pred_data_group %>%
        mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
        group_by(x) %>%
        summarize(
          name = "Type-1",
          mean = mean(p, na.rm = T),
          q5 = quantile(p, 0.05),
          q95 = quantile(p, 0.95),
          .groups = "drop"
        ),
      pred_data_group %>%
        mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
        group_by(x, Correct) %>%
        summarize(name = "Confidence",
                  mean = mean(conf_mu, na.rm = T),
                  q5 = quantile(conf_mu, 0.05),
                  q95 = quantile(conf_mu, 0.95),
                  .groups = "drop")
    ) 
    
    
    
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
        select(-X_bin)
      
    }
    
    
    behplot = rbind(
      bin = df %>%
        mutate(Correct = NA) %>%
        group_by(X,Correct) %>%
        summarize(
          name = "Type-1",
          mean = mean(Y, na.rm = T),
          q5 = mean(Y, na.rm = T) - 2* (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
          q95 = mean(Y, na.rm = T) + 2 * (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
          .groups = "drop"
        ),
        df %>%
          mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
          group_by(X, Correct) %>%
          summarize(name = "Confidence",
                    mean = mean(Confidence, na.rm = T),
                    q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                    q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                    .groups = "drop")
      ) 
      
      
    group_plot = behplot %>% filter(name != "RT") %>% 
      ggplot() +
      geom_pointrange(data = behplot%>% filter(name != "RT"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                      shape = 21, color = "black", alpha = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
      theme_classic(base_size = 14) +
      labs(color = "Correct", fill = "Correct",
           y = "Value") +
      geom_vline(xintercept = 0, linetype = 2) +
      scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
      # theme(legend.position = "top")+
      geom_line(data = df1, aes(x = x, y = mean, col = Correct))+
      geom_ribbon(data = df1, aes(x = x, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.5)
    
    
    
    pred_subj = as_draws_df(fit$draws(c("alpha1","beta1","meta_un_cor1","meta_un_inc1"))) %>% select(-contains(".")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw) %>% 
      mutate(
        sub_idx = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
        param   = str_remove(name, "\\[\\d+\\]"),
        ID  = subs[sub_idx],
        name = NULL
      ) %>% 
      pivot_wider(names_from = "param", values_from = "value") %>% 
      mutate(x = list(seq(min(df$X)-0.2,max(df$X)+0.2,by = 0.05))) %>% unnest() %>% 
      group_by(draw, ID) %>% 
      mutate(p = psycho(x, (brms::inv_logit_scaled(alpha1)-0.5)*2,exp(beta1)),
             resp = rbinom(n(),1,p)) %>% 
      mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
      mutate(conf_mu = ifelse(ACC == 1, psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1) + meta_un_cor1), psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1) + meta_un_inc1)))
    
    
    df1_sub = bind_rows(
      pred_subj %>%
        mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
        group_by(x,ID) %>%
        summarize(
          name = "Type-1",
          mean = mean(p, na.rm = T),
          q5 = quantile(p, 0.05),
          q95 = quantile(p, 0.95),
          .groups = "drop"
        ),
      pred_subj %>%
        mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
        group_by(x,ID, Correct) %>%
        summarize(name = "Confidence",
                  mean = mean(conf_mu, na.rm = T),
                  q5 = quantile(conf_mu, 0.05),
                  q95 = quantile(conf_mu, 0.95),
                  .groups = "drop")
    )  %>% mutate(variable = name, name = NULL)
    
    
    
    behplot = rbind(
      bin = df %>%
        mutate(Correct = NA) %>%
        group_by(X,ID,Correct) %>%
        summarize(
          name = "Type-1",
          mean = mean(Y, na.rm = T),
          q5 = mean(Y, na.rm = T) - 2* (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
          q95 = mean(Y, na.rm = T) + 2 * (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
          .groups = "drop"
        ),
      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X,ID, Correct) %>%
        summarize(name = "Confidence",
                  mean = mean(Confidence, na.rm = T),
                  q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                  q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
                  .groups = "drop")
    ) %>% mutate(variable = name, name = NULL)
    
    sub_plots = list()
    for(name in unique(df1_sub$ID)){
      
      sub_plot = behplot  %>% 
        filter(variable != "RT" & ID == name) %>% 
        ggplot() +
        geom_pointrange(data = behplot%>% filter(variable != "RT"& ID == name), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                        shape = 21, color = "black", alpha = 0.5) +
        facet_wrap(~variable, scales = "free_y", ncol = 1) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
        theme_classic(base_size = 14) +
        labs(color = "Correct", fill = "Correct",
             y = "Value") +
        geom_vline(xintercept = 0, linetype = 2) +
        labs(subtitle = name)+
        scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
        # theme(legend.position = "top")+
        geom_line(data = df1_sub %>% filter(ID == name), aes(x = x, y = mean, col = Correct))+
        geom_ribbon(data = df1_sub%>% filter(ID == name), aes(x = x, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.5)
      
      sub_plots[[name]] = sub_plot
      
    }

    
    
    
    
    return(list(group_plot,sub_plots))
    
  }
  
}