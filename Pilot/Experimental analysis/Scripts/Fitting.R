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
                  interval = df$interval)
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 100,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.95,
    max_treedepth = 12,
    # init  = 0,
    parallel_chains = 4)
  
  if(grepl("MCMC",model$stan_file())){
    fit$save_object(here::here("Pilot","Experimental analysis","MCMC","saved models",paste0("fit_",unique(df$ID),".rds")))
  }
  if(grepl("CANSANDRE",model$stan_file())){
    fit$save_object(here::here("Pilot","Experimental analysis","CANSANDRE","saved models",paste0("fit_",unique(df$ID),".rds")))
  }
  if(grepl("Heurestic",model$stan_file())){
    fit$save_object(here::here("Pilot","Experimental analysis","Heurestic models","saved models","single subject","x1",paste0("fit_",unique(df$ID),".rds")))
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
  parameters = c("sigma_choice_log","mean_choice","prec_conf_log","sigma_m_log","sigma_e_log", "meta_bias")
  
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
  }
  parameters = c("c0","c11","alpha1","beta1","alpha","beta","conf_prec1","meta_un_cor1","meta_un_inc1","meta_bias")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    
  }
  
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
  parameters = c("sigma_choice","mean_choice","prec_conf","sigma_m","sigma_e", "meta_bias")
  
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    
    
    behplot = plot_beh_data(df,NULL,F, r_data = T)
    
    preds = fit$draws(c("p_action","mu_conf","evidence")) %>% as_draws_df() %>% select(-contains(".")) %>% pivot_longer(everything()) %>% 
      separate(
        name,
        into = c("param_base", "trial"),
        sep = "\\[|\\]",
        extra = "drop"
      ) %>% 
      mutate(trial = as.integer(trial)) %>% 
      group_by(param_base,trial) %>% 
      summarize(mean = mean(value, na.rm =T),
                q5 = quantile(value, 0.05, na.rm =T),
                q95 = quantile(value, 0.95, na.rm =T)) %>% 
      pivot_wider(names_from = "param_base", values_from = c(mean,q5,q95))
    
    
    dq = inner_join(behplot,preds)
    
    
    Confidence = dq %>% ggplot(aes(x = Confidence, y = mean_mu_conf, ymin = q5_mu_conf, ymax = q95_mu_conf, col = as.factor(Correct)))+
      geom_pointrange()+
      geom_abline()+
      theme_classic()+
      labs(y= "Predicted confidence",
           x = "Reported confidence")+
      theme(legend.position = "top")+
      ggtitle(unique(behplot$ID))
    
    
    Actions = dq %>% 
      ggplot(aes(x = (Y), y = mean_p_action, ymin = q5_p_action, ymax = q95_p_action))+
      geom_pointrange(position = position_dodge2(width = 0.2), alpha = 0.3)+
      geom_abline()+
      theme_classic()+
      labs(y= "Predicted Probability of action",
           x = "Reported Action")+
      theme(legend.position = "top")+
      ggtitle(unique(behplot$ID))
    
    
    RT_evi = dq %>% 
      ggplot(aes(y = RT, x = mean_evidence, xmin = q5_evidence, xmax = q95_evidence))+
      geom_pointrange(position = position_dodge2(width = 0.2), alpha = 0.3)+
      theme_classic()+
      labs(y= "Response time (s)",
           x = "Estimated Evidence")+
      theme(legend.position = "top")+
      ggtitle(unique(behplot$ID))
    
    
    return(list(Confidence, Actions, RT_evi))
    
  }
  
  ## old
  parameters = c("c0","c11","alpha1","beta1","conf_prec1","meta_un_cor1","meta_un_inc1","meta_bias")
  if(length(intersect(parameters, available)) > 5){
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
      mutate(conf_mu = ifelse(ACC == 1,
                              psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1 - exp(meta_un_cor1))),
                              psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2,meta_un_inc1))) %>% 
      mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias)) %>% 
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
    
    
    
    pp_plot = behplot %>% 
      filter(name != "RT") %>%
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
    
    
    
    psycho = function(x,beta,sigma_e, sigma_k,lapse){
      return((brms::inv_logit_scaled(lapse) / 2) + (1-2*brms::inv_logit_scaled(lapse) / 2) * pnorm((x - (brms::inv_logit_scaled(beta)-0.5)*2) / (sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))))
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
      mutate(p = psycho(x, beta,sigma_e,sigma_k,lapse)) %>% 
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



dia_hier = function(fit,df){
  
  
  available <- names(as_draws_df(fit$draws()))
  parameters = c("sigma_choice[1]","mean_choice[1]","conf_prec[1]","sigma_m[1]","sigma_e[1]")
  
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
    
    
    n_subj = length(unique(df$ID))
    subs = rep(unique(df$ID), length(params))
    
    subj_parameters = str_sub(parameters, 1, -4)
    
    subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
    
  }
  
  ######## old
  parameters = c("c0[1]","c11[1]","alpha1[1]","beta1[1]","conf_prec1[1]","meta_un_cor1[1]","meta_un_inc1[1]","meta_bias[1]","lapse[1]","meta_un_beta[1]","meta_bias_beta[1]")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    n_subj = length(unique(df$ID))
    subs = rep(unique(df$ID), length(params))
    
    subj_parameters = str_sub(parameters, 1, -4)
    
    subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
  }
  
  parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","lapse[1]","sigmam_beta[1]","meta_bias_beta[1]")
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


pp_hier = function(fit,df,n_bins){
  
  
  available <- names(as_draws_df(fit$draws()))
  parameters = c("sigma_choice[1]","mean_choice[1]","prec_conf[1]","sigma_m[1]","sigma_e[1]")
  
  if(length(intersect(parameters, available)) != 0){
    params <- intersect(parameters, available)
  }
  
  ##############
  ## real OLD
  ##############
  
  # parameters = c("c0[1]","c11[1]","alpha1[1]","beta1[1]","conf_prec1[1]","meta_un_cor1[1]","meta_un_inc1[1]","meta_bias[1]","lapse[1]","meta_un_beta[1]","meta_bias_beta[1]")
  # if(length(intersect(parameters, available)) > 5){
  #   params <- intersect(parameters, available)
  #   
  #   n_subj = length(unique(df$ID))
  #   subs = rep(unique(df$ID), length(params))
  #   
  #   subj_parameters = str_sub(parameters, 1, -4)
  #   
  #   subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
  #   
  #   
  #   psycho = function(x,alpha,beta,lapse){
  #     return(pnorm(beta * (x-alpha)))
  #   }
  #   
  #   psycho_ACC = function(x,alpha,beta){
  #     return(pnorm(beta * abs(x-alpha)))
  #   }
  #   
  #   sum = fit$summary(c("gm","tau_u")) %>% dplyr::select(-c(mad,median))
  #   sub = fit$summary(c(subj_parameters)) %>% dplyr::select(-c(mad,median)) %>% mutate(ID = subs)
  #   
  #   pred_data_group = as_draws_df(fit$draws("gm")) %>% select(-contains(".")) %>% 
  #     rename_with(~c(  "alpha1",
  #                      "beta1",
  #                      "conf_prec1",
  #                      "meta_un_cor1",
  #                      "meta_un_inc1",
  #                      "meta_bias",
  #                      "lapse",
  #                      "meta_un_beta",
  #                      "meta_bias_beta")) %>%
  #     mutate(draw = 1:n()) %>% 
  #     mutate(x = list(seq(min(df$X)-0.2,max(df$X)+0.2,by = 0.05))) %>% 
  #     unnest() %>% 
  #     mutate(interval = list(seq(1,5,by = 0.1))) %>% 
  #     unnest() %>% 
  #     group_by(draw) %>% 
  #     mutate(p = psycho(x, (brms::inv_logit_scaled(alpha1)-0.5)*2,exp(beta1) , brms::inv_logit_scaled(lapse) / 2),
  #            resp = rbinom(n(),1,p)) %>% 
  #     mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
  #     mutate(conf_mu = ifelse(ACC == 1,
  #                             psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1 - exp(meta_un_cor1) + interval * meta_un_beta)),
  #                             psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, meta_un_inc1))) %>% 
  #     mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + interval * meta_bias_beta))
  #   
  #   
  #   
  #   
  #   df1 = bind_rows(
  #     pred_data_group %>%
  #       mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
  #       group_by(x) %>%
  #       summarize(
  #         name = "Type-1",
  #         mean = mean(p, na.rm = T),
  #         q5 = quantile(p, 0.05),
  #         q95 = quantile(p, 0.95),
  #         .groups = "drop"
  #       ),
  #     pred_data_group %>%
  #       mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
  #       group_by(x, Correct) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(conf_mu, na.rm = T),
  #                 q5 = quantile(conf_mu, 0.05),
  #                 q95 = quantile(conf_mu, 0.95),
  #                 .groups = "drop")
  #   ) 
  #   
  #   
  #   
  #   # Prepare observed data
  #   if (!is.null(n_bins)) {
  #     # Create common bin boundaries based on the range of both datasets
  #     all_X <- c(df$X)
  #     X_range <- range(all_X, na.rm = TRUE)
  #     bin_breaks <- seq(X_range[1], X_range[2], length.out = n_bins + 1)
  #     
  #     # Calculate bin centers (midpoints)
  #     bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2
  #     
  #     # Bin the data using the common breaks and assign bin centers
  #     df <- df %>%
  #       mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
  #              X = bin_centers[X_bin]) %>%
  #       select(-X_bin)
  #     
  #   }
  #   
  #   
  #   behplot = rbind(
  #     bin = df %>%
  #       mutate(Correct = NA) %>%
  #       group_by(X,Correct) %>%
  #       summarize(
  #         name = "Type-1",
  #         mean = mean(Y, na.rm = T),
  #         q5 = mean(Y, na.rm = T) - 2* (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
  #         q95 = mean(Y, na.rm = T) + 2 * (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
  #         .groups = "drop"
  #       ),
  #     df %>%
  #       mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
  #       group_by(X, Correct) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(Confidence, na.rm = T),
  #                 q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 .groups = "drop")
  #   ) 
  #   
  #   
  #   group_plot = behplot %>% filter(name != "RT") %>% 
  #     ggplot() +
  #     geom_pointrange(data = behplot%>% filter(name != "RT"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
  #                     shape = 21, color = "black", alpha = 0.5) +
  #     facet_wrap(~name, scales = "free_y", ncol = 1) +
  #     scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  #     theme_classic(base_size = 14) +
  #     labs(color = "Correct", fill = "Correct",
  #          y = "Value") +
  #     geom_vline(xintercept = 0, linetype = 2) +
  #     scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
  #     # theme(legend.position = "top")+
  #     geom_line(data = df1, aes(x = x, y = mean, col = Correct))+
  #     geom_ribbon(data = df1, aes(x = x, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.5)
  #   
  #   
  #   
  #   pred_subj = as_draws_df(fit$draws(c("alpha1","beta1","meta_un_cor1",
  #                                       "meta_un_inc1","meta_bias",
  #                                       "lapse",
  #                                       "meta_un_beta",
  #                                       "meta_bias_beta"
  #   ))) %>% select(-contains(".")) %>% 
  #     mutate(draw = 1:n()) %>% 
  #     pivot_longer(-draw) %>% 
  #     mutate(
  #       sub_idx = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
  #       param   = str_remove(name, "\\[\\d+\\]"),
  #       ID  = subs[sub_idx],
  #       name = NULL
  #     ) %>% 
  #     pivot_wider(names_from = "param", values_from = "value") %>% 
  #     mutate(x = list(seq(min(df$X)-0.2,max(df$X)+0.2,by = 0.05))) %>% 
  #     unnest() %>% 
  #     mutate(interval = list(seq(1,5,by = 0.1))) %>% 
  #     unnest() %>% 
  #     group_by(draw, ID) %>% 
  #     mutate(p = psycho(x, (brms::inv_logit_scaled(alpha1)-0.5)*2,exp(beta1), brms::inv_logit_scaled(lapse) / 2),
  #            resp = rbinom(n(),1,p)) %>% 
  #     mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
  #     mutate(conf_mu = ifelse(ACC == 1,
  #                             psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, exp(beta1 - exp(meta_un_cor1) + interval * meta_un_beta)),
  #                             psycho_ACC(x,(brms::inv_logit_scaled(alpha1)-0.5)*2, meta_un_inc1))) %>% 
  #     mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias+ interval * meta_bias_beta))
  #   
  #   
  #   
  #   df1_sub = bind_rows(
  #     pred_subj %>%
  #       mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
  #       group_by(x,ID) %>%
  #       summarize(
  #         name = "Type-1",
  #         mean = mean(p, na.rm = T),
  #         q5 = quantile(p, 0.05),
  #         q95 = quantile(p, 0.95),
  #         .groups = "drop"
  #       ),
  #     pred_subj %>%
  #       mutate(Correct = ifelse(ACC == 1, "Correct", "Incorrect")) %>%
  #       group_by(x,ID, Correct) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(conf_mu, na.rm = T),
  #                 q5 = quantile(conf_mu, 0.05),
  #                 q95 = quantile(conf_mu, 0.95),
  #                 .groups = "drop")
  #   )  %>% mutate(variable = name, name = NULL)
  #   
  #   
  #   
  #   behplot = rbind(
  #     bin = df %>%
  #       mutate(Correct = NA) %>%
  #       group_by(X,ID,Correct) %>%
  #       summarize(
  #         name = "Type-1",
  #         mean = mean(Y, na.rm = T),
  #         q5 = mean(Y, na.rm = T) - 2* (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
  #         q95 = mean(Y, na.rm = T) + 2 * (mean(Y, na.rm = T) * (1-mean(Y, na.rm = T)) / sqrt(n())),
  #         .groups = "drop"
  #       ),
  #     df %>%
  #       mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
  #       group_by(X,ID, Correct) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(Confidence, na.rm = T),
  #                 q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 .groups = "drop")
  #   ) %>% mutate(variable = name, name = NULL)
  #   
  #   sub_plots = list()
  #   for(name in unique(df1_sub$ID)){
  #     
  #     sub_plot = behplot  %>% 
  #       filter(variable != "RT" & ID == name) %>% 
  #       ggplot() +
  #       geom_pointrange(data = behplot%>% filter(variable != "RT"& ID == name), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
  #                       shape = 21, color = "black", alpha = 0.5) +
  #       facet_wrap(~variable, scales = "free_y", ncol = 1) +
  #       scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  #       theme_classic(base_size = 14) +
  #       labs(color = "Correct", fill = "Correct",
  #            y = "Value") +
  #       geom_vline(xintercept = 0, linetype = 2) +
  #       labs(subtitle = name)+
  #       scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
  #       # theme(legend.position = "top")+
  #       geom_line(data = df1_sub %>% filter(ID == name), aes(x = x, y = mean, col = Correct))+
  #       geom_ribbon(data = df1_sub%>% filter(ID == name), aes(x = x, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.5)
  #     
  #     sub_plots[[name]] = sub_plot
  #     
  #   }
  #   
  #   
  #   
  #   
  #   
  #   return(list(group_plot,sub_plots))
  #   
  # }
  # 
  # ## Old
  # parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","lapse[1]","sigmam_beta[1]","meta_bias_beta[1]")
  # if(length(intersect(parameters, available)) > 5){
  #   params <- intersect(parameters, available)
  #   
  #   n_subj = length(unique(df$ID))
  #   subs = rep(unique(df$ID), length(params))
  #   
  #   subj_parameters = str_sub(parameters, 1, -4)
  #   
  #   subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
  #   
  #   
  #   psycho = function(x,beta,sigma_e, sigma_k,lapse){
  #     return((brms::inv_logit_scaled(lapse) / 2) + (1-2*brms::inv_logit_scaled(lapse) / 2) * pnorm((x - (brms::inv_logit_scaled(beta)-0.5)*2) / (sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))))
  #   }
  #   
  #   confs = function(action,X,beta,sigma_e,sigma_k,sigma_m){
  #     
  #     mu_e = X - (brms::inv_logit_scaled(beta)-0.5)*2
  #     sigma_total = sqrt(exp(sigma_k)^2 + exp(sigma_e)^2)
  #     z = mu_e / sigma_total
  #     
  #     # Correction factor
  #     correction_factor = exp(sigma_e)^2 / sigma_total
  #     
  #     # Conditional expectation (different Mills ratios for each choice)
  #     e_cond = ifelse(action == 1,
  #                     mu_e + correction_factor * dnorm(z) / pnorm(z),
  #                     mu_e - correction_factor * dnorm(z) / pnorm(-z))  # Note: pnorm(-z) here!
  #     
  #     # Now use e_cond in the confidence formula
  #     Cc = 2/1.7
  #     sigma2 = exp(sigma_m)^2 + exp(sigma_e)^2
  #     
  #     c_mu = ifelse(action == 1,
  #                   pnorm((1/(sqrt(1 + (Cc^2*abs(X)^2) / sigma2))) * (e_cond * abs(X) * Cc) / sigma2),
  #                   1-pnorm((1/(sqrt(1 + (Cc^2*abs(X)^2) / sigma2))) * (e_cond * abs(X) * Cc) / sigma2))
  #     
  #     return(c_mu)
  #     
  #   }
  #   
  #   sum = fit$summary(c("gm","tau_u")) %>% dplyr::select(-c(mad,median))
  #   sub = fit$summary(c(subj_parameters)) %>% dplyr::select(-c(mad,median)) %>% mutate(ID = subs)
  #   
  #   pred_data_group = as_draws_df(fit$draws("gm")) %>% select(-contains(".")) %>% 
  #     rename_with(~c(  "beta",
  #                      "sigma_e",
  #                      "sigma_k",
  #                      "sigma_m",
  #                      "meta_bias",
  #                      "lapse",
  #                      "confprec",
  #                      "sigmam_beta",
  #                      "meta_bias_beta")) %>%
  #     mutate(draw = 1:n()) %>% 
  #     mutate(x = list(seq(-1,1,by = 0.1))) %>% 
  #     unnest() %>% 
  #     mutate(interval = mean(df$interval)) %>% 
  #     unnest() %>% 
  #     group_by(draw) %>% 
  #     mutate(p = psycho(x, beta,sigma_e,sigma_k,lapse)) %>% 
  #     # resp = rbinom(n(),1,p)) %>% 
  #     # mutate(ACC = ifelse(resp == 0 & x < 0,1, ifelse(resp == 1 & x > 0, 1, 0))) %>% 
  #     rowwise() %>% 
  #     mutate(action = rbinom(1,1,p)) %>% 
  #     unnest() %>% 
  #     rowwise() %>% 
  #     mutate(conf_mu = ifelse(action == 1,
  #                             confs(1,x,beta, sigma_e, sigma_k,sigma_m + interval * sigmam_beta ),
  #                             confs(0,x,beta, sigma_e, sigma_k,sigma_m+ interval * sigmam_beta ))) %>% 
  #     mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + interval * meta_bias_beta)) %>% 
  #     ungroup()
  #   
  #   
  #   
  #   
  #   df1 = bind_rows(
  #     pred_data_group %>%
  #       mutate(action = ifelse(action == 1, "up", "down")) %>%
  #       group_by(x) %>%
  #       summarize(
  #         name = "Type-1",
  #         mean = mean(p, na.rm = T),
  #         q5 = quantile(p, 0.05),
  #         q95 = quantile(p, 0.95),
  #         .groups = "drop"
  #       ),
  #     pred_data_group %>%
  #       mutate(action = ifelse(action == 1, "up", "down")) %>%
  #       group_by(x, action) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(conf_mu, na.rm = T),
  #                 q5 = quantile(conf_mu, 0.05),
  #                 q95 = quantile(conf_mu, 0.95),
  #                 .groups = "drop")
  #   ) 
  #   
  #   
  #   
  #   # Prepare observed data
  #   if (!is.null(n_bins)) {
  #     # Create common bin boundaries based on the range of both datasets
  #     df_g <- df %>%
  #       mutate(
  #         X_bin = cut(
  #           X,
  #           breaks = seq(min(X, na.rm = TRUE),
  #                        max(X, na.rm = TRUE),
  #                        length.out = n_bins + 1),
  #           labels = FALSE,
  #           include.lowest = TRUE
  #         ),
  #         # compute bin centers per subject
  #         X = {
  #           breaks_i <- seq(min(X, na.rm = TRUE),
  #                           max(X, na.rm = TRUE),
  #                           length.out = n_bins + 1)
  #           centers_i <- (breaks_i[-1] + breaks_i[-length(breaks_i)]) / 2
  #           centers_i[X_bin]
  #         }
  #       ) %>%
  #       ungroup() %>%
  #       select(-X_bin)
  #   }
  #   
  #   
  #   behplot = rbind(
  #     bin = df_g %>%
  #       # mutate(action = ifelse(Y == 1, "up", "down")) %>%
  #       mutate(action = NA) %>%
  #       group_by(X,action) %>%
  #       summarize(
  #         name = "Type-1",
  #         k = sum(Y),
  #         n = n(),
  #         mean = k / n,
  #         
  #         a_post = k + 1,
  #         b_post = (n - k) + 1,
  #         
  #         q5 = qbeta(0.025, a_post, b_post),
  #         q95 = qbeta(0.975, a_post, b_post),
  #         .groups = "drop") %>% mutate(k = NULL, n = NULL, a_post = NULL, b_post = NULL),
  #     df_g %>%
  #       mutate(action = ifelse(Y == 1, "up", "down")) %>%
  #       group_by(X, action) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(Confidence, na.rm = T),
  #                 q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 .groups = "drop")
  #   ) 
  #   
  #   
  #   group_plot = behplot %>% filter(name != "RT") %>% 
  #     ggplot() +
  #     geom_pointrange(data = behplot%>% filter(name != "RT"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
  #                     shape = 21, color = "black", alpha = 0.5) +
  #     facet_wrap(~name, scales = "free_y", ncol = 1) +
  #     scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  #     theme_classic(base_size = 14) +
  #     labs(color = "action", fill = "action",
  #          y = "Value") +
  #     geom_vline(xintercept = 0, linetype = 2) +
  #     scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
  #     geom_line(data = df1, aes(x = x, y = mean, col = action))+
  #     geom_ribbon(data = df1, aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)
  #   
  #   
  #   
  #   pred_subj = as_draws_df(fit$draws(c(  "beta",
  #                                         "sigma_e",
  #                                         "sigma_k",
  #                                         "sigma_m",
  #                                         "meta_bias",
  #                                         "lapse",
  #                                         "confprec",
  #                                         "sigmam_beta",
  #                                         "meta_bias_beta"
  #   ))) %>% select(-contains(".")) %>% 
  #     mutate(draw = 1:n()) %>% 
  #     pivot_longer(-draw) %>% 
  #     mutate(
  #       sub_idx = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
  #       param   = str_remove(name, "\\[\\d+\\]"),
  #       ID  = subs[sub_idx],
  #       name = NULL
  #     ) %>% 
  #     pivot_wider(names_from = "param", values_from = "value") %>% 
  #     mutate(x = list(seq(-1,1,by = 0.1))) %>% 
  #     unnest() %>% 
  #     mutate(interval = mean(df$interval)) %>% 
  #     unnest() %>% 
  #     group_by(draw, ID) %>% 
  #     mutate(p = psycho(x, beta,sigma_e,sigma_k,lapse)) %>% 
  #     rowwise() %>% 
  #     mutate(action = rbinom(1,1,p)) %>% 
  #     unnest() %>% 
  #     rowwise() %>% 
  #     mutate(conf_mu = ifelse(action == 1,
  #                             confs(1,x,beta, sigma_e, sigma_k,sigma_m + interval * sigmam_beta ),
  #                             confs(0,x,beta, sigma_e, sigma_k,sigma_m+ interval * sigmam_beta ))) %>% 
  #     mutate(conf_mu = brms::inv_logit_scaled(brms::logit_scaled(conf_mu) + meta_bias + interval * meta_bias_beta)) %>% 
  #     ungroup()
  #   
  #   
  #   df1_sub = bind_rows(
  #     pred_subj %>%
  #       mutate(action = ifelse(action == 1, "up", "down")) %>%
  #       group_by(x,ID) %>%
  #       summarize(
  #         name = "Type-1",
  #         mean = mean(p, na.rm = T),
  #         q5 = quantile(p, 0.05),
  #         q95 = quantile(p, 0.95),
  #         .groups = "drop"
  #       ),
  #     pred_subj %>%
  #       mutate(action = ifelse(action == 1, "up", "down")) %>%
  #       group_by(x,ID, action) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(conf_mu, na.rm = T),
  #                 q5 = quantile(conf_mu, 0.05),
  #                 q95 = quantile(conf_mu, 0.95),
  #                 .groups = "drop")
  #   )  %>% mutate(variable = name, name = NULL)
  #   
  #   
  #   
  #   
  #   # Prepare observed data
  #   if (!is.null(n_bins)) {
  #     # Create common bin boundaries based on the range of both datasets
  #     df_s <- df %>%
  #       group_by(ID) %>% 
  #       mutate(
  #         X_bin = cut(
  #           X,
  #           breaks = seq(min(X, na.rm = TRUE),
  #                        max(X, na.rm = TRUE),
  #                        length.out = n_bins + 1),
  #           labels = FALSE,
  #           include.lowest = TRUE
  #         ),
  #         # compute bin centers per subject
  #         X = {
  #           breaks_i <- seq(min(X, na.rm = TRUE),
  #                           max(X, na.rm = TRUE),
  #                           length.out = n_bins + 1)
  #           centers_i <- (breaks_i[-1] + breaks_i[-length(breaks_i)]) / 2
  #           centers_i[X_bin]
  #         }
  #       ) %>%
  #       ungroup() %>%
  #       select(-X_bin)
  #   }
  #   
  #   
  #   
  #   behplot = rbind(
  #     bin = df_s %>%
  #       mutate(action = NA) %>%
  #       group_by(X,ID,action) %>%
  #       summarize(
  #         name = "Type-1",
  #         k = sum(Y),
  #         n = n(),
  #         mean = k / n,
  #         
  #         a_post = k + 1,
  #         b_post = (n - k) + 1,
  #         
  #         q5 = qbeta(0.025, a_post, b_post),
  #         q95 = qbeta(0.975, a_post, b_post),
  #         .groups = "drop") %>% mutate(k = NULL, n = NULL, a_post = NULL, b_post = NULL),
  #     df_s %>%
  #       mutate(action = ifelse(Y == 1, "up", "down")) %>%
  #       group_by(X,ID, action) %>%
  #       summarize(name = "Confidence",
  #                 mean = mean(Confidence, na.rm = T),
  #                 q5 = mean(Confidence, na.rm = T) - 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 q95 = mean(Confidence, na.rm = T) + 2 * (sd(Confidence, na.rm = T) / sqrt(n())),
  #                 .groups = "drop")
  #   ) %>% mutate(variable = name, name = NULL)
  #   
  #   sub_plots = list()
  #   qq = 0
  #   for(name in unique(df1_sub$ID)){
  #     
  #     parameters = pred_subj %>% filter(ID == name) %>% 
  #       pivot_longer(cols = c("sigma_e","sigma_m","sigma_k")) %>%
  #       select(ID,name,value) %>% distinct() %>% 
  #       mutate(name = ifelse(name == "sigma_e","S_e",ifelse(name == "sigma_m","S_m",ifelse(name == "sigma_k","S_c",name)))) %>% 
  #       group_by(name) %>% summarize(mean = mean(value, na.rm = T),
  #                                    q5 = quantile(value,0.05, na.rm = T),
  #                                    q95 = quantile(value,0.95, na.rm = T)) %>% 
  #       mutate(label = paste0(name," = ",round(mean,2)," [",round(q5,2)," ; ",round(q95,2),"]")) %>% 
  #       mutate(variable = "Confidence") %>% mutate(param = c(0,0.15,0.3))
  #     
  #     qq = qq+1
  #     sub_plot = behplot  %>% 
  #       filter(variable != "RT" & ID == name) %>% 
  #       ggplot() +
  #       geom_pointrange(data = behplot%>% filter(variable != "RT"& ID == name), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
  #                       shape = 21, color = "black", alpha = 0.5) +
  #       facet_wrap(~variable, scales = "free_y", ncol = 1) +
  #       scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  #       theme_classic(base_size = 14) +
  #       labs(color = "action", fill = "action",
  #            y = "P(a=1) Or P(a=D)", x = "Signed Stimulus (XD)") +
  #       geom_vline(xintercept = 0, linetype = 2) +
  #       # labs(subtitle = name)+
  #       labs(subtitle = paste0("ID =", qq))+
  #       scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
  #       # geom_text(data = parameters, aes(x = 0, y = param,label = label), size = 4)+
  #       geom_line(data = df1_sub %>% filter(ID == name), aes(x = x, y = mean, col = action))+
  #       geom_ribbon(data = df1_sub%>% filter(ID == name), aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)
  #     
  #     
  #     sub_plots[[name]] = sub_plot
  #     
  #   }
  #   
  #   
  #   
  #   
  #   
  #   return(list(group_plot,sub_plots))
  #   
  # }
  # 
  
  
  #########
  ## New ##
  #########
  
  
  parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","lapse[1]","sigmam_beta[1]","meta_bias_beta[1]")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    n_subj = length(unique(df$ID))
    subs = rep(unique(df$ID), length(params))
    
    subj_parameters = str_sub(parameters, 1, -4)
    
    subj_parameters = paste0(rep(subj_parameters, each = n_subj),"[",rep(seq_len(n_subj), times = length(subj_parameters)),"]")
    
    
    psycho = function(x,beta,sigma_e, sigma_k,lapse){
      return((brms::inv_logit_scaled(lapse) / 2) + (1-2*brms::inv_logit_scaled(lapse) / 2) * pnorm((x - (brms::inv_logit_scaled(beta)-0.5)*2) / (sqrt(exp(sigma_k)^2 + exp(sigma_e)^2))))
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
                       "lapse",
                       "confprec",
                       "sigmam_beta",
                       "meta_bias_beta")) %>%
      mutate(draw = 1:n()) %>% 
      mutate(x = list(seq(-1,1,by = 0.1))) %>% 
      unnest() %>% 
      mutate(interval = mean(df$interval)) %>% 
      unnest() %>% 
      group_by(draw) %>% 
      mutate(p = psycho(x, beta,sigma_e,sigma_k,lapse)) %>% 
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
    
    
    group_plot = behplot %>% filter(name != "RT") %>% 
      ggplot() +
      geom_pointrange(data = behplot%>% filter(name != "RT"), aes(x = X, y = mean, ymin = q5, ymax = q95, fill = action),
                      shape = 21, color = "black", alpha = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
      theme_classic(base_size = 14) +
      labs(color = "action", fill = "action",
           y = "Value") +
      geom_vline(xintercept = 0, linetype = 2) +
      scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
      geom_line(data = df1, aes(x = x, y = mean, col = action))+
      geom_ribbon(data = df1, aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)
    
    
    
    pred_subj = as_draws_df(fit$draws(c(  "beta",
                                          "sigma_e",
                                          "sigma_k",
                                          "sigma_m",
                                          "meta_bias",
                                          "lapse",
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
      ) %>% 
      pivot_wider(names_from = "param", values_from = "value") %>% 
      mutate(x = list(seq(-1,1,by = 0.1))) %>% 
      unnest() %>% 
      mutate(interval = mean(df$interval)) %>% 
      unnest() %>% 
      group_by(draw, ID) %>% 
      mutate(p = psycho(x, beta,sigma_e,sigma_k,lapse)) %>% 
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
      
      parameters = pred_subj %>% filter(ID == name) %>% 
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
        filter(variable != "RT" & ID == name) %>% 
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
        scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks(n = 5))+
        # geom_text(data = parameters, aes(x = 0, y = param,label = label), size = 4)+
        geom_line(data = df1_sub %>% filter(ID == name), aes(x = x, y = mean, col = action))+
        geom_ribbon(data = df1_sub%>% filter(ID == name), aes(x = x, y = mean, ymin = q5, ymax = q95, fill = action), alpha = 0.5)
      
      
      sub_plots[[name]] = sub_plot
      
    }
    
    
    
    
    
    return(list(group_plot,sub_plots))
    
  }
  
  
}


correlation_matrix = function(fit){
  
  available <- names(as_draws_df(fit$draws()))
  
  
  ## Old

  parameters = c("c0[1]","c11[1]","alpha1[1]","beta1[1]","conf_prec1[1]",
                 "meta_un_cor1[1]","meta_un_inc1[1]","meta_bias[1]", 
                 "lapse[1]",
                 "meta_un_beta[1]",
                 "meta_bias_beta[1]")
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    param_names <- c(
      "alpha1",
      "beta1",
      "conf_prec1",
      "meta_un_cor1",
      "meta_un_inc1",
      "meta_bias",
      "lapse",
      "meta_un_beta",
      "meta_bias_beta"
    )
    
    corr_df <- fit$draws("correlation_matrix") %>% 
      as_draws_df() %>% 
      pivot_longer(
        cols = everything(),
        names_to = "param",
        values_to = "value"
      ) %>% 
      separate(
        param,
        into = c("matrix", "row", "col"),
        sep = "\\[|,|\\]",
        extra = "drop"
      ) %>% 
      mutate(
        row = as.integer(row),
        col = as.integer(col),
        row_name = param_names[row],
        col_name = param_names[col]
      )
    
    
    corr_hdi_df <- corr_df %>% 
      mutate(lower_triangle = row > col) %>% 
      # filter(lower_triangle) %>% 
      group_by(row_name, col_name,lower_triangle) %>% 
      summarise(
        mean = mean(value),
        q5 = hdi(value, prob = 0.95)[,1],
        q95 = hdi(value, prob = 0.95)[,2],
        .groups = "drop"
      ) %>% drop_na() %>%
      mutate(
        label = sprintf(
          "%.2f \n [%.2f, %.2f]",
          mean, q5, q95
        )
      )
    
    plot = corr_hdi_df %>% 
      mutate(
        row_name = factor(row_name, levels = param_names),
        col_name = factor(col_name, levels = param_names),
        fill_plot  = ifelse(lower_triangle, mean, NA_real_),
        label_plot = ifelse(lower_triangle, label, "")
      ) %>% 
      ggplot(aes(x = col_name, y = row_name, fill = fill_plot)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        limits = c(-1, 1)
      ) +
      geom_text(
        aes(label = label_plot),
        size = 2.5,
        colour = "black"
      ) +
      coord_fixed() +
      theme_minimal() +
      labs(
        x = NULL,
        y = NULL,
        fill = "Mean corr"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    
    table_df <- corr_hdi_df %>% 
      select(row_name, col_name, mean, q5, q95) %>% 
      arrange(row_name, col_name) %>%
      mutate(sig = (q5 > 0) | (q95 < 0))
    
    
    
    table = flextable::flextable(corr_hdi_df %>% 
                                   filter(mean != 1) %>% 
                                   select(row_name, col_name, mean, q5, q95) %>% 
                                   arrange(row_name, col_name) %>%
                                   mutate(sig = (q5 > 0) | (q95 < 0))) %>% 
      bg(i  =~  (q5 > 0) | (q95 < 0), bg = "darkgreen")
    
    
    
    # parameter by parameter scatter plots:
    
    # Extract subject-level draws
    draws_df <- fit$draws(param_names) %>%  # replace "theta" with your per-subject parameters
      as_draws_df() %>% 
      select(-contains(".")) %>% 
      pivot_longer(
        cols = everything(),
        names_to = "param",
        values_to = "value"
      ) %>% 
      separate(
        param,
        into = c("param_base", "subject"),
        sep = "\\[|\\]",
        extra = "drop"
      ) %>% 
      mutate(
        subject = as.integer(subject),
      ) %>% 
      pivot_wider(names_from = param_base, values_from = value) %>% drop_na()
    
    # Create all combinations of parameters
    param_combos <- expand.grid(
      row_name = param_names,
      col_name = param_names,
      stringsAsFactors = FALSE
    ) %>% 
      mutate(lower_triangle = match(row_name, param_names) > match(col_name, param_names))
    
    # Join data for plotting
    scatter_df <- param_combos %>% 
      filter(lower_triangle) %>%  # only lower triangle
      rowwise() %>% 
      mutate(
        data_list = list(draws_df %>% select(subject, row_val = all_of(row_name), col_val = all_of(col_name)))
      ) %>% 
      unnest(cols = c(data_list)) %>% 
      group_by(row_name,col_name,subject) %>% 
      summarize(mean_row = mean(unlist(row_val)),
                q5_row = hdi(unlist(row_val))[,1],
                q95_row = hdi(unlist(row_val))[,2],
                mean_col = mean(unlist(col_val)),
                q5_col = hdi(unlist(col_val))[,1],
                q95_col = hdi(unlist(col_val))[,2]
      )
    
    scatter_plot <- scatter_df %>% 
      mutate(
        row_name = factor(row_name, levels = rev(param_names)),
        col_name = factor(col_name, levels = param_names),
      ) %>% 
      ggplot(aes(x = mean_col, y = mean_row)) +
      geom_pointrange(aes(ymin = q5_row,ymax = q95_row),alpha = 0.5, size = 0.05) +
      geom_pointrange(aes(xmin = q5_col,xmax = q95_col),alpha = 0.5, size = 0.05) +
      facet_grid(row_name ~ col_name, scales = "free") +
      theme_minimal() +
      labs(x = NULL, y = NULL)
    # geom_smooth(method = "lm", se = F)
    
    
    
    
    return(list(plot,scatter_plot,table))
    
  }
  
  
  parameters = c("c0[1]","c11[1]","beta[1]","sigma_e[1]","sigma_k[1]","sigma_m[1]","confprec[1]","meta_bias[1]","lapse[1]","sigmam_beta[1]","meta_bias_beta[1]")

  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    param_names <- c(
      "beta",
      "sigma_e",
      "sigma_k",
      "sigma_m",
      "confprec",
      "meta_bias",
      "lapse",
      "sigmam_beta",
      "meta_bias_beta"
    )
    
    corr_df <- fit$draws("correlation_matrix") %>% 
      as_draws_df() %>% 
      pivot_longer(
        cols = everything(),
        names_to = "param",
        values_to = "value"
      ) %>% 
      separate(
        param,
        into = c("matrix", "row", "col"),
        sep = "\\[|,|\\]",
        extra = "drop"
      ) %>% 
      mutate(
        row = as.integer(row),
        col = as.integer(col),
        row_name = param_names[row],
        col_name = param_names[col]
      )
    
    
    corr_hdi_df <- corr_df %>% 
      mutate(lower_triangle = row > col) %>% 
      # filter(lower_triangle) %>% 
      group_by(row_name, col_name,lower_triangle) %>% 
      summarise(
        mean = mean(value),
        q5 = hdi(value, prob = 0.95)[,1],
        q95 = hdi(value, prob = 0.95)[,2],
        .groups = "drop"
      ) %>% drop_na() %>%
      mutate(
        label = sprintf(
          "%.2f \n [%.2f, %.2f]",
          mean, q5, q95
        )
      )
    
    plot = corr_hdi_df %>% 
      mutate(
        row_name = factor(row_name, levels = param_names),
        col_name = factor(col_name, levels = param_names),
        fill_plot  = ifelse(lower_triangle, mean, NA_real_),
        label_plot = ifelse(lower_triangle, label, "")
      ) %>% 
      ggplot(aes(x = col_name, y = row_name, fill = fill_plot)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        limits = c(-1, 1)
      ) +
      geom_text(
        aes(label = label_plot),
        size = 2.5,
        colour = "black"
      ) +
      coord_fixed() +
      theme_minimal() +
      labs(
        x = NULL,
        y = NULL,
        fill = "Mean corr"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    
    table_df <- corr_hdi_df %>% 
      select(row_name, col_name, mean, q5, q95) %>% 
      arrange(row_name, col_name) %>%
      mutate(sig = (q5 > 0) | (q95 < 0))
    
    
    
    table = flextable::flextable(corr_hdi_df %>% 
                                   filter(mean != 1) %>% 
                                   select(row_name, col_name, mean, q5, q95) %>% 
                                   arrange(row_name, col_name) %>%
                                   mutate(sig = (q5 > 0) | (q95 < 0))) %>% 
      bg(i  =~  (q5 > 0) | (q95 < 0), bg = "darkgreen")
    
    
    
    # parameter by parameter scatter plots:
    
    # Extract subject-level draws
    draws_df <- fit$draws(param_names) %>%  # replace "theta" with your per-subject parameters
      as_draws_df() %>% 
      select(-contains(".")) %>% 
      pivot_longer(
        cols = everything(),
        names_to = "param",
        values_to = "value"
      ) %>% 
      separate(
        param,
        into = c("param_base", "subject"),
        sep = "\\[|\\]",
        extra = "drop"
      ) %>% 
      mutate(
        subject = as.integer(subject),
      ) %>% 
      pivot_wider(names_from = param_base, values_from = value) %>% drop_na()
    
    # Create all combinations of parameters
    param_combos <- expand.grid(
      row_name = param_names,
      col_name = param_names,
      stringsAsFactors = FALSE
    ) %>% 
      mutate(lower_triangle = match(row_name, param_names) > match(col_name, param_names))
    
    # Join data for plotting
    scatter_df <- param_combos %>% 
      filter(lower_triangle) %>%  # only lower triangle
      rowwise() %>% 
      mutate(
        data_list = list(draws_df %>% select(subject, row_val = all_of(row_name), col_val = all_of(col_name)))
      ) %>% 
      unnest(cols = c(data_list)) %>% 
      group_by(row_name,col_name,subject) %>% 
      summarize(mean_row = mean(unlist(row_val)),
                q5_row = hdi(unlist(row_val))[,1],
                q95_row = hdi(unlist(row_val))[,2],
                mean_col = mean(unlist(col_val)),
                q5_col = hdi(unlist(col_val))[,1],
                q95_col = hdi(unlist(col_val))[,2]
      )
    
    scatter_plot <- scatter_df %>% 
      mutate(
        row_name = factor(row_name, levels = rev(param_names)),
        col_name = factor(col_name, levels = param_names),
      ) %>% 
      ggplot(aes(x = mean_col, y = mean_row)) +
      geom_pointrange(aes(ymin = q5_row,ymax = q95_row),alpha = 0.5, size = 0.05) +
      geom_pointrange(aes(xmin = q5_col,xmax = q95_col),alpha = 0.5, size = 0.05) +
      facet_grid(row_name ~ col_name, scales = "free") +
      theme_minimal() +
      labs(x = NULL, y = NULL)
    # geom_smooth(method = "lm", se = F)
    
    
    
    return(list(plot,scatter_plot,table))
    
  }
  
  
  
  
  
}

get_marginal_estimates = function(fit){
  
  
  available <- names(as_draws_df(fit$draws()))
  parameters = c("c0[1]","c11[1]","alpha1[1]","beta1[1]","conf_prec1[1]",
                 "meta_un_cor1[1]","meta_un_inc1[1]","meta_bias[1]", 
                 "lapse[1]","meta_un_beta[1]","meta_bias_beta[1]")
  
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    param_names <- c(
      "alpha1",
      "beta1",
      "conf_prec1",
      "meta_un_cor1",
      "meta_un_inc1",
      "meta_bias",
      "lapse",
      "meta_un_beta",
      "meta_bias_beta")
    
    
    grouplevel_posterior <- fit$draws(c("gm","tau_u")) %>% 
      as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("mean_",param_names),paste0("tau_",param_names))) %>% 
      pivot_longer(
        cols = everything(),
        names_to = "variable",
        values_to = "value"
      ) %>% mutate(Posterior = "Posterior")
    
    
    
    n_draws = nrow(as_draws_df(fit$sampler_diagnostics()))
    priors <- tibble(
      variable = c(
        paste0("mean_", param_names),
        paste0("tau_",  param_names)
      ),
      mean = c(
        # gm (global means)
        0,  1,  3,  0,  0,  0, -4,0,0,
        # tau_u
        0,  0,  0,  0,  0,  0,0,0,0
      ),
      sd = c(
        # gm (global means)
        0.5, 2, 2, 2, 2, 0.5, 2, 0.5,0.5,
        # tau_u
        1,   2, 2, 2, 2, 1,2,0.5,0.5
      ))
    
    grouplevel_prior <- priors %>%
      mutate(value = map2(mean, sd, ~ rnorm(n_draws, .x, .y))) %>%
      unnest(value) %>% 
      mutate(Posterior = "Prior",
             mean = NULL,
             sd = NULL)
    
    priorposterior_mean = rbind(grouplevel_prior,grouplevel_posterior) %>% 
      filter(grepl("mean_",variable)) %>% 
      ggplot(aes(x = value,fill = Posterior))+
      geom_histogram(col = "black", position = "identity", alpha = 0.5)+
      facet_wrap(~variable, scales = "free")+
      theme_classic()+
      theme(legend.position = "top")
    
    priorposterior_sd = rbind(grouplevel_prior,grouplevel_posterior) %>% 
      filter(grepl("tau_",variable)) %>% 
      ggplot(aes(x = value,fill = Posterior))+
      geom_histogram(col = "black", position = "identity", alpha = 0.5)+
      facet_wrap(~variable, scales = "free")+
      theme_classic()+
      theme(legend.position = "top")
    
    
    return(list(priorposterior_mean,priorposterior_sd))
    
    
  }
  
  
  
  parameters = c("sigma_choice[1]","mean_choice[1]","conf_prec[1]","sigma_m[1]","sigma_e[1]","meta_bias[1]")
  
  if(length(intersect(parameters, available)) > 5){
    params <- intersect(parameters, available)
    
    param_names <- c(
      "mean_choice",
      "sigma_e",
      "conf_prec",
      "sigma_m",
      "sigma_choice",
      "meta_bias")
    
    
    grouplevel_posterior <- fit$draws(c("gm","tau_u")) %>% 
      as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("mean_",param_names),paste0("tau_",param_names))) %>% 
      pivot_longer(
        cols = everything(),
        names_to = "variable",
        values_to = "value"
      ) %>% mutate(Posterior = "Posterior")
    
    
    n_draws = nrow(as_draws_df(fit$sampler_diagnostics()))
    
    priors <- tibble(
      variable = c(
        paste0("mean_", param_names),
        paste0("tau_",  param_names)
      ),
      mean = c(
        # gm (global means)
        0, -2,  3,  -1,  -2,  0,
        # tau_u
        0,  0,  0,  0,  0,  0
      ),
      sd = c(
        # gm (global means)
        0.3, 1, 2, 1, 1, 0.5,
        # tau_u
        0.3, 2, 2, 2, 2, 0.5
      ))
    
    grouplevel_prior <- priors %>%
      mutate(value = map2(mean, sd, ~ rnorm(n_draws, .x, .y))) %>%
      unnest(value) %>% 
      mutate(Posterior = "Prior",
             mean = NULL,
             sd = NULL)
    
    priorposterior_mean = rbind(grouplevel_prior,grouplevel_posterior) %>% 
      filter(!grepl("tau_",variable)) %>% 
      ggplot(aes(x = value,fill = Posterior))+
      geom_histogram(col = "black", position = "identity", alpha = 0.5)+
      facet_wrap(~variable, scales = "free")+
      theme_classic()+
      theme(legend.position = "top")
    
    priorposterior_sd = rbind(grouplevel_prior,grouplevel_posterior) %>% 
      filter(grepl("tau_",variable)) %>% 
      ggplot(aes(x = value,fill = Posterior))+
      geom_histogram(col = "black", position = "identity", alpha = 0.5)+
      facet_wrap(~variable, scales = "free")+
      theme_classic()+
      theme(legend.position = "top")
    
    
    return(list(priorposterior_mean,priorposterior_sd))
    
  }
  
  
}


