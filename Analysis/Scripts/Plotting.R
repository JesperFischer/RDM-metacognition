

plot_beh_data = function(df,n_bins,ACC = F){
  
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

  if(ACC){
    bin = df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        mean = mean(ACC),
        q5 = mean(ACC) - 2* (mean(ACC) * (1-mean(ACC)) / sqrt(n())),
        q95 = mean(ACC) + 2 * (mean(ACC) * (1-mean(ACC)) / sqrt(n())),
        .groups = "drop"
      )
  }else{
    bin = df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        mean = mean(Y),
        q5 = mean(Y) - 2* (mean(Y) * (1-mean(Y)) / sqrt(n())),
        q95 = mean(Y) + 2 * (mean(Y) * (1-mean(Y)) / sqrt(n())),
        .groups = "drop"
      )
  }
  
  
  
    
  df1 = bind_rows(
    bin,
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT),
                q5 = mean(RT) - 2 * (sd(RT) / sqrt(n())),
                q95 = mean(RT) + 2 * (sd(RT) / sqrt(n())),
                .groups = "drop"),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence),
                q5 = mean(Confidence) - 2 * (sd(Confidence) / sqrt(n())),
                q95 = mean(Confidence) + 2 * (sd(Confidence) / sqrt(n())),
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
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())

  
  
  
  
  return(plot_mean)
  
  
  }

Plot_group_predictive_psycho = function(predictions, df, n_bins = NULL) {
  
  cutoff = 2
  
  # Prepare observed data
  if (!is.null(n_bins)) {
    # Create common bin boundaries based on the range of both datasets
    all_X <- c(df$X, predictions$X)
    X_range <- range(all_X, na.rm = TRUE)
    bin_breaks <- seq(X_range[1], X_range[2], length.out = n_bins + 1)
    
    # Calculate bin centers (midpoints)
    bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2
    
    # Bin the data using the common breaks and assign bin centers
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
             X = bin_centers[X_bin]) %>%
      select(-X_bin)
    
    # Bin the predictions using the same common breaks and assign same bin centers
    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
             X = bin_centers[X_bin]) %>%
      select(-X_bin)
  }
  
  dataq = bind_rows(
    # df %>%
    #   mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
    #   group_by(X) %>%
    #   summarize(
    #     name = "Type-1",
    #     n = n(),
    #     k = sum(Y),
    #     mean = (1 + k) / (2 + n),
    #     q5  = qbeta(0.05, 1 + k, 1 + n - k),
    #     q10 = qbeta(0.10, 1 + k, 1 + n - k),
    #     q20 = qbeta(0.20, 1 + k, 1 + n - k),
    #     q80 = qbeta(0.80, 1 + k, 1 + n - k),
    #     q90 = qbeta(0.90, 1 + k, 1 + n - k),
    #     q95 = qbeta(0.95, 1 + k, 1 + n - k),
    #     .groups = "drop"
    #   ),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        mean = mean(ACC),
        q5 = mean(ACC) - 2* (mean(ACC) * (1-mean(ACC)) / sqrt(n())),
        q95 = mean(ACC) + 2 * (mean(ACC) * (1-mean(ACC)) / sqrt(n())),
        # n = n(),
        # k = sum(Y),
        # mean = (1 + k) / (2 + n),
        # q5  = qbeta(0.05, 1 + k, 1 + n - k),
        # q10 = qbeta(0.10, 1 + k, 1 + n - k),
        # q20 = qbeta(0.20, 1 + k, 1 + n - k),
        # q80 = qbeta(0.80, 1 + k, 1 + n - k),
        # q90 = qbeta(0.90, 1 + k, 1 + n - k),
        # q95 = qbeta(0.95, 1 + k, 1 + n - k),
        .groups = "drop"
      ),
    
    
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT),
                q5 = mean(RT) - 2 * (sd(RT) / sqrt(n())),
                q95 = mean(RT) + 2 * (sd(RT) / sqrt(n())),
                .groups = "drop"),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence),
                q5 = mean(Confidence) - 2 * (sd(Confidence) / sqrt(n())),
                q95 = mean(Confidence) + 2 * (sd(Confidence) / sqrt(n())),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff)
  
  # Prepare predicted data (using expected means)
  predictionsq_mean = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob ),
                q5 = quantile(prob , 0.05),
                q10 = quantile(prob , 0.1),
                q20 = quantile(prob , 0.2),
                q95 = quantile(prob , 0.95),
                q90 = quantile(prob , 0.90),
                q80 = quantile(prob , 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(rt_mu),
                q5 = quantile(rt_mu, 0.05),
                q10 = quantile(rt_mu, 0.1),
                q20 = quantile(rt_mu, 0.2),
                q95 = quantile(rt_mu, 0.95),
                q90 = quantile(rt_mu, 0.90),
                q80 = quantile(rt_mu, 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(conf_mu_actual),
                q5 = quantile(conf_mu_actual, 0.05),
                q10 = quantile(conf_mu_actual, 0.1),
                q20 = quantile(conf_mu_actual, 0.2),
                q95 = quantile(conf_mu_actual, 0.95),
                q90 = quantile(conf_mu_actual, 0.90),
                q80 = quantile(conf_mu_actual, 0.80),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))
  
  # Calculate trial-level residuals properly
  # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
  pred_mean_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X) %>%
    summarize(
      pred_prob_faster = mean(prob),
      pred_rt = mean(rt_mu),
      .groups = "drop"
    )
  
  # Join predictions to actual trial-level data
  df_with_pred = df %>%
    mutate(X = round(X, 2)) %>%
    left_join(pred_mean_per_trial, by = "X") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))
  
  # Confidence predictions need to be split by Correct
  pred_conf_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X, Correct) %>%
    summarize(pred_conf = mean(conf_mu_actual), .groups = "drop") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))
  
  df_with_pred = df_with_pred %>%
    left_join(pred_conf_per_trial %>% select(X, Correct_label, pred_conf),
              by = c("X", "Correct_label"))
  
  # Calculate residuals at trial level, then aggregate
  residuals_data = bind_rows(
    df_with_pred %>%
      mutate(residual = ACC - pred_prob_faster,
             name = "Type-1") %>%
      group_by(X) %>%
      summarize(residual_mean = mean(residual, na.rm = TRUE),
                residual_se = sd(residual, na.rm = TRUE) / sqrt(n()),
                name = first(name),
                .groups = "drop"),
    
    df_with_pred %>%
      mutate(residual = RT - pred_rt,
             name = "RT") %>%
      group_by(X) %>%
      summarize(residual_mean = mean(residual, na.rm = TRUE),
                residual_se = sd(residual, na.rm = TRUE) / sqrt(n()),
                name = first(name),
                .groups = "drop"),
    
    df_with_pred %>%
      mutate(residual = Confidence - pred_conf,
             name = "Confidence") %>%
      group_by(X, Correct_label) %>%
      summarize(residual_mean = mean(residual, na.rm = TRUE),
                residual_se = sd(residual, na.rm = TRUE) / sqrt(n()),
                name = first(name),
                .groups = "drop") %>%
      rename(Correct = Correct_label)
  )
  
  # Plot 1: Expected means (main plot)
  plot_mean = predictionsq_mean %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_pointrange(data = dataq, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                    shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct",
         # title = "Group predictions (expected means)",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  # Plot 2: Residuals
  plot_residuals = residuals_data %>%
    ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = residual_mean - 2*residual_se,
                        ymax = residual_mean + 2*residual_se),
                    alpha = 0.5, size = 0.3) +
    # geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 14) +
    labs(x = "Stimulus strength (X)", y = "(Obs - Pred)") +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none")
  
  library(patchwork)
  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))
  
  combined_plot
  
}

Plot_group_predictive_psycho = function(predictions, df, n_bins = NULL) {
  
  cutoff = 2
  
  dataq = bind_rows(
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        mean = (Y),
        .groups = "drop"
      ),
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = (RT),
                .groups = "drop"),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = (Confidence),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff)
  
  # Prepare predicted data (using expected means)
  predictionsq_mean = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob_faster ),
                q5 = quantile(prob_faster , 0.05),
                q10 = quantile(prob_faster , 0.1),
                q20 = quantile(prob_faster , 0.2),
                q95 = quantile(prob_faster , 0.95),
                q90 = quantile(prob_faster , 0.90),
                q80 = quantile(prob_faster , 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(rt_mu),
                q5 = quantile(rt_mu, 0.05),
                q10 = quantile(rt_mu, 0.1),
                q20 = quantile(rt_mu, 0.2),
                q95 = quantile(rt_mu, 0.95),
                q90 = quantile(rt_mu, 0.90),
                q80 = quantile(rt_mu, 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(conf_mu_actual),
                q5 = quantile(conf_mu_actual, 0.05),
                q10 = quantile(conf_mu_actual, 0.1),
                q20 = quantile(conf_mu_actual, 0.2),
                q95 = quantile(conf_mu_actual, 0.95),
                q90 = quantile(conf_mu_actual, 0.90),
                q80 = quantile(conf_mu_actual, 0.80),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))
  
  # Calculate trial-level residuals properly
  # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
  pred_mean_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X) %>%
    summarize(
      pred_prob_faster = mean(prob_faster),
      pred_rt = mean(rt_mu),
      .groups = "drop"
    )
  
  # Join predictions to actual trial-level data
  df_with_pred = df %>% select(Y,X,Confidence,RT,Correct) %>%
    mutate(X = round(X, 2)) %>%
    filter(abs(X) < cutoff) %>%
    left_join(pred_mean_per_trial, by = "X") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))
  
  # Confidence predictions need to be split by Correct
  pred_conf_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X, Correct) %>%
    summarize(pred_conf = mean(conf_mu_actual), .groups = "drop") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))
  
  df_with_pred = df_with_pred %>%
    left_join(pred_conf_per_trial %>% select(X, Correct_label, pred_conf),
              by = c("X", "Correct_label"))
  
  # Calculate residuals at trial level, then aggregate
  residuals_data = bind_rows(
    df_with_pred %>%
      mutate(residual = Y - pred_prob_faster,
             name = "Type-1") %>%
      group_by(X) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop"),
    
    df_with_pred %>%
      mutate(residual = RT - pred_rt,
             name = "RT") %>%
      group_by(X) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop"),
    
    df_with_pred %>%
      mutate(residual = Confidence - pred_conf,
             name = "Confidence") %>%
      group_by(X, Correct_label) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop") %>%
      rename(Correct = Correct_label)
  )
  
  
  
  
  # Plot 1: Expected means (main plot)
  plot_mean = predictionsq_mean %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    # mutate(name = ifelse(name == "RT","Response time",ifelse(name == "Type-1","Binary choice","Confidence"))) %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_point(data = dataq %>%
                 mutate(name = ifelse(name == "RT","Response time",
                                      ifelse(name == "Type-1","Binary choice","Confidence")),
                        name = factor(name, levels = c("Binary choice",
                                                       "Response time",
                                                       "Confidence")))
               , aes(x = X, y = mean, fill = Correct),
               shape = 21, color = "black", alpha = 0.5, size = 3) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    scale_color_manual(values = c("darkgreen","darkred","grey"))+
    scale_fill_manual(values = c("darkgreen","darkred","grey"))+
    theme_classic(base_size = 20) +
    labs(color = "Correct", fill = "Correct",
         # title = "Group predictions (expected means)",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top",
          legend.text = element_text(size = 20),      # text of legend items
          legend.title = element_text(size = 20),      # title of legend
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  
  # Plot 2: Residuals
  plot_residuals = residuals_data  %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_point(alpha = 0.5, size = 3) +
    # geom_smooth(method = "loess", se = F, alpha = 0.2) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 20) +
    scale_color_manual(values = c("darkgreen","darkred","grey"))+
    scale_fill_manual(values = c("darkgreen","darkred","grey"))+
    labs(x = "Coherence", y = "(Obs - Pred)") +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none")
  
  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))
  
  combined_plot
  
  
  ggsave(here::here("All_Siebe_data.tiff"),combined_plot ,dpi = 300,width = 24, height = 16, units = "cm")
  
  cop_cor = as_draws_df(fit$draws(c("rho_rt_conf","rho_p_rt","rho_p_conf"))) %>% 
    select(-contains(".")) %>% 
    pivot_longer(everything()) %>% 
    ggplot(aes(x = value))+geom_histogram(col = "black")+facet_wrap(~name)
  
  ggsave(here::here("Copulas_full.tiff"),cop_cor ,dpi = 300,width = 16, height = 10, units = "cm")
  
  
  #################### Probability of responding crorectly:
  
  cutoff = 2
  
  dataq = bind_rows(
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        mean = (ACC),
        .groups = "drop"
      ),
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = (RT),
                .groups = "drop"),
    
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = (Confidence),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff)
  
  # Prepare predicted data (using expected means)
  predictionsq_mean = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob    ),
                q5 = quantile(prob    , 0.05),
                q10 = quantile(prob , 0.1),
                q20 = quantile(prob , 0.2),
                q95 = quantile(prob , 0.95),
                q90 = quantile(prob , 0.90),
                q80 = quantile(prob , 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(rt_mu),
                q5 = quantile(rt_mu, 0.05),
                q10 = quantile(rt_mu, 0.1),
                q20 = quantile(rt_mu, 0.2),
                q95 = quantile(rt_mu, 0.95),
                q90 = quantile(rt_mu, 0.90),
                q80 = quantile(rt_mu, 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(conf_mu_actual),
                q5 = quantile(conf_mu_actual, 0.05),
                q10 = quantile(conf_mu_actual, 0.1),
                q20 = quantile(conf_mu_actual, 0.2),
                q95 = quantile(conf_mu_actual, 0.95),
                q90 = quantile(conf_mu_actual, 0.90),
                q80 = quantile(conf_mu_actual, 0.80),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))
  
  # Calculate trial-level residuals properly
  # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
  pred_mean_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X) %>%
    summarize(
      pred_prob_faster = mean(prob),
      pred_rt = mean(rt_mu),
      .groups = "drop"
    )
  
  # Join predictions to actual trial-level data
  df_with_pred = df %>% select(Y,X,Confidence,RT,Correct) %>%
    mutate(X = round(X, 2)) %>%
    filter(abs(X) < cutoff) %>%
    left_join(pred_mean_per_trial, by = "X") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))
  
  # Confidence predictions need to be split by Correct
  pred_conf_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X, Correct) %>%
    summarize(pred_conf = mean(conf_mu_actual), .groups = "drop") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))
  
  df_with_pred = df_with_pred %>%
    left_join(pred_conf_per_trial %>% select(X, Correct_label, pred_conf),
              by = c("X", "Correct_label"))
  
  # Calculate residuals at trial level, then aggregate
  residuals_data = bind_rows(
    df_with_pred %>%
      mutate(residual = Correct - pred_prob_faster,
             name = "Type-1") %>%
      group_by(X) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop"),
    
    df_with_pred %>%
      mutate(residual = RT - pred_rt,
             name = "RT") %>%
      group_by(X) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop"),
    
    df_with_pred %>%
      mutate(residual = Confidence - pred_conf,
             name = "Confidence") %>%
      group_by(X, Correct_label) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop") %>%
      rename(Correct = Correct_label)
  )
  
  
  
  
  # Plot 1: Expected means (main plot)
  plot_mean = predictionsq_mean %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    # mutate(name = ifelse(name == "RT","Response time",ifelse(name == "Type-1","Binary choice","Confidence"))) %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_point(data = dataq %>%
                 mutate(name = ifelse(name == "RT","Response time",
                                      ifelse(name == "Type-1","Binary choice","Confidence")),
                        name = factor(name, levels = c("Binary choice",
                                                       "Response time",
                                                       "Confidence"))) %>%
                 mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct)))
               , aes(x = X, y = mean, fill = Correct),
               shape = 21, color = "black", alpha = 0.5, size = 3) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    scale_color_manual(values = c("darkgreen","darkred","grey"))+
    scale_fill_manual(values = c("darkgreen","darkred","grey"))+
    theme_classic(base_size = 20) +
    labs(color = "Correct", fill = "Correct",
         # title = "Group predictions (expected means)",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top",
          legend.text = element_text(size = 20),      # text of legend items
          legend.title = element_text(size = 20),      # title of legend
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  
  # Plot 2: Residuals
  plot_residuals = residuals_data  %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    mutate(Correct = ifelse(name == "Binary choice" & residual_mean  < 0, "Incorrect",ifelse(name == "Binary choice" & residual_mean  > 0, "Correct",Correct))) %>%
    ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_point(alpha = 0.5, size = 3) +
    # geom_smooth(method = "loess", se = F, alpha = 0.2) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 20) +
    scale_color_manual(values = c("darkgreen","darkred","grey"))+
    scale_fill_manual(values = c("darkgreen","darkred","grey"))+
    labs(x = "Coherence", y = "(Obs - Pred)") +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none")
  
  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))
  
  combined_plot
  # ggsave(here::here("Siebe","Siebe_results_bin.tiff"),combined_plot ,dpi = 300,width = 30, height = 22, units = "cm")
  
  
  #######################################################
  
  
  
  # Prepare predicted data (using actual samples)
  predictionsq_preds = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob_faster ),
                q5 = quantile(prob_faster , 0.05),
                q10 = quantile(prob_faster , 0.1),
                q20 = quantile(prob_faster , 0.2),
                q95 = quantile(prob_faster , 0.95),
                q90 = quantile(prob_faster , 0.90),
                q80 = quantile(prob_faster , 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT_pred),
                q5 = quantile(RT_pred, 0.05),
                q10 = quantile(RT_pred, 0.1),
                q20 = quantile(RT_pred, 0.2),
                q95 = quantile(RT_pred, 0.95),
                q90 = quantile(RT_pred, 0.90),
                q80 = quantile(RT_pred, 0.80),
                .groups = "drop"),
    
    predictions %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence),
                q5 = quantile(Confidence, 0.05),
                q10 = quantile(Confidence, 0.1),
                q20 = quantile(Confidence, 0.2),
                q95 = quantile(Confidence, 0.95),
                q90 = quantile(Confidence, 0.90),
                q80 = quantile(Confidence, 0.80),
                .groups = "drop")
  ) %>%
    filter(abs(X) < cutoff) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))
  
  
  # Plot 2: Actual predictions
  plot_preds = predictionsq_preds %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_point(data = dataq, aes(x = X, y = mean, fill = as.factor(Correct)),
               shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free", ncol = 3) +
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct") +
    # ggtitle("Group predictions (posterior predictive samples)")+
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")
  
  return(list(
    plot_combined = combined_plot,
    plot_mean = plot_mean,
    plot_residuals = plot_residuals,
    plot_preds = plot_preds
  ))
}
