library(future)
library(future.apply)
library(ordbetareg)

pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes,
               brms, patchwork, cowplot,ggpubr,flextable)



files = list.files((here::here("Simulations","Heurestic","Sequential Sampling","results","new")), full.names = T)
files = files[1:length(files)]
results_list = list()
q = 0
for(i in 1:length(files)){
  q = q+1
  results_list[[q]] = readRDS(files[i])
}

# results_list = results_list[-which(results_list == "Error")]

simulated_means = purrr::map_dfr(
  results_list[1:length(results_list)],
  ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 3))
) %>% filter(variable %in% c("gm[8]")) %>% 
  select(variable,simulated,draw_id) %>% 
  # select(ess_bulk, ess_tail, rhat)
  distinct()



renames = c("gm[1]" = "beta",
            "gm[2]" = "sigma_e",
            "gm[3]" = "sigma_k",
            "gm[4]" = "sigma_m",
            "gm[5]" = "meta_bias",
            "gm[6]" = "lapse",
            "gm[7]" = "conf_prec",
            "gm[8]" = "beta_sigma_m",
            "gm[9]" = "beta_meta_bias")


group_param_recov = purrr::map_dfr(
  results_list[1:length(results_list)],
  ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 3))
) %>%   
  filter(grepl("gm",variable)) %>%
  mutate(variable = recode(variable, !!!renames, .default = variable)) %>% 
  filter(div == 0) %>% 
  mutate(prior_sd = as.factor(prior_sd)) %>% 
  mutate(n_subjects = as.factor(n_subs)) %>% 
  ggplot(aes(x = simulated, y = mean,ymin = q5,ymax = q95))+
  geom_pointrange(aes(group = variable),alpha = 0.5)+
  facet_wrap(~variable, scales = "free")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 2))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
  theme_classic()+
  labs(x = "Simulated",y = "Parameter value, mean/95HDI")+
  ggtitle("Group level parameter recovery")+
  geom_abline(col = "red")



subject_param_recov = purrr::map_dfr(
  results_list[1:length(results_list)],
  ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 2))
) %>%   
  mutate(variable = ifelse(variable == "sigmam_beta","beta_sigma_m",ifelse(variable == "meta_bias_beta","beta_meta_bias",variable))) %>% 
  # filter(n_subs == 25) %>% 
  filter(div == 0, ess_bulk > 400, ess_tail > 400, rhat < 1.03, tree == 0) %>% 
  ggplot(aes(x = simulated, y = mean,ymin = q5,ymax = q95, group = draw_id))+
  geom_pointrange(alpha = 0.25)+
  facet_wrap(~variable, scales = "free")+
  theme_classic()+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 2))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
  labs(x = "Simulated",y = " ")+
  ggtitle("Subject level parameter recovery")+
  geom_abline(col = "red")



group_param = purrr::map_dfr(
  results_list[1:length(results_list)],
  ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 3))
) %>%   
  filter(grepl("gm",variable)) %>% 
  filter(variable %in% c("gm[8]","gm[9]")) %>% 
  mutate(variable = recode(variable, !!!renames, .default = variable)) %>% 
  mutate(Diagnostics = as.factor(ifelse((div == 0 |  ess_bulk < 400 | ess_tail < 400 | rhat > 1.03),"Good","Bad"))) %>% 
  select(variable,rhat,ess_bulk,ess_tail,n_subs,sim_id,draw_id,Diagnostics)



draws05 = c(683,420,1599,931,1228)
draws0 = c(932,1728,1077,1084,961)
drawsm05 = c(1664,642,1351,1105,1207)



main_parameters = purrr::map_dfr(
  results_list[1:length(results_list)],
  ~ dplyr::bind_rows(purrr::map(.x, purrr::pluck, 3))
) %>%   
  filter(grepl("gm",variable)) %>% 
  filter(variable %in% c("gm[8]","gm[9]")) %>% 
  mutate(variable = recode(variable, !!!renames, .default = variable)) %>% 
  mutate(Diagnostics = as.factor(ifelse((div == 0 |  ess_bulk < 400 | ess_tail < 400 | rhat > 1.03),"Good","Bad"))) %>% 
  mutate(prior_sd = as.factor(prior_sd)) %>% 
  filter(draw_id %in% c(draws05,draws0,drawsm05)) %>%
  mutate(effect = ifelse(draw_id %in% draws0, paste0("effect = 0","\n id = ", draw_id),
  ifelse(draw_id %in% draws05, paste0("effect = 0.05","\n id = ", draw_id),
  ifelse(draw_id %in% drawsm05, paste0("effect = -0.05","\n id = ", draw_id),NA)))) %>%
  ggplot(aes(x = n_subs, y = mean,ymin = q5,ymax = q95, shape = Diagnostics, col = variable))+
  geom_pointrange(position = position_dodge(width = 1), show.legend = F)+
  geom_hline(aes(yintercept = simulated), col = "black",linetype = 2)+
  scale_shape_manual(
    name = "Diagnostics",
    values = c("Bad" = 4,  # filled circle
               "Good" = 1)   # triangle (bad)
  )+
  labs(y = "Parameter value, mean/95HDI", x = "N_subjects")+
  scale_colour_manual(
    name = "Parameter",
    values = c("beta_sigma_m" = "orange",
               "beta_meta_bias" = "darkgreen")
  ) +
  scale_x_continuous(limits = c(0,30), breaks = seq(0,30,by = 10), labels = seq(0,30,by = 10))+
  facet_wrap(~effect, nrow = 3)+
  # facet_wrap(~draw_id, nrow = 4)+
  theme(legend.position = "top", scales = "free")+
  theme_classic()


Sequential_sampling = inner_join(purrr::map_dfr(
  results_list,
  function(x) {
    valid <- purrr::keep(x, ~ !identical(.x, "Error"))
    dplyr::bind_rows(purrr::map(valid, purrr::pluck, 1))
  }
),group_param) %>% rename("beta_sigma_m" = "bf8",
                          "beta_meta_bias" = "bf9") %>% 
  mutate(variable = recode(variable, !!!renames, .default = variable)) %>% 
  pivot_longer(cols = c("beta_sigma_m","beta_meta_bias")) %>% 
  # filter(draw_id != "1530") %>% 
  filter(draw_id %in% c(draws05,draws0,drawsm05)) %>%
  mutate(effect = ifelse(draw_id %in% draws0, paste0("effect = 0","\n id = ", draw_id),
  ifelse(draw_id %in% draws05, paste0("effect = 0.05","\n id = ", draw_id),
  ifelse(draw_id %in% drawsm05, paste0("effect = -0.05","\n id = ", draw_id),NA)))) %>%
  ggplot()+
  geom_line(aes(x = n_subs, y =log(value), group = name, col = name))+
  geom_point(aes(x = n_subs, y = log(value), group = name, shape = Diagnostics))+
  scale_colour_manual(
    name = "Parameter",
    values = c("beta_sigma_m" = "orange",
               "beta_meta_bias" = "darkgreen")
  ) +
  scale_shape_manual(
    name = "Diagnostics",
    values = c("Bad" = 4,  # filled circle
               "Good" = 1)   # triangle (bad)
  ) +
  labs(y = "Log(BF)", x = "N_subjects")+
  facet_wrap(~effect, nrow = 3)+
  # facet_wrap(~draw_id, nrow = 4)+
  geom_hline(yintercept = c(log(1/30),log(30)), linetype = 2)+
  scale_y_continuous(breaks = c(round(log(1/30),2),0,round(log(30),2)))+
  scale_x_continuous(limits = c(0,30), breaks = seq(0,30,by = 10), labels = seq(0,30,by = 10))+
  theme_classic()+
  theme(legend.position = "top")


library(patchwork)

sequential_samping_plot = (main_parameters + 
                             theme(
                               legend.position = "top",
                               legend.justification = "center",
                               legend.text = element_text(size = 18),
                               legend.title = element_text(size = 20),
                               strip.text = element_text(size = 20),
                               axis.title = element_text(size = 18),
                               axis.text = element_text(size = 16, color = "black"),
                               axis.ticks.length = unit(0.2, "cm")
                             ) | 
                             Sequential_sampling + 
                             theme(
                               legend.position = "top",
                               legend.justification = "center",
                               legend.text = element_text(size = 18),
                               legend.title = element_text(size = 20),
                               strip.text = element_text(size = 20),
                               axis.title = element_text(size = 18),
                               axis.text = element_text(size = 16, color = "black"),
                               axis.ticks.length = unit(0.2, "cm")
                             ))


ggsave(here::here("Preregistration","Figures","sequential_samping_plot.PNG"),sequential_samping_plot,dpi = 300, height = 12, width = 20)


parameterrecovery = group_param_recov +                                
  theme(
  legend.position = "top",
  legend.justification = "center",
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 20),
  strip.text = element_text(size = 20),
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 16, color = "black"),
  axis.ticks.length = unit(0.2, "cm")
) | 
  subject_param_recov+  
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )


ggsave(here::here("Preregistration","Figures","parameterrecovery_hier.PNG"),parameterrecovery,dpi = 300, height = 12, width = 20)

