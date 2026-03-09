
library(furrr)
library(shiny)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tibble)
pacman::p_load(ordbetareg)


sim_data <- function(df){
  
  
  sigma_e = df$sigma_e
  sigma_k = df$sigma_k
  sigma_m = df$sigma_m
  sigma_x = df$sigma_x
  replicate = df$replicate
  
  n_trials <- 5000
  conf_prec <- 500
  c0 <- -3
  c1 <- 3
  
  # Generate data
  X <- abs(rnorm(n_trials, 10, 5))
  D <- sample(c(-1, 1), n_trials, replace = TRUE)
  e <- rnorm(n_trials, mean = D * X, sd = sigma_e)
  ehat <- rnorm(n_trials, e, sigma_m)
  k_gauss <- rnorm(n_trials, 0, sigma_k)
  p_gauss <- pnorm(e / sigma_k)
  a <- rbinom(n_trials, 1, p_gauss)
  Xhat = abs(rnorm(n_trials, X,sigma_x))
  
  c_mu <- ifelse(a == 1, 
                 plogis(2 * e * Xhat / (sigma_e^2)),
                 1 - plogis(2 * e * Xhat / (sigma_e^2)))

  ACC <- ifelse(a == 0 & D == -1, 1, ifelse(a == 1 & D == 1, 1, 0))
  SS <- D * X
  
  
  df <- data.frame(
    e = e,
    ehat = ehat,
    C = C,
    c_mu = c_mu,
    a = a,
    X = X,
    D = D,
    ACC = ACC,
    SS = SS
  )
  
  # m1 = summary(lm(c_mu~ACC+X+ACC:X, data = df %>% mutate(SS = D*X)))
  # 
  # df = data.frame(m1$coefficients) %>% 
  #   rownames_to_column("variable") %>% 
  #   rename_with(~c("variable","Esti","std","t","p")) %>% 
  #   mutate(sigma_e = sigma_e,
  #          sigma_k = sigma_k,
  #          sigma_m = sigma_m,
  #          sigma_x = sigma_x,
  #          replicate = replicate)
  # 
  return(df)
  
  
  
}

q = sim_data(data.frame(sigma_e = 0.5,
                    sigma_k = 0.25,
                    sigma_m = 0,
                    sigma_x = 0.5,
                    replicate = 1))

to = 5
by = 1

df = expand.grid(sigma_e = seq(0.1,to+0.1,by = by),
                 sigma_k = seq(0,to,by = by),
                 sigma_m = seq(0,to,by = by),
                 sigma_x = seq(0,to,by = by),
                 replicate = 1:20) %>% 
  filter(sigma_e != 0 | sigma_k != 0 | sigma_m != 0 | sigma_x != 0)

df <- split(df, seq(nrow(df)))


plan(multisession, workers = 2)  # Windows-friendly

results_list <- future_map(
  df,
  sim_data,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

library(colorspace)


cols_5x5 <- c(
  # BLUE
  hcl(h = 220, c = c(90), l = 60),
  
  # ORANGE
  hcl(h = 40,  c = c(90), l = 60),
  
  # GREEN
  hcl(h = 140, c = c(90), l = 60),
  
  # PURPLE
  hcl(h = 300, c = c(90), l = 60),
  
  # RED
  hcl(h = 0,   c = c(90), l = 60)
)


bind_rows(results_list) %>% 
  filter(sigma_m != 0 & variable == "(Intercept)") %>% 
  group_by(variable,sigma_e,sigma_k,sigma_m,sigma_x) %>% 
  summarize(mean = mean(t),
            q5 = mean(t) - 2 * sd(t) / sqrt(n()),
            q95 = mean(t) + 2 * sd(t) / sqrt(n())) %>% 
  ggplot(aes(x = variable, y = mean, ymin = q5, ymax = q95, fill = interaction(as.factor(sigma_k))))+
  geom_pointrange(position = position_dodge(width = .4), shape = 21,col = "black")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_color_manual(values = cols_5x5)+
  theme_classic()+
  # scale_y_continuous(limits = c(-35,70))+
  facet_grid(sigma_e~sigma_x, labeller = label_both)


bind_rows(results_list) %>% 
  filter(sigma_m != 0 & variable == "X") %>% 
  group_by(variable,sigma_e,sigma_k,sigma_m,sigma_x) %>% 
  summarize(mean = mean(t),
            q5 = mean(t) - 2 * sd(t) / sqrt(n()),
            q95 = mean(t) + 2 * sd(t) / sqrt(n())) %>% 
  ggplot(aes(x = variable, y = mean, ymin = q5, ymax = q95, fill = interaction(as.factor(sigma_k))))+
  geom_pointrange(position = position_dodge(width = .4), shape = 21,col = "black")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_color_manual(values = cols_5x5)+
  theme_classic()+
  # scale_y_continuous(limits = c(-35,70))+
  facet_grid(sigma_e~sigma_x, labeller = label_both)



bind_rows(results_list) %>% 
  filter(sigma_m != 0 & variable == "ACC") %>% 
  group_by(variable,sigma_e,sigma_k,sigma_m,sigma_x) %>% 
  summarize(mean = mean(t),
            q5 = mean(t) - 2 * sd(t) / sqrt(n()),
            q95 = mean(t) + 2 * sd(t) / sqrt(n())) %>% 
  ggplot(aes(x = variable, y = mean, ymin = q5, ymax = q95, fill = interaction(as.factor(sigma_k))))+
  geom_pointrange(position = position_dodge(width = .4), shape = 21,col = "black")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_color_manual(values = cols_5x5)+
  theme_classic()+
  # scale_y_continuous(limits = c(-35,70))+
  facet_grid(sigma_e~sigma_x, labeller = label_both)



bind_rows(results_list) %>% 
  filter(sigma_m != 0 & variable == "ACC:X") %>% 
  group_by(variable,sigma_e,sigma_k,sigma_m,sigma_x) %>% 
  summarize(mean = mean(t),
            q5 = mean(t) - 2 * sd(t) / sqrt(n()),
            q95 = mean(t) + 2 * sd(t) / sqrt(n())) %>% 
  ggplot(aes(x = variable, y = mean, ymin = q5, ymax = q95, fill = interaction(as.factor(sigma_k))))+
  geom_pointrange(position = position_dodge(width = .4), shape = 21,col = "black")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_color_manual(values = cols_5x5)+
  theme_classic()+
  # scale_y_continuous(limits = c(-35,70))+
  facet_grid(sigma_e~sigma_x, labeller = label_both)


