
set.seed(199742)

n_trials <- 5000

conf_prec <- 100
c0 <- -3
c1 <- 3

simulate_interaction <- function(sigma_m, sigma_k,sigma_e) {
  X <- sample(c(0.01,1,2,5), n_trials, replace = TRUE)
  D <- sample(c(-1, 1), n_trials, replace = TRUE)
  e <- rnorm(n_trials, mean = D * X, sd = sigma_e)
  ehat <- rnorm(n_trials, e, sigma_m)
  k_gauss <- rnorm(n_trials, 0, sigma_k)
  p_gauss <- pnorm(e / sigma_k)
  a <- rbinom(n_trials, 1, p_gauss)
  C_mu <- ifelse(a == 1, 
                 plogis(2 * ehat * abs(ehat) / (sigma_e^2 + sigma_m^2)),
                 1 - plogis(2 * ehat * abs(ehat) / (sigma_e^2 + sigma_m^2)))
  
  
  c_mu <- pmax(0.001, pmin(0.999, C_mu))
  C <- ordbetareg::rordbeta(n_trials, c_mu, conf_prec, c(c0, c1))
  
  data.frame(e = e,
             C = C,
             C_mu = C_mu,
             # C_mu2 = C_mu2,
             ehat = ehat,
             a = a,
             X = X,
             D=D,
             ACC = ifelse(a == 0 & D == -1,1,ifelse(a == 1 & D == 1,1,0)))
}

df_interaction <- expand.grid(
  sigma_m = c(0.1, 0.5 ,1,3, 5),
  sigma_k = c(0.1, 0.5 ,1,3, 5),
  sigma_e = c(0.1, 0.5 ,1,3, 5)
  
) %>%
  rowwise() %>%
  mutate(data = list(simulate_interaction(sigma_m, sigma_k,sigma_e))) %>% 
  unnest(data)
