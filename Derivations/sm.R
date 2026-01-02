# Load libraries
library(ggplot2)

# -----------------------------
# 1. Set parameters
# -----------------------------
n_trials <- 500        # number of trials
sigma_e <- 1             # sensory noise
mu_beta <- 0.5               # Beta criterion alpha
prec_beta <- 1.8                # Beta criterion beta
mu_k <- 0.2
sigma_k <- 1           # Gaussian threshold SD
rt_ndt = 0.2
rt_sigma = 0.5
rt_e = 0.3

X <- abs(rnorm(n_trials,1,2))


# -----------------------------
# 2. Simulate stimuli
# -----------------------------
D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus

# -----------------------------
# 3. Simulate evidence
# -----------------------------
e <- rnorm(n_trials, mean = D * X, sd = sigma_e)

# -----------------------------
# 4. Compute posterior P(D=1|e)
# -----------------------------
p_D1 <- 1 / (1 + exp(-2 * X * e / sigma_e^2))

# -----------------------------
# 5a. Posterior-based stochastic policy
# -----------------------------
c_beta <- extraDistr::rprop(n_trials,prec_beta,mu_beta)           # random threshold
a_beta <- ifelse(p_D1 > c_beta, 1, -1)         # choices
# -----------------------------
# 5b. Evidence-based stochastic policy
# -----------------------------
k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # random threshold on e
a1_gauss <- ifelse(e > k_gauss, 1, -1)          # choices
p_gauss = pnorm((e-mu_k)/sigma_k)
a_gauss = rbinom(n_trials,1,p_gauss)

# -----------------------------
# 6. Plot induced psychometric curves
# -----------------------------
df <- data.frame(
  e = e,
  D = D,
  X =X,
  rt = rlnorm(n_trials,rt_e * -abs(e),rt_sigma),
  rt2 = rlnorm(n_trials,-abs(p_D1-0.5),rt_sigma) - rt_ndt,
  a_beta = a_beta,
  a_gauss = a_gauss
)

# Plot
plot_df <- data.frame(
  e = rep(x_bin, 2),
  p_choice = c(p_beta, p_gauss),
  policy = rep(c("Posterior Beta", "Evidence Gaussian"), each = length(x_bin))
)



df %>% pivot_longer(cols = c("a_gauss","a_beta")) %>% mutate(value = ifelse(value < 0, 0, value)) %>% 
  ggplot(aes(x = e, y = value, col = name)) +
  geom_point()+
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE)+
  labs(x = "Evidence (e)", y = "P(choice = 1)",
       title = "Psychometric curves under different stochastic policies") +
  theme_minimal() +
  theme(text = element_text(size = 14))


df %>%  
  ggplot(aes(x = e, y = rt)) +
  geom_point()+
  labs(x = "Evidence (e)", y = "P(choice = 1)",
       title = "Psychometric curves under different stochastic policies") +
  theme_minimal() +
  theme(text = element_text(size = 14))




mod = cmdstanr::cmdstan_model(here::here("Derivations","test.stan"))

datastan = list(a = df$a_gauss,
                RT = df$rt,
                D = df$D,
                N = nrow(df),
                X = df$X)


fit = mod$sample(data = datastan,
                 iter_warmup = 500,
                 iter_sampling = 500,
           parallel_chains = 4)

fit
dia = as_draws_df(fit$sampler_diagnostics())
params = as_draws_df(fit$draws(c("sigma_choice","sigma_rt")))
plot(dia$energy__, params$sigma_rt)
# plot(dia$energy__, params$sigma_e)
plot(dia$energy__, params$sigma_choice)

mcmc_pairs(params)

data.frame(fit$summary("evidence")) %>% mutate(simulated = df$e) %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5, ymax = q95))+
  geom_pointrange()

qq = data.frame(fit$summary("evidence"))

qq %>% mutate(dif = mean - df$e) %>% ggplot(aes(x = mean, y = ess_bulk))+geom_point()

qq %>% mutate(dif = mean - df$e, e = df$e) %>% ggplot(aes(x = e, y = ess_bulk))+geom_point()

qq %>% mutate(dif = mean - df$e, e = df$e) %>% ggplot(aes(x = dif, y = ess_bulk))+geom_point()


fit$summary(c("RT_e","sigma_choice","sigma_rt","mean_choice"))
