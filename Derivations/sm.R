# Load libraries
library(ggplot2)

# -----------------------------
# 1. Set parameters
# -----------------------------
n_trials <- 500        # number of trials
mu <- 1                  # signal strength
sigma_e <- 1             # sensory noise
mu_beta <- 0.5               # Beta criterion alpha
prec_beta <- 1.8                # Beta criterion beta
mu_k <- 0
sigma_k <- 1           # Gaussian threshold SD
rt_ndt = 0.2
rt_sigma = 0.5


# -----------------------------
# 2. Simulate stimuli
# -----------------------------
D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus

# -----------------------------
# 3. Simulate evidence
# -----------------------------
e <- rnorm(n_trials, mean = D * mu, sd = sigma_e)

# -----------------------------
# 4. Compute posterior P(D=1|e)
# -----------------------------
p_D1 <- 1 / (1 + exp(-2 * mu * e / sigma_e^2))

# -----------------------------
# 5a. Posterior-based stochastic policy
# -----------------------------
c_beta <- extraDistr::rprop(n_trials,prec_beta,mu_beta)           # random threshold
a_beta <- ifelse(p_D1 > c_beta, 1, -1)         # choices
# -----------------------------
# 5b. Evidence-based stochastic policy
# -----------------------------
k_gauss <- rnorm(n_trials, mu_k, sigma_k)          # random threshold on e
a_gauss <- ifelse(e > k_gauss, 1, -1)          # choices

# -----------------------------
# 6. Plot induced psychometric curves
# -----------------------------
df <- data.frame(
  e = e,
  D = D,
  rt = rlnorm(n_trials,-abs(e),rt_sigma) - rt_ndt,
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


df %>%  
  ggplot(aes(x = e, y = rt2)) +
  geom_point()+
  labs(x = "Evidence (e)", y = "P(choice = 1)",
       title = "Psychometric curves under different stochastic policies") +
  theme_minimal() +
  theme(text = element_text(size = 14))

