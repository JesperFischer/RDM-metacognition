# Parameters from your Python code
intensity_vals <- seq(0.01, 1, length.out = 50)
threshold_vals <- seq(0.01, 1, length.out = 50)  # Your prior has 50 values
slope_vals <- seq(1, 10, length.out = 10)
lower_asymptote <- 0.5  # chance level
lapse_rate <- 0


weibull_psychometric <- function(x, alpha, beta, gamma = 0.5, lambda = 0) {
  
  1-lambda - (1 - gamma - lambda) * (exp(-(x/alpha)^beta))
}

# Sample from uniform priors
set.seed(123)  # for reproducibility
alpha_sample <- sample(threshold_vals, 1)  # or: runif(1, 5, 50)
alpha_sample =0.3
beta_sample <- sample(slope_vals, 1)        # or: runif(1, 1, 10)
beta_sample = 8

# Generate psychometric curve
intensities <- seq(0.01, 1, length.out = 100)
probabilities <- weibull_psychometric(
  x = intensities,
  alpha = alpha_sample,
  beta = beta_sample,
  gamma = lower_asymptote,
  lambda = lapse_rate
)

# Single curve plot
df_single <- data.frame(intensity = intensities, probability = probabilities)

ggplot(df_single, aes(x = intensity, y = probability)) +
  geom_line(linewidth = 1, color = "blue") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
  labs(
    x = "Stimulus Intensity (coherence)",
    y = "P(Correct)",
    title = paste0("Weibull: α=", round(alpha_sample, 2), ", β=", round(beta_sample, 2))
  ) +
  theme_minimal() +
  ylim(0.5, 1)

# Sample multiple parameter sets
n_samples <- 100
alpha_samples <- runif(n_samples, min = 0.01, max = 1)
beta_samples <- runif(n_samples, min = 1, max = 10)

# Generate curves
intensities <- seq(0.01, 1, length.out = 100)
df_multi <- data.frame()

for(i in 1:n_samples) {
  probs <- weibull_psychometric(intensities, alpha_samples[i], beta_samples[i], 0.5, 0)
  df_temp <- data.frame(
    intensity = intensities,
    probability = probs,
    curve_id = i
  )
  df_multi <- rbind(df_multi, df_temp)
}

ggplot(df_multi, aes(x = intensity, y = probability, group = curve_id)) +
  geom_line(alpha = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(
    x = "Stimulus Intensity (coherence)",
    y = "P(Correct)",
    title = "Simulated Psychometric Functions from Prior"
  ) +
  theme_minimal() +
  ylim(0.5, 1)








library(ggplot2)
library(truncnorm)
# Parameters from your Python code
intensity_vals <- seq(0.01, 1, length.out = 50)
threshold_vals <- seq(0.01, 1, length.out = 50)  # Your prior has 50 values
slope_vals <- seq(1, 10, length.out = 10)


# Create non-uniform priors
# Example 1: Normal prior for threshold (centered at 25)
threshold_prior <- dnorm(threshold_vals, mean = 0.3, sd = 0.2)
threshold_prior <- threshold_prior / sum(threshold_prior)  # normalize

# Example 2: Normal prior for slope (centered at 5)
slope_prior <- dnorm(slope_vals, mean = 2, sd = 2)
slope_prior <- slope_prior / sum(slope_prior)  # normalize

# Sample from these priors
n_samples <- 1000
alpha_samples <- sample(threshold_vals, n_samples, replace = TRUE, prob = threshold_prior)
beta_samples <- sample(slope_vals, n_samples, replace = TRUE, prob = slope_prior)

# Generate curves
intensities <- seq(0.01,1, length.out = 100)
df_multi <- data.frame()

for(i in 1:n_samples) {
  probs <- weibull_psychometric(intensities, alpha_samples[i], beta_samples[i], 0.5, 0)
  df_temp <- data.frame(
    intensity = intensities,
    probability = probs,
    curve_id = i
  )
  df_multi <- rbind(df_multi, df_temp)
}

# Plot

ggplot(df_multi, aes(x = intensity, y = probability, group = curve_id)) +
  geom_line(alpha = 0.1, color = "gray50") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(
    x = "Stimulus Intensity (coherence)",
    y = "P(Correct)"
  ) +
  theme_minimal() +
  ylim(0.5, 1)


