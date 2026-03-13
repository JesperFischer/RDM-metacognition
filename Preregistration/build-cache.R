# build-cache.R
#
# Run this script locally to generate all cached simulation results.
# These cached files are committed to the repo so that CI builds do not
# need heavy dependencies (brms, ordbetareg, tinytex, rstan).
#
# Usage: Rscript build-cache.R
#   (run from the repo root, or from the Preregistration/ directory)

library(here)
library(tidyverse)
library(truncnorm)
library(cowplot)

cat("=== Building cached results for Quarto preregistration ===\n\n")

cache_dir <- here("Preregistration", "cached-results")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Plate notation ---
cat("1/7 Plate notation...\n")
cache_file <- file.path(cache_dir, "platenotation.png")
if (!file.exists(cache_file)) {
  library(tinytex)
  library(magick)
  invisible(tinytex::latexmk(here("Preregistration", "Figure scripts", "Platenotations", "platenotation.tex")))
  img <- image_read_pdf(here("Preregistration", "Figure scripts", "Platenotations", "platenotation.pdf"), density = 600)
  img <- image_trim(img)
  img <- image_resize(img, "1200x")
  image_write(img, cache_file)
  cat("  -> saved\n")
} else {
  cat("  -> already cached\n")
}

# --- 2. ordbetareg simulation ---
cat("2/7 ordbetareg simulation (df_interaction)...\n")
cache_file <- file.path(cache_dir, "df_interaction_ordbeta.rds")
if (!file.exists(cache_file)) {
  library(ordbetareg)

  set.seed(199742)
  n_trials <- 5000
  conf_prec <- 100
  c0 <- -3
  c1 <- 3

  simulate_interaction <- function(sigma_m, sigma_c, sigma_e) {
    X <- sample(c(0.01, 1, 2, 5), n_trials, replace = TRUE)
    D <- sample(c(-1, 1), n_trials, replace = TRUE)
    e <- rnorm(n_trials, mean = D * X, sd = sigma_e)
    ehat <- rnorm(n_trials, e, sigma_m)
    p_gauss <- pnorm(e / sigma_c)
    a <- rbinom(n_trials, 1, p_gauss)
    C_mu <- ifelse(a == 1,
                   plogis(2 * ehat * X / (sigma_e^2 + sigma_m^2)),
                   1 - plogis(2 * ehat * X / (sigma_e^2 + sigma_m^2)))
    c_mu <- pmax(0.001, pmin(0.999, C_mu))
    C <- ordbetareg::rordbeta(n_trials, c_mu, conf_prec, c(c0, c1))
    data.frame(e = e, C = C, C_mu = C_mu, ehat = ehat, a = a,
               X = X, D = D,
               ACC = ifelse(a == 0 & D == -1, 1, ifelse(a == 1 & D == 1, 1, 0)))
  }

  df_interaction <- expand.grid(
    sigma_m = c(0.1, 0.5, 1, 3, 5),
    sigma_c = c(0.1, 0.5, 1, 3, 5),
    sigma_e = c(0.1, 0.5, 1, 3, 5)
  ) %>%
    rowwise() %>%
    mutate(data = list(simulate_interaction(sigma_m, sigma_c, sigma_e))) %>%
    unnest(data)

  saveRDS(df_interaction, cache_file)
  cat("  -> saved\n")
} else {
  cat("  -> already cached\n")
}

# --- 3. Simulation.R component sims ---
cat("3/7 Simulation.R component sims (df_grid, df_A, df_B)...\n")
cache_grid <- file.path(cache_dir, "df_grid_components.rds")
cache_A <- file.path(cache_dir, "df_A_components.rds")
cache_B <- file.path(cache_dir, "df_B_components.rds")
if (!all(file.exists(cache_grid, cache_A, cache_B))) {
  source(here("Preregistration", "scripts", "Simulation.R"))

  n_trials <- 50000
  sigma_e_vals <- c(0.3, 0.4, 0.5)
  sigma_m_vals <- c(0.3, 0.4, 0.5)
  sigma_c_fixed <- 0.5

  df_grid <- expand.grid(se = sigma_e_vals, sm = sigma_m_vals) %>%
    rowwise() %>%
    mutate(data = list(simulate_components(se, sigma_c_fixed, sm, beta = 0, n = 20000))) %>%
    unnest(data)

  sigma_m_vals_A <- c(0.2, 0.5, 1.0, 2.0)
  df_A <- expand.grid(sm = sigma_m_vals_A) %>%
    rowwise() %>%
    mutate(data = list(simulate_components(0.5, 0.5, sm, beta = 0, n = 20000))) %>%
    unnest(data) %>%
    select(-sm)

  sigma_e_vals_B <- c(0.3, 0.5, 1.0)
  sigma_c_vals_B <- c(0.3, 0.5, 1.0)
  df_B <- expand.grid(se = sigma_e_vals_B, sc = sigma_c_vals_B) %>%
    rowwise() %>%
    mutate(data = list(simulate_components(se, sc, sigma_m = 0.2, beta = 0, n = 20000))) %>%
    unnest(data) %>%
    select(-se, -sc)

  saveRDS(df_grid, cache_grid)
  saveRDS(df_A, cache_A)
  saveRDS(df_B, cache_B)
  cat("  -> saved\n")
} else {
  cat("  -> already cached\n")
}

# --- 4. Validation sims 1 (500K trials) ---
cat("4/7 Validation sims 1 (500K trials)...\n")
cache_file <- file.path(cache_dir, "validation_sims_1.rds")
if (!file.exists(cache_file)) {
  if (!exists("simulate_data_p_bias")) {
    source(here("Preregistration", "scripts", "Simulation.R"))
  }

  trials <- 500000
  sigma_e_vals <- c(-2, -1.5)
  sigma_k_vals <- c(-3, -2)
  sigma_m_vals <- c(-1, 0)

  param_grid <- expand.grid(
    sigma_e = sigma_e_vals, sigma_k = sigma_k_vals,
    sigma_m = sigma_m_vals, bias = 0
  ) %>% mutate(
    id = row_number(),
    facet_label = paste0("S_e=", sigma_e, ", S_k=", sigma_k,
                         ", S_m=", sigma_m, ", B=", bias)
  )

  breaks <- seq(-1, 1, by = 0.05)

  df_mech <- param_grid %>%
    rowwise() %>%
    mutate(sims = list(simulate_data_p_bias(500000, sigma_e, sigma_k, sigma_m, bias, seed = 1))) %>%
    ungroup() %>% unnest() %>%
    mutate(XD = X * D,
           bin = cut(XD, breaks = breaks, include.lowest = TRUE, labels = FALSE),
           XD = (breaks[bin] + breaks[bin + 1]) / 2,
           Action = as.factor(a_gauss), Accuracy = as.factor(ACC))

  df_heu <- param_grid %>%
    rowwise() %>%
    mutate(sims = list(simulate_data_margin_bias(500, sigma_e, sigma_k, sigma_m, bias, seed = 1))) %>%
    ungroup() %>% unnest() %>%
    mutate(XD = X * D, bin = NA,
           Action = as.factor(a_gauss), Accuracy = as.factor(ACC))

  saveRDS(list(df_mech = df_mech, df_heu = df_heu, trials = trials), cache_file)
  cat("  -> saved\n")
} else {
  cat("  -> already cached\n")
}

# --- 5. Validation sims 2 (bias variation) ---
cat("5/7 Validation sims 2 (bias variation)...\n")
cache_file <- file.path(cache_dir, "validation_sims_2.rds")
if (!file.exists(cache_file)) {
  if (!exists("simulate_data_p_bias")) {
    source(here("Preregistration", "scripts", "Simulation.R"))
  }

  param_grid <- expand.grid(
    sigma_e = -1.5, sigma_k = -2, sigma_m = -1.2,
    bias = seq(-0.7, 0.7, length.out = 9)
  ) %>% mutate(
    id = row_number(),
    facet_label = paste0("S_e=", sigma_e, ", S_k=", sigma_k,
                         ", S_m=", sigma_m, ", B=", bias)
  )

  breaks <- seq(-1, 1, by = 0.05)

  df_mech <- param_grid %>%
    rowwise() %>%
    mutate(sims = list(simulate_data_p_bias(500000, sigma_e, sigma_k, sigma_m, bias, seed = 1))) %>%
    ungroup() %>% unnest() %>%
    mutate(XD = X * D,
           bin = cut(XD, breaks = breaks, include.lowest = TRUE, labels = FALSE),
           XD = (breaks[bin] + breaks[bin + 1]) / 2,
           Action = as.factor(a_gauss), Accuracy = as.factor(ACC))

  df_heu <- param_grid %>%
    rowwise() %>%
    mutate(sims = list(simulate_data_margin_bias(5000, sigma_e, sigma_k, sigma_m, bias, seed = 1))) %>%
    ungroup() %>% unnest() %>%
    mutate(XD = X * D, bin = NA,
           Action = as.factor(a_gauss), Accuracy = as.factor(ACC))

  saveRDS(list(df_mech = df_mech, df_heu = df_heu), cache_file)
  cat("  -> saved\n")
} else {
  cat("  -> already cached\n")
}

# --- 6. loocomparison ---
cat("6/7 loocomparison...\n")
cache_file <- file.path(cache_dir, "loocomparison.rds")
results_file <- here("Preregistration", "Results", "loocomparison.rds")
if (!file.exists(cache_file)) {
  if (file.exists(results_file)) {
    file.copy(results_file, cache_file)
    cat("  -> copied from Results/\n")
  } else {
    cat("  -> WARNING: loocomparison.rds not found in Results/. You need to generate this from the model fits.\n")
  }
} else {
  cat("  -> already cached\n")
}

# --- 7. Prior predictive + effect sims ---
cat("7/7 Prior predictive + effect sims...\n")
cache_prior <- file.path(cache_dir, "prior_predictive.rds")
cache_effects <- file.path(cache_dir, "effect_sims.rds")
if (!all(file.exists(cache_prior, cache_effects))) {
  source(here("Simulations", "Heurestic", "scripts", "Simulate_mcmc.R"))

  n_sims <- 50
  df_list <- lapply(1:n_sims, function(i) {
    df <- sim_trials(300, seed = i)
    df$sim_id <- i
    df
  })
  df_all <- do.call(rbind, df_list)
  saveRDS(df_all, cache_prior)

  effects <- c(-0.05, 0, 0.05)
  grid <- expand.grid(effect1 = effects, effect2 = effects)
  df_list <- vector("list", nrow(grid))
  for (i in 1:nrow(grid)) {
    df_list[[i]] <- sim_trials_fix(300, 1997, grid$effect1[i], grid$effect2[i])
    df_list[[i]]$sigma_m_beta <- grid$effect1[i]
    df_list[[i]]$meta_bias_beta <- grid$effect2[i]
  }
  df_effects <- do.call(rbind, df_list)
  saveRDS(df_effects, cache_effects)
  cat("  -> saved\n")
} else {
  cat("  -> already cached\n")
}

cat("\n=== Cache build complete ===\n")
cat("Cached files in:", cache_dir, "\n")
cat(paste(" ", list.files(cache_dir), collapse = "\n"), "\n")
