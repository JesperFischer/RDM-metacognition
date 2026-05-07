
pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes,
               brms, patchwork, cowplot,ggpubr,flextable)
source(here::here("Sequential sampling", "utility.R"))

df = load_data()

dd = df %>% 
  mutate(subject = as.numeric(as.factor(subject))) %>% 
  filter(Trialtype == "Main", scale == "conf") %>% 
  filter(
    abs(X) < case_when(
      ID == "RDM_reportz_sub001.csv" ~ 0.6,
      ID == "RDM_reportz_sub010.csv" ~ 0.4,
      TRUE ~ Inf
    )) %>% 
  group_by(subject) %>% 
  mutate(
    trial_n = 1:n(), n_trials = n()
  ) %>% 
  ungroup() %>% 
  mutate(
    across(-c("ID","subject"),
           ~ ifelse(
             ID == "RDM_reportz_sub030.csv" & trial_n > (n_trials - 36),
             NA,
             .
           )
    )
  ) 

exclusions_correct = dd %>% group_by(ID) %>% summarize(n = sum(Correct)/n()) %>% filter(n < 0.60)
exclusions_ntrials = dd %>% 
  select(interTrial.interval, interval,coherence,X,Correct,Y,RT,Confidence,ID) %>% 
  mutate(
    no_response = is.na(Correct) |
      is.na(Confidence) |
      is.na(RT) |
      RT < 0.150) %>% 
  group_by(ID) %>%
  summarize(noresp_rate = mean(no_response)) %>%
  filter(noresp_rate > 0.3)

n_excluded_correctness = length(unique(exclusions_correct$ID))
n_excluded_trials = length(unique(exclusions_ntrials$ID))
n_excluded = length(unique(c(exclusions_correct$ID, exclusions_ntrials$ID)))
dd = dd %>% filter(!ID %in% c(exclusions_correct$ID, exclusions_ntrials$ID)) %>% drop_na()


fit_and_save_model = function(dd, model,name){

  samples = 1000
  
  datastan = list(N = nrow(dd),
                  S = length(unique(dd$ID)),
                  S_id = as.numeric(as.factor(dd$ID)),
                  a = dd$Y,
                  a_vec = dd$Y,
                  RT = dd$RT,
                  D = dd$D,
                  C = dd$Confidence,
                  X = dd$X * dd$D,
                  XD = dd$X,
                  minRT = min(dd$RT),
                  ACC = dd$Correct,
                  interval = dd$interval)
  
  # pf = model$pathfinder(data = datastan, psis_resample = F, calculate_lp = F)
  
  
  fit <-model$sample(
    data = datastan,
    refresh = 100,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.99,
    init = 0,
    max_treedepth = 12,
    parallel_chains = 4)
  
  fit$save_object(here::here("Sequential sampling","Fits","Hierarchical","Hierarchical_nolapse_",length(unique(dd$ID)),"_clamped_",name,".rds"))
  
  
}

model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Hierarchical_nolapse_clamped_standard.stan"))
n_ids <- c(15, 20, 25, 30, 35, 40, 45, 50)
all_ids <- unique(dd$ID)

for (n in n_ids) {
  
  ids_use <- all_ids[1:n]   # first n IDs
  
  dd_sub <- dd %>%
    filter(ID %in% ids_use)
  
  fit_and_save_model(
    dd = dd_sub,
    model = model,
    name = "standard"
  )
}



model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Hierarchical_nolapse_clamped_wide.stan"))
n_ids <- c(15, 20, 25, 30, 35, 40, 45, 50)
all_ids <- unique(dd$ID)

for (n in n_ids) {
  
  ids_use <- all_ids[1:n]   # first n IDs
  
  dd_sub <- dd %>%
    filter(ID %in% ids_use)
  
  fit_and_save_model(
    dd = dd_sub,
    model = model,
    name = "wide"
  )
}


model = cmdstanr::cmdstan_model(here::here("Sequential sampling","Stanmodels","Hierarchical_nolapse_clamped_narrow.stan"))
n_ids <- c(15, 20, 25, 30, 35, 40, 45, 50)
all_ids <- unique(dd$ID)

for (n in n_ids) {
  
  ids_use <- all_ids[1:n]   # first n IDs
  
  dd_sub <- dd %>%
    filter(ID %in% ids_use)
  
  fit_and_save_model(
    dd = dd_sub,
    model = model,
    name = "narrow"
  )
}
