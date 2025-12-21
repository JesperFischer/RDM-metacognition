
pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes,
               brms, patchwork, cowplot,ggpubr,flextable)


invisible(lapply(
  list.files(here::here("Analysis", "Scripts"), full.names = TRUE),
  function(f) source(f, echo = FALSE, print.eval = FALSE)
))

df <- read_csv(here::here("Data","Pilots","Kelly_2025_17_12","RDM_reportz_sub3_1.csv")) %>%
  mutate(resp = ifelse(resp == "down",0,1),
         interval = as.numeric(interval),
         stim = ifelse(`dots direction` == 270, -coherence,coherence)) %>%
  mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
                      ifelse(resp == 0 & stim < 0, 1,0))) %>%
  filter(Trialtype == "Main") %>%
  filter(SR_conf != "None") %>%
  filter(scale == "conf") %>%
  mutate(SR_conf = as.numeric(SR_conf)) %>%
  rename(X = stim, RT = RTdec, Correct = cor, Y = resp) %>%
  mutate(Confidence = (SR_conf+1)/2)


kelly = prelim_ss_plot(df, "kelly")

df <- read_csv(here::here("Data","Pilots","Jesper_2025_17_12","RDM_reportz_sub5_1.csv")) %>%
  mutate(resp = ifelse(resp == "down",0,1),
         interval = as.numeric(interval),
       stim = ifelse(`dots direction` == 270, -coherence,coherence)) %>%
  mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
                      ifelse(resp == 0 & stim < 0, 1,0))) %>%
  filter(Trialtype == "Main") %>%
  filter(SR_conf != "None") %>%
    filter(scale == "conf") %>%
  mutate(SR_conf = as.numeric(SR_conf)) %>%
  rename(X = stim, RT = RTdec, Correct = cor, Y = resp) %>%
  mutate(Confidence = (SR_conf+1)/2)

jesper = prelim_ss_plot(df, "jesper")


df <- read_csv(here::here("Data","Pilots","Siebe_2025_03_12","data.csv")) %>%
  mutate(resp = ifelse(resp == "down",0,1),
         interval = as.numeric(interval),
       stim = ifelse(`dots direction` == 270, -coherence,coherence)) %>%
  mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
                      ifelse(resp == 0 & stim < 0, 1,0))) %>%
  filter(Trialtype == "Main") %>%
  filter(SR_conf != "None") %>%
    # filter(scale == "conf") %>%
  mutate(SR_conf = as.numeric(SR_conf)) %>%
  rename(X = stim, RT = RTdec, Correct = cor, Y = resp) %>%
  mutate(Confidence = (SR_conf+1)/2)

siebe = prelim_ss_plot(df, "siebe")



df <- read_csv(here::here("Data","Pilots","Siebe_2025_17_12","RDM_reportz_sub6.csv")) %>%
  mutate(resp = ifelse(resp == "down",0,1),
         interval = as.numeric(interval),
         stim = ifelse(`dots direction` == 270, -coherence,coherence)) %>%
  mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
                      ifelse(resp == 0 & stim < 0, 1,0))) %>%
  filter(Trialtype == "Main") %>%
  filter(SR_conf != "None") %>%
  filter(scale == "conf") %>%
  mutate(SR_conf = as.numeric(SR_conf)) %>%
  rename(X = stim, RT = RTdec, Correct = cor, Y = resp) %>%
  mutate(Confidence = (SR_conf+1)/2)

siebe2 = prelim_ss_plot(df, "siebe2")


siebe[[2]][[1]] / 
  (siebe2[[2]][[1]]+ theme(legend.position = "none")) / 
  (jesper[[2]][[1]]+ theme(legend.position = "none")) / 
  (kelly[[2]][[1]] + theme(legend.position = "none"))



siebe[[1]] / 
  (siebe2[[1]]+ theme(legend.position = "none")) / 
  (jesper[[1]]+ theme(legend.position = "none")) / 
  (kelly[[1]] + theme(legend.position = "none"))



