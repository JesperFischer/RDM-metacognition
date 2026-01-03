
## in progress
qq = function(){

library(tidyverse)
source("C:/Users/au645332/Desktop/RDM metacognition/Coverter.R")
convert_edf()

folder = here::here("data","Siebe_full_v1")

asc_lines <- readLines(grep(".asc",list.files(folder, full.names = T), value = T))

# View some sample lines
samples <- grep("^\\d+", asc_lines, value = TRUE)  # lines starting with timestamps
messages <- grep("^MSG", asc_lines, value = TRUE)  # all event messages


trials_start <- grep("start_trialID", messages, value = TRUE)  # all event messages

stim_start <- grep("start_stimulus", messages, value = TRUE)  # all event messages

stim_end <- grep("end_stimulus", messages, value = TRUE)  # all event messages

conf_start = grep("start_confidence",messages, value = TRUE)

conf_end = grep("end_confidence",messages, value = TRUE)

trial_end = grep("end_trialID",messages, value = TRUE)



trial_starts = match(asc_lines, trials_start)
stim_start = match(asc_lines, stim_start)
stim_end = match(asc_lines, stim_end)
conf_start = match(asc_lines, conf_start)
conf_end = match(asc_lines, conf_end)
trial_end = match(asc_lines, trial_end)

idx_start = which(!is.na(trial_starts))
idx_stim_start = which(!is.na(stim_start))
idx_stim_end = which(!is.na(stim_end))
idx_conf_start = which(!is.na(conf_start))
idx_conf_end = which(!is.na(conf_end))
idx_end = which(!is.na(trial_end))




before_preprocessing = data.frame(T = asc_lines) %>% 
  separate(T, into = c("time", "x", "y", "pupil"), sep = "\t", convert = TRUE, fill = "right") %>% 
  mutate(n = 1:n()) %>% filter(n > 65) %>%  
  mutate(pupil = as.numeric(pupil)) %>% 
  mutate(flag = ifelse(pupil > 2000,1,0)) %>% 
  drop_na()

before_preprocessing_plot = before_preprocessing %>% 
  ggplot(aes(x = n, y = pupil))+
  geom_point(alpha = 0.05)+
  facet_wrap(~flag, scales = "free", ncol = 1)

before_preprocessing_plot


blink_index = 500

data = data.frame()
for(i in 1:length(idx_start)){
  print(i)
  tjek = data.frame(T = asc_lines[idx_stim_start[i]:idx_end[i]]) %>% 
    separate(T, into = c("time", "x", "y", "pupil"), sep = "\t", convert = TRUE, fill = "right")
  
  # skip trial if there is a blink of over 500 ms
  if(nrow(tjek %>% filter(pupil == "    0.0")) > blink_index){
    print(paste0("@@@@@@@@@@@@@@@@@@@@@@@@, Big Blink detected in ", folder, " at trial = ", i))
    next
  }
  
  ref = data.frame(T = asc_lines[(idx_start[i]):idx_stim_start[i]]) %>%
    mutate(
      num_fields = lengths(strsplit(T, "\t")),
      is_valid_row = num_fields == 5
    ) %>%
    separate(T, into = c("time", "x", "y", "pupil"), sep = "\t", convert = TRUE, fill = "right")  %>% 
    filter(is_valid_row) %>% mutate(pupil = as.numeric(pupil),
                                    time = as.numeric(time)) %>% 
    mutate(pupil = pad_and_interpolate_blinks(time, pupil, pad_ms = 100, sample_rate = 1000)) %>% 
    mutate(
      real_row = NA,
      conf_start = NA,
      stim_start = NA,
      conf_end = NA,
      stim_end = NA, 
      trial = i,
      pre_stim = T
    )
  
  ref_sum = ref %>% 
    summarise(mean = mean(pupil), sd = sd(pupil))
  
  
  bb <- data.frame(T = asc_lines[idx_start[i]:idx_end[i]]) %>%
    mutate(
      num_fields = lengths(strsplit(T, "\t")),
      is_valid_row = num_fields == 5
    ) %>%
    separate(T, into = c("time", "x", "y", "pupil"), sep = "\t", convert = TRUE, fill = "right")  %>% 
    mutate(real_row = idx_start[i]:idx_end[i]) %>% 
    mutate(conf_start = real_row == idx_conf_start[i],
           stim_start = real_row == idx_stim_start[i],
           conf_end = real_row == idx_conf_end[i],
           stim_end = real_row == idx_stim_end[i]) %>% 
    filter(is_valid_row | conf_start | conf_end | stim_end |stim_start)%>% mutate(pupil = as.numeric(pupil),
                                    time = as.numeric(time)) %>% 
    mutate(pupil = pad_and_interpolate_blinks(time, pupil, pad_ms = 100, sample_rate = 1000)) %>% 
    mutate(
      trial = i,
      pre_stim = F
    )
    
    
    
    
  bb = rbind(ref,bb)
  bb$ref_mean = ref_sum$mean
  bb$ref_sd = ref_sum$sd
  
  data = rbind(data,bb)
  
  
}

testdata = read.csv(grep("data.csv",list.files(folder,recursive = T,full.names = T),value = T)) %>% filter(Trialtype == "Main") %>% 
  mutate(trial = 1:n())

qq = inner_join(data,testdata)%>% 
  mutate(n = 1:n()) %>% 
  mutate(pupil = as.numeric(pupil))
  # mutate(pupil = interpolate_artifacts(time,pupil))

# write.csv(qq,here::here("full_siebe_data.csv"))

qq <- read_csv("data/Siebe_full_v1/full_siebe_data.csv")


# qq <- qq[seq(1, nrow(qq), by = 100), ]
# write.csv(qq,here::here("full_siebe_data_downsampled.csv"))


after_preprocessing_plot <- qq %>% 
  filter(trial < 10) %>% 
  mutate(
    n = row_number(),
    pupil = as.numeric(pupil)
  ) %>%
  ggplot() +
  geom_point(aes(x = n, y = pupil), col = "black", alpha = 0.05) +
  geom_vline(
    data = . %>% filter(conf_start == TRUE),
    linewidth = 2,
    aes(xintercept = n),
    color = "red"
  ) +
  geom_vline(
    data = . %>% filter(conf_end == TRUE),
    linewidth = 2,
    aes(xintercept = n),
    color = "darkgreen"
  ) +
  geom_vline(
    data = . %>% filter(stim_end == TRUE),
    linewidth = 2,
    aes(xintercept = n),
    color = "black"
  ) +
  geom_vline(
    data = . %>% filter(stim_start == TRUE),
    linewidth = 2,
    aes(xintercept = n),
    color = "yellow"
  ) +
  facet_wrap(~trial, scales = "free")

after_preprocessing_plot


qq   %>% 
  mutate(coherence_bin = cut(coherence,5)) %>% 
  group_by(trial) %>% 
  mutate(
  n = 1:n(),
  pupil = as.numeric(pupil)) %>%
  filter(n < 50) %>% 
  group_by(n,dots.direction,coherence_bin) %>% 
  summarize(mean = mean(pupil),
            se = sd(pupil)/sqrt(n())) %>% 
  ggplot() +
  geom_pointrange(aes(x = n, y = mean,ymin = mean-2*se,ymax = mean+2*se, col = as.factor(dots.direction)))+
    facet_wrap(~coherence_bin)
  
}


