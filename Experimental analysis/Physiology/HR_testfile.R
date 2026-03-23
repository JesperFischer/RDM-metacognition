pacman::p_load("jsonlite", "tidyverse", "ggplot2", "purrr")

file_path <- "C:/Users/User/Downloads/trials.json"

file_path = "C:\\Users\\User\\OneDrive\\Documenten\\trials_kelly.json"

# Read JSON
trials <- fromJSON(file_path)

ROR = c()
C=c()
# Check correlation reorientation response
for (i in 1:length(trials)){
  HR = trials[[i]][["physio"]][["ECG_Rate"]]
  RO = which.min(HR)
  conf = which(trials[[i]][["physio"]][["trigger"]]==50)
  
  ROR= c(ROR, RO)
  C = c(C, conf)
}

data = data.frame(ROR, C)
result = cor.test(ROR, C)
result
minC= min(C)
maxC= max(C)
# Plot single trials
df_plot <- bind_rows(trials[1:50], .id = "trial") %>%
  group_by(trial) %>%
  mutate(trigger_idx = which(trigger == 50)) %>% mutate(idx = row_number()) %>% unnest()

ggplot(df_plot, aes(x = idx, y = ECG_Rate)) +
  geom_smooth(alpha = 0.4)+ facet_wrap(~trial, scales = "free")+
  geom_vline(aes(xintercept=trigger_idx))

# Averaged signal 1 PP
var = data.frame(trial = seq_along(trials)-1, var= map_chr(trials, ~ .x[["RDM"]][["cor"]]))

data_av = map(trials, "physio") %>% bind_rows(.id = "trial") %>%
  mutate(trial = as.integer(trial))  %>% left_join(var, by="trial") %>% group_by(trial)%>%filter(!is.na(var)) %>% 
  mutate(idx = row_number()) %>% mutate(base = mean(ECG_Rate[idx <= 500], na.rm = TRUE)) %>% filter(idx<=7000) %>% 
  mutate(base_cor = ECG_Rate - base) %>% ungroup() %>% group_by(idx, var) %>%         
  summarise(HR_avg = mean(ECG_Rate, na.rm = TRUE), HR_cor_avg = mean(base_cor, na.rm = TRUE)) %>% ungroup()

trials_data = bind_rows(trials[1:length(trials)], .id = "trial") %>%group_by(trial) %>%
  mutate(idx = row_number()) %>% mutate(base = mean(ECG_Rate[idx <= 1000], na.rm = TRUE)) %>% 
  mutate(base_cor = ECG_Rate - base)

ggplot(data_av)+geom_smooth(aes(x=idx, y= HR_cor_avg))+geom_vline(aes(xintercept = mean(C))) + 
  geom_vline(aes(xintercept = 1000), color = "red", linetype = "dashed")+labs(title = "stimulus locked Jesper")+ 
  annotate("rect",xmin = min(C), xmax = max(C),  ymin = -Inf, ymax = Inf,   alpha = 0.2, fill = "grey") +theme_classic()
  

ggplot(data_av)+geom_smooth(aes(x=idx, y= HR_cor_avg, color = var))+geom_vline(aes(xintercept = mean(C))) + 
  geom_vline(aes(xintercept = 2900), color = "red", linetype = "dashed")+labs(title = "stimulus locked Jesper")+ 
  annotate("rect",xmin = min(C), xmax = max(C),  ymin = -Inf, ymax = Inf,   alpha = 0.2, fill = "grey") +theme_classic()

# Averaged around conf
data_conf = map(trials, "physio") %>% bind_rows(.id = "trial") %>%mutate(trial = as.integer(trial))  %>% 
  left_join(var, by="trial") %>% group_by(trial) %>%  mutate(idx = row_number()) %>% filter(!is.na(var)) %>% 
  mutate(trigger_idx = which(trigger == 50)[1]) %>% filter(idx >= (trigger_idx - 4000) & idx <= (trigger_idx + 7000)) %>% 
  ungroup() %>% select(-idx, -trigger_idx) %>%group_by(trial) %>%  mutate(idx = row_number()) %>% group_by(idx, var) %>% summarise(HR_avg = mean(ECG_Rate, na.rm = TRUE))

ggplot(data_conf)+geom_smooth(aes(x=idx, y= HR_avg, color = var))+geom_vline(aes(xintercept = 4001))+labs(title= "Confidence locked Jesper")+ theme_classic()

ggplot(data_conf)+geom_line(aes(x=idx, y= HR_avg))+geom_vline(aes(xintercept = 4001))+labs(title= "Confidence locked Jesper")+ theme_classic()