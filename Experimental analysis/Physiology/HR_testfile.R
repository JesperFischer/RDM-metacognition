pacman::p_load("jsonlite", "tidyverse", "ggplot2")

file_path <- "C:/Users/User/Downloads/trials.json"

file_path = "C:\\Users\\User\\OneDrive\\Documenten\\trials_kelly.json"

# Read JSON
trials <- fromJSON(file_path)

ROR = c()
C=c()
# CHeck correlation reorientation response
for (i in 1:length(trials)){
  HR = trials[[i]][["ECG_Rate"]]
  RO = which.min(HR)
  conf = which(trials[[i]][["trigger"]]==50)
  
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
data_av = bind_rows(trials[1:length(trials)], .id = "trial") %>%group_by(trial) %>%
  mutate(idx = row_number()) %>% 
  ungroup() %>% filter(idx <= 5000) %>% group_by(idx) %>%         
  summarise(HR_avg = mean(ECG_Rate, na.rm = TRUE))

ggplot(data_av)+geom_smooth(aes(x=idx, y= HR_avg))+geom_vline(aes(xintercept = mean(C)), color = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = 1000))+labs(title = "stimulus locked Kelly")+ 
  annotate("rect",xmin = min(C), xmax = max(C),  ymin = -Inf, ymax = Inf,   alpha = 0.2, fill = "grey") +theme_classic()

# Averaged around conf
data_conf = bind_rows(trials[1:length(trials)], .id = "trial") %>%group_by(trial) %>%  mutate(idx = row_number()) %>% 
  mutate(trigger_idx = which(trigger == 50)[1]) %>% filter(idx >= (trigger_idx - 4000) & idx <= (trigger_idx + 1000)) %>% 
  ungroup() %>% select(-idx, -trigger_idx) %>%group_by(trial) %>%  mutate(idx = row_number()) %>% group_by(idx) %>% summarise(HR_avg = mean(ECG_Rate, na.rm = TRUE))

ggplot(data_conf)+geom_smooth(aes(x=idx, y= HR_avg))+geom_vline(aes(xintercept = 4001))+labs(title= "Confidence locked Jesper")+ theme_classic()