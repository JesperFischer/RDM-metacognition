library(tidyverse)
library(here)
library(cmdstanr)
library(ggplot2)
library(ordbetareg)
library(posterior)

###############################################################################
######## load / manipulate data ###############################################
###############################################################################

subfolders = list.dirs(here("Data", "Pilots"), recursive = FALSE)

data <- lapply(subfolders, function(f) {
  csv_file <- list.files(f,pattern = "RDM_reportz.*\\.csv$",full.names = TRUE,ignore.case = TRUE)
  if (length(csv_file) == 0) return(NULL)
  read.csv(csv_file[[1]])
})

data <- map(data, ~ .x %>%mutate(
                SR_conf = as.numeric(SR_conf),
                SR_conf = if_else(SR_conf < 0, (SR_conf + 1) / 2, SR_conf),
                up = ifelse(resp == "up", 1, 0),
                cohe) %>% 
              filter(scale == "conf") %>% filter(Trialtype == "Main") %>% filter(!is.na(SR_conf))
              
)



