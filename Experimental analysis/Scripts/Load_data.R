

load_data = function(){
  
  
  paths = list.files(here::here("Data","pilots"), recursive = T, full.names = T)[grepl("\\d\\.csv$",list.files(here::here("Data","pilots"), recursive = T, full.names = T))]
  
  df = data.frame()
  for(path in paths){
    print(path)
    subjectID = strsplit(path, '/')[[1]][[8]]
    
    dd = read.csv(path) %>% mutate(X = NULL) %>% 
      mutate(ID = subjectID) %>% 
        mutate(resp = ifelse(resp == "down",0,1),
               interval = as.numeric(interval),
               D = ifelse(dots.direction == 270,-1,1),
               stim = ifelse(dots.direction == 270, -coherence,coherence)) %>%
      mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
                          ifelse(resp == 0 & stim < 0, 1,0)))
      # filter(Trialtype == "Main") %>%
      # filter(SR_conf != "None") %>%
      # filter(scale == "conf") %>%
      
      if ("sub" %in% names(dd)) {
        names(dd)[names(dd) == "sub"] = "subject"
      }
      
    dd = dd %>% mutate(SR_conf = as.numeric(SR_conf)) %>%
      rename(X = stim, RT = RTdec, Correct = cor, Y = resp, Confidence = SR_conf) %>% 
      select(X,Correct,Y,RT,Confidence,scale,D, coherence,Trialtype,
             interval,interTrial.interval,
             subject,age,gender,ID)
    
  
  df = rbind(df,dd)
  }
  
  return(df)
  
}