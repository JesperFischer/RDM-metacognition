fit_model_ss = function(df, model,samples){
  
  

  datastan = list(N = nrow(df),
                a = df$Y,
                RT = df$RT,
                D = df$D,
                C = df$Confidence,
                X = df$X,
                minRT = min(df$RT),
                ACC = df$Correct,
                starts = 1,
                ends = nrow(df),
                interval = df$interTrial.interval)


  fit <-model$sample(
    data = datastan,
    refresh = 10,
    iter_sampling = samples,
    iter_warmup = samples,
    adapt_delta = 0.95,
    max_treedepth = 12,
    init  = 0,
    parallel_chains = 4)
  
  if(grepl("MCMC",model1$stan_file())){
    fit$save_object(here::here("Experimental analysis","MCMC","saved models",paste0("fit_",unique(df$ID),".rds")))
  }
  if(grepl("CANSANDRE",model1$stan_file())){
    fit$save_object(here::here("Experimental analysis","CANSANDRE","saved models",paste0("fit_",unique(df$ID),".rds")))
  }
  return(fit)
  
}


dia = function(fit){
  
  parameters = c("sigma_choice","mean_choice","prec_conf","sigma_m","sigma_e")
  sum = fit$summary(c(parameters))
  
  plot = mcmc_pairs(fit$draws(c(parameters)),np = nuts_params(fit))
  plot2 = mcmc_trace(fit$draws(c(parameters)),np = nuts_params(fit))
  
  return(list(sum,plot,plot2, fit$diagnostic_summary()))
  
}