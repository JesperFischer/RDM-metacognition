simulate_trial <- function(D, X, prior_choice, A, prior_correct, delta,conf_prec = 10,cut0 = -4,cut1 = 4,
                           model_type = "independent_biased",
                           mixture_weight = 0.5) {
  
  
  sigma = 1
  # noise
  e <- rnorm(length(D), mean = 0, sd = sigma)
  
  # evidence
  ai <- (X*D-prior_choice) + A * e
  
  # response time
  rt = rnorm(length(D),-abs(ai),0.5)

  a_det = ifelse(ai > 0, 1, -1)
  
  return(data.frame(X = X,
                    D = D,
                    a = a_det,
                    ai = ai,
                    rt = rt))
  # Correct or incorrect?
  correct <- (a == D)
  
  # Asymmetric priors for confidence bias i.e. the prior probability of being correct (according to the agent.)
  prior_correct = prior_correct

  prob_correct = brms::inv_logit_scaled(2*ai/sigma + log(prior_choice / (1-prior_choice)))
  
  
  
  confidence <- prob_correct / (prob_correct + prob_incorrect)
  
  # confidence = ordbetareg::rordbeta(n = 1, mu = confidence,conf_prec,cutpoints = c(cut0,cut1))
  # confidence = rnorm(1,confidence,0.05)
  
  
  return(data.frame(
    correct = correct,
    confidence = confidence,
    action = a_i,
    e = e,
    rt = rt,
    # m = m,
    model_type = model_type
  ))
}

library(tidyverse)

simulate_trial(rep(-1,100),rep(1,100),prior_choice = 0,A = 1,prior_correct = 0.5) %>% 
  ggplot(aes(x = ai,y = a))+
  geom_point()

simulate_trial(rep(-1,100),rep(1,100),prior_choice = 0,A = 1,prior_correct = 0.5) %>% 
  ggplot(aes(x = ai,y = rt))+
  geom_point()




