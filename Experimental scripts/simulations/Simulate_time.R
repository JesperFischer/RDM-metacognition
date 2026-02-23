RT_mean = -0.0811
SD_RT = 0.3859683

mean_conf_RT = 0.1
SD_conf_RT = 0.3697234

n_tutorial = 6
n_training_2 = 10
n_staircase_D = 40
n_SC_coherence = 80

n_t = 36
n_b = 7
  
time_exp = function(){
  intro_T = 10
  RT = exp(rnorm(n_tutorial, RT_mean, SD_RT))
  RT_1 = pmin(RT, 2) + 0.5
  
  tutorial = 10
  
  RT = exp(rnorm(n_staircase_D+n_SC_coherence, RT_mean, SD_RT))
  RT_2 = pmin(RT, 2) + sample(c(0.5,1),1)
  
  instructions = 10*6
  
  RT = exp(rnorm(n_training_2, RT_mean, SD_RT))
  RT_3 = pmin(RT, 2) + 1
  RT_conf =  pmin(exp(rnorm(n_training_2, mean_conf_RT, SD_conf_RT)), 5) + 1
  RT_3 = RT_3 + RT_conf
  
  inst = 10
  
  RT = exp(rnorm(n_t*n_b, RT_mean, SD_RT))
  RT_4 = pmin(RT, 2) + 1
  RT_conf =  pmin(exp(rnorm(n_t*n_b, mean_conf_RT, SD_conf_RT)), 3) + 1 +0.1
  man = runif(n_t*n_b,2,5)
  interT = exp(rnorm(n_t*n_b, 1,0.1))
  
  breakT = 60*n_b + 60+5+5
  
  return( sum(RT_1) + sum(RT_2) + sum(RT_3) + intro_T + tutorial + 
            instructions + sum(RT_4) + inst + sum(RT_conf) + sum(man) + sum(interT) + breakT )
  
}

tot_time = c()
for (i in 1:1000){
  cur = time_exp()
  tot_time = c(tot_time, cur)
}

