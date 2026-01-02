data {
  int<lower=1> N;

  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=0, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  vector[N] RT;                         // stimulus strength
  
  
}


parameters {
  vector[N] evidence_std;
  real mean_choice;
  real sigma_choice_log;              // decision bias
  // real sigma_e;              // decision bias
  real sigma_rt_log;              // decision bias
  real RT_e;              // decision bias
  
  real <lower=0, upper = min(RT)> RT_ndt;
}


transformed parameters{
  vector[N] evidence;
  vector[N] p_action;
  
  for(i in 1:N){
    // evidence[i] = X[i] * D[i] + exp(sigma_e) * evidence_std[i];
    evidence[i] = X[i] * D[i] + evidence_std[i];

    p_action[i] = Phi((evidence[i]-mean_choice)/exp(sigma_choice_log));
  }
}

model {
  sigma_choice_log ~ normal(0, 1);
  // sigma_e ~ normal(0, 1);
  sigma_rt_log ~ normal(0, 1);
  mean_choice ~ normal(0, 1);
  RT_e ~ normal(0, 1);
  
  // RT_ndt ~ normal(0.2, 0.1);
  evidence_std ~ std_normal();
  
  
  a ~ bernoulli( p_action);
  (RT - RT_ndt) ~ lognormal(RT_e * -abs(evidence), exp(sigma_rt_log));
  // (RT) ~ lognormal(RT_e * -abs(evidence), exp(sigma_rt_log));  
  

}

generated quantities{
  
  real sigma_choice = exp(sigma_choice_log);
  real sigma_rt = exp(sigma_rt_log);

}
