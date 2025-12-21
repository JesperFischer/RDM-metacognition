data {
  int<lower=1> N;

  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=-1, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  vector[N] RT;                         // stimulus strength
  
  
}


parameters {
  real A;              // decision bias
  vector[N] evidence_std;
  real sigma_choice;              // decision bias
  real sigma_rt;              // decision bias
  
}


transformed parameters{
  vector[N] evidence;
  vector[N] choice;
  
  for(i in 1:N){
    evidence[i] = exp(A) * (D[i] * X[i]) + evidence_std[i];
    if(evidence[i] > 0){
      choice[i] = 1;
    } else if(evidence[i] <= 0){
      choice[i] = 0;
    }
  }
}

model {
  A ~ normal(0, 1);
  sigma_choice ~ normal(0, 1);
  sigma_rt ~ normal(0, 1);
  evidence_std ~ std_normal();
  
  for (i in 1:N) {

    target += normal_lpdf(a[i] | choice[i], exp(sigma_choice));
    target += normal_lpdf(RT[i] | abs(evidence[i]), exp(sigma_rt));


  }


}

