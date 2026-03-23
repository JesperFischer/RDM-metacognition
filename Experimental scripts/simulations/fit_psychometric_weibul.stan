
data {
  int N;
  array[N] int y;
  vector[N] x;
  real gamma;
}


parameters {
  real a;
  real b;
  real lapse_inv;
}

transformed parameters {
  real alpha = exp(a);
  real beta = exp(b);
  real lapse = inv_logit(lapse_inv);
  vector[N] psi;
  
  for (i in 1:N){
    psi[i] = gamma + (1-lapse-gamma)*(1-exp(-pow(x[i]/alpha, beta)));}
}

model {
  target += normal_lpdf(a|1.4,0.5);
  target += normal_lpdf(b|1.2,0.7);
  target += normal_lpdf(lapse_inv|-3,0.7);

  for (i in 1:N){y[i] ~ bernoulli(psi[i]);}
 
}
