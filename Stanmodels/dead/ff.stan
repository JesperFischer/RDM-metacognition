functions {
  // Posterior probability of being correct given evidence
  real p_correct_given_e(real evi, real mu_c, real mu_ic) {
    real lp_c  = normal_lpdf(evi | mu_c, 1);
    real lp_ic = normal_lpdf(evi | mu_ic, 1);
    return (0.5 * exp(lp_c) / (0.5 * exp(lp_c) + 0.5 * exp(lp_ic)));
  }
  
  real conf_integrand(real evi, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real mu_c = theta[1];
    real mu_ic = theta[2];
    real conf_i = x_r[1];  // pass observed confidence via x_r
  
    // posterior probability given evi
    real pc = p_correct_given_e(evi, mu_c, mu_ic);
  
    // likelihood of observed confidence given posterior
    return exp(normal_lpdf(conf_i | pc, 0.05)) * exp(normal_lpdf(evi | mu_c, 1)); 
    // multiply by density of evidence
  }
}


data {
  int<lower=1> N;

  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=-1, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  array[N] real C;    // reported confidence
}

transformed data {
  array[0] int x_i;
  array[N] real x_r;
  x_r = C;  // observed confidence
    
}


parameters {
  real alpha;              // decision bias
  real beta;               // sensitivity (on evidence scale)
  // real<lower=-0.5, upper 0.5> delta;              // confidence prior bias
  // real<lower=0> phi_conf;  // confidence precision
}
model {
  alpha ~ normal(0, 1);
  beta  ~ normal(1, 1);
  // phi_conf ~ normal(0, 1);
  for (i in 1:N) {
    // Choice likelihood (psychometric)
    real p_choice = Phi(exp(beta) * (D[i] * X[i] - alpha));
    target += bernoulli_lpmf((a[i] + 1) / 2 | p_choice);

    // Confidence likelihood (integrated over evidence)
    array[2] real theta;
    theta[1] = exp(beta) * (a[i]*X[i] - alpha);   // mu_c
    theta[2] = exp(beta) * (-a[i]*X[i] - alpha);  // mu_ic

    // print(exp(integrate_1d(conf_integrand, -6, 6, theta, x_r, x_i)));
    target += (
      log(integrate_1d(conf_integrand, -6, 6, theta, x_r, x_i))
    );
  }
}


generated quantities{
  real beta1 = exp(beta);
}
