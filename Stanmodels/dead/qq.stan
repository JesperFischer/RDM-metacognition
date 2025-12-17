functions {
  // Posterior probability of being correct given evidence
  real p_correct_given_e(real evi, real mu_c, real mu_ic) {
    real lp_c  = normal_lpdf(evi | mu_c, 1);
    real lp_ic = normal_lpdf(evi | mu_ic, 1);
    return (0.5 * exp(lp_c) / (0.5 * exp(lp_c) + 0.5 * exp(lp_ic)));
  }

  // Integrand for confidence likelihood
  real conf_integrand(real evi, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    // Unpack theta
    real mu_c = theta[1];
    real mu_ic = theta[2];

    // Posterior correctness given e
    real pc = p_correct_given_e(evi, mu_c, mu_ic);

    // Likelihood contribution (deterministic, e.g., delta-function)
    // Here we treat observed confidence as matching posterior correctness
    return pc;  // tiny variance to approximate delta
  }
}


data {
  int<lower=1> N;

  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=-1, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  vector<lower=0, upper=1>[N] C;    // reported confidence
}


parameters {
  real alpha;              // decision bias
  real beta;               // sensitivity (on evidence scale)
  // real<lower=-0.5, upper 0.5> delta;              // confidence prior bias
  real<lower=0> phi_conf;  // confidence precision
}
model {
  alpha ~ normal(0, 1);
  beta  ~ normal(0, 1);

  for (i in 1:N) {
    // Choice likelihood (psychometric)
    real p_choice = Phi(beta * (D[i] * X[i] - alpha));
    target += bernoulli_lpmf((a[i] + 1) / 2 | p_choice);

    // Confidence likelihood (integrated over evidence)
    array[2] real theta;
    theta[1] = beta * (a[i]*X[i] - alpha);   // mu_c
    theta[2] = beta * (-a[i]*X[i] - alpha);  // mu_ic

    // print(integrate_1d(conf_integrand,-6, 6,theta, rep_array(0.0,0),rep_array(0,0)));
    real pc;
    pc = integrate_1d(conf_integrand,-6, 6,theta, rep_array(0.0,0),rep_array(0,0));
    print(pc);
    target += normal_lpdf(C[i]| pc, 0.1);
  
  }
}
