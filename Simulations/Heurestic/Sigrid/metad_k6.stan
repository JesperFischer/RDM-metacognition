// generated with brms 2.23.0
functions {
  
  // Convert a binary int x from {0, 1} to {-1, 1}
  int to_signed(int x) {
    return 2*x - 1;
  }

  // P(response, confidence | stimulus) given as simplex
  // [P(resp=0, conf=K), .... P(resp=0, conf=1), P(resp=1, conf=1), ... P(resp=1, conf=K)]
  vector metad_normal_pmf(int stimulus, real dprime, real c, real meta_dprime, real meta_c, vector meta_c2_0, vector meta_c2_1) {
    // number of confidence levels
    int K = size(meta_c2_0) + 1;

    // type-1 response probabilities
    real lp_1 = std_normal_lccdf(c - to_signed(stimulus)*dprime/2);
    real lp_0 = std_normal_lcdf(c - to_signed(stimulus)*dprime/2);

    // means of type-2 distributions
    real meta_mu = to_signed(stimulus) * meta_dprime/2;

    vector[K] lp2_1;         // CDFs (response == 1)
    vector[K] lp2_0;         // CDFs (response == 0)
    vector[2*K] log_theta;   // joint (type-1 x type-2) response probabilities

    lp2_1[1] = std_normal_lccdf(meta_c - meta_mu);
    lp2_0[1] = std_normal_lcdf(meta_c - meta_mu);
    for (k in 2:K) {
      lp2_1[k] = std_normal_lccdf(meta_c2_1[k-1] - meta_mu);
      lp2_0[k] = std_normal_lcdf(meta_c2_0[k-1] - meta_mu);

      log_theta[K-k+2] = log_diff_exp(lp2_0[k-1], lp2_0[k]);
      log_theta[K+k-1] = log_diff_exp(lp2_1[k-1], lp2_1[k]);
    }
    log_theta[1] = lp2_0[K];
    log_theta[2*K] = lp2_1[K];

    // weight by P(response|stimulus) and normalize
    log_theta[1:K] += lp_0 - lp2_0[1];
    log_theta[(K+1):(2*K)] += lp_1 - lp2_1[1];
    return log_theta;
}

  real metad__6__normal__absolute__multinomial_lpmf(array[] int Y, real M, real dprime, real c, real z_meta_c2_0_1, real z_meta_c2_0_2, real z_meta_c2_0_3, real z_meta_c2_0_4, real z_meta_c2_0_5, real z_meta_c2_1_1, real z_meta_c2_1_2, real z_meta_c2_1_3, real z_meta_c2_1_4, real z_meta_c2_1_5) {
  int K = 6; // number of confidence levels

  real meta_dprime = M * dprime;
  real meta_c = c;
  vector[K-1] meta_c2_0 = meta_c - cumulative_sum([z_meta_c2_0_1, z_meta_c2_0_2, z_meta_c2_0_3, z_meta_c2_0_4, z_meta_c2_0_5]');
  vector[K-1] meta_c2_1 = meta_c + cumulative_sum([z_meta_c2_1_1, z_meta_c2_1_2, z_meta_c2_1_3, z_meta_c2_1_4, z_meta_c2_1_5]');

  // use multinomial likelihood
  return multinomial_logit_lpmf(Y[1:(2*K)] | metad_normal_pmf(0, dprime, c,
                          meta_dprime, meta_c, meta_c2_0, meta_c2_1)) +
    multinomial_logit_lpmf(Y[(2*K+1):(4*K)] |  metad_normal_pmf(1, dprime, c,
                      meta_dprime, meta_c, meta_c2_0, meta_c2_1));
}
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=2> ncat;  // number of categories
  array[N, ncat] int Y;  // response array
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real dprime;
  real c;
  real<lower=0> metac2zero1diff;
  real<lower=0> metac2zero2diff;
  real<lower=0> metac2zero3diff;
  real<lower=0> metac2zero4diff;
  real<lower=0> metac2zero5diff;
  real<lower=0> metac2one1diff;
  real<lower=0> metac2one2diff;
  real<lower=0> metac2one3diff;
  real<lower=0> metac2one4diff;
  real<lower=0> metac2one5diff;
}
transformed parameters {
  // prior contributions to the log posterior
  real lprior = 0;
  lprior += normal_lpdf(Intercept | 0, 1);
  lprior += normal_lpdf(dprime | 0, 1);
  lprior += normal_lpdf(c | 0, 1);
  lprior += lognormal_lpdf(metac2zero1diff | 0, 1);
  lprior += lognormal_lpdf(metac2zero2diff | 0, 1);
  lprior += lognormal_lpdf(metac2zero3diff | 0, 1);
  lprior += lognormal_lpdf(metac2one1diff | 0, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    mu = exp(mu);
    for (n in 1:N) {
      target += metad__6__normal__absolute__multinomial_lpmf(Y[n] | mu[n], dprime, c, metac2zero1diff, metac2zero2diff, metac2zero3diff, metac2zero4diff, metac2zero5diff, metac2one1diff, metac2one2diff, metac2one3diff, metac2one4diff, metac2one5diff);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
  
  real dp = dprime;
  real mdp = exp(Intercept) * dprime;
  real mratio = exp(Intercept);
  
  
}