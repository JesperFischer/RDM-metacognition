functions {

  // Vectorized psychometric function for accuracy with lapse
  vector psycho_ACC_lapse_vec(vector x, real beta, real alpha, real lapse){
    return lapse + (1-2*lapse) .* Phi(beta .* (x - alpha));
  }
  
  real psycho_ACC_lapse(real x, real beta, real alpha, real lapse){
    return (lapse + (1-2*lapse) * Phi(beta * (x-alpha)));
  }
   
  real psycho_conf(real x, real beta, real alpha){
    return (Phi(beta * abs(x-alpha)));
  }
   
  // ordered beta function - optimized with early returns
  real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {
    real thresh1 = cutzero;
    real thresh2 = cutzero + exp(cutone);

    if(y == 0) {
      return log1m_inv_logit(mu - thresh1);
    } else if(y == 1) {
      return log_inv_logit(mu - thresh2);
    } else {
      real log_mu = log_inv_logit(mu);
      real log1m_mu = log1m_inv_logit(mu);
      return log_diff_exp(log_inv_logit(mu - thresh1), log_inv_logit(mu - thresh2)) +
             beta_lpdf(y | exp(log_mu + log(phi)), exp(log1m_mu + log(phi)));
    }
  }

  real induced_dirichlet_lpdf(real nocut, vector alpha, real phi, int cutnum, real cut1, real cut2) {
    if(cutnum != 1) return 0;
    
    int K = num_elements(alpha);
    vector[K-1] c = [cut1, cut1 + exp(cut2)]';
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    return dirichlet_lpdf(p | alpha) + log_determinant(J) + cut2;
  }

}


data {
  int<lower=0> N;
  int<lower=0> S;
  array[N] int<lower=1,upper=S> S_id;

  array[S] int starts;
  array[S] int ends;


  vector[N] C;
  vector[N] X;
  vector[N] XD;
  vector[N] interval;

  vector[N] ACC;
  array[N] int<lower=0,upper=1> a;
}

transformed data{
  int P = 9;
  vector[3] dirichlet_alpha = [1, 10, 1]';
}

parameters {
  vector[P] gm;
  vector<lower=0>[P] tau_u;
  cholesky_factor_corr[P] L_u;
  matrix[P, S] z_expo;

  vector[S] c0;
  vector[S] c11;
}

transformed parameters{
  // Vectorized parameter computation
  matrix[S, P] param = rep_matrix(gm', S) + (diag_pre_multiply(tau_u, L_u) * z_expo)';

  // Extract transformed parameters once
  vector[S] alpha1 = param[,1];
  vector[S] beta1 = param[,2];
  vector[S] conf_prec1 = param[,3];
  vector[S] meta_un_cor1 = param[,4];
  vector[S] meta_un_inc1 = param[,5];
  vector[S] meta_bias = param[,6];
  vector[S] lapse = param[,7];
  vector[S] meta_un_beta = param[,8];
  vector[S] meta_bias_beta = param[,9];

  // Pre-compute commonly used transformations
  vector[S] exp_beta1 = exp(beta1);
  vector[S] inv_logit_alpha1 = inv_logit(alpha1);
  vector[S] transformed_alpha = (inv_logit_alpha1 - 0.5) * 2;
  vector[S] inv_logit_lapse = inv_logit(lapse) * 0.5;
  vector[S] exp_conf_prec1 = exp(conf_prec1);

  vector[N] theta;
  vector[N] theta_conf;
  vector[N] logit_theta_conf_bias;

  // Vectorized computation where possible
  for (s in 1:S) {
    theta[starts[s]:ends[s]] = psycho_ACC_lapse_vec(XD[starts[s]:ends[s]], exp_beta1[s], transformed_alpha[s], inv_logit_lapse[s]);
  }

  for(n in 1:N){
    int sid = S_id[n];

    if(ACC[n] == 1){
      theta_conf[n] = psycho_conf(XD[n], 
                                   exp(beta1[sid] - exp(meta_un_cor1[sid]) + meta_un_beta[sid] * interval[n]), 
                                   transformed_alpha[sid]);
    } else {
      theta_conf[n] = psycho_conf(XD[n], meta_un_inc1[sid], transformed_alpha[sid]);
    }
    
    logit_theta_conf_bias[n] = logit(theta_conf[n]) + meta_bias[sid] + meta_bias_beta[sid] * interval[n];
  }
}


model {
  // Priors
  gm[1] ~ normal(0, 0.5);
  gm[2] ~ normal(1, 2);
  gm[3] ~ normal(3, 2);
  gm[4] ~ normal(0, 2);
  gm[5] ~ normal(0, 2);
  gm[6] ~ normal(0, 0.5);
  gm[7] ~ normal(-4, 2);
  gm[8] ~ normal(0, 0.5);
  gm[9] ~ normal(0, 0.5);

  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(0, 1);
  tau_u[2] ~ normal(0, 2);
  tau_u[3] ~ normal(0, 2);
  tau_u[4] ~ normal(0, 2);
  tau_u[5] ~ normal(0, 2);
  tau_u[6] ~ normal(0, 1);
  tau_u[7] ~ normal(0, 2);
  tau_u[8] ~ normal(0, 0.5);
  tau_u[9] ~ normal(0, 0.5);
  
  L_u ~ lkj_corr_cholesky(2);

  // Vectorized likelihood for binary outcomes
  target += bernoulli_lpmf(a | theta);

  // Confidence likelihood
  for (n in 1:N) {
    int sid = S_id[n];
    target += ord_beta_reg_lpdf(C[n] | logit_theta_conf_bias[n], exp_conf_prec1[sid], c0[sid], c11[sid]);
  }

  // Cutpoint priors
  for(s in 1:S){
    target += induced_dirichlet_lpdf(c0[s] | dirichlet_alpha, 0, 1, c0[s], c11[s]);
    target += induced_dirichlet_lpdf(c11[s] | dirichlet_alpha, 0, 2, c0[s], c11[s]);
  }
}

generated quantities {
  vector[S] c1 = c0 + exp(c11);
  matrix[P,P] correlation_matrix = L_u * L_u';

  vector[N] log_lik_bin;
  vector[N] log_lik_conf;
  vector[N] log_lik;

  for(n in 1:N){
    int sid = S_id[n];
    log_lik_bin[n] = bernoulli_lpmf(a[n] | theta[n]);
    log_lik_conf[n] = ord_beta_reg_lpdf(C[n] | logit_theta_conf_bias[n], exp_conf_prec1[sid], c0[sid], c11[sid]);
    log_lik[n] = log_lik_bin[n] + log_lik_conf[n];
  }
}
