functions {

  real lambda(real z){
      return exp(std_normal_lpdf(z)) / Phi(z);
  }

  // ordered beta function
  real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

  if(y==0) {
      return log1m_inv_logit(mu - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu  - thresh[2]);
    } else {
      return log_diff_exp(log_inv_logit(mu - thresh[1]), log_inv_logit(mu - thresh[2])) +
                beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
    }
  }

  real induced_dirichlet_lpdf(real nocut, vector alpha, real phi, int cutnum, real cut1, real cut2) {
    int K = num_elements(alpha);
    vector[K-1] c = [cut1, cut1 + exp(cut2)]';
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    if(cutnum==1) {

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

    // divide in half for the two cutpoints

    // don't forget the ordered transformation

      return   dirichlet_lpdf(p | alpha)
           + log_determinant(J) + cut2;
    } else {
      return(0);
    }
  }

}


data {
  
  int<lower=0> N;
  int<lower=0> S;
  array[N] int S_id;
  vector[N] interval;
  vector[N] C;
  vector[N] X;
  vector[N] XD;
  vector[N] ACC;
  array[N] int a;


}

transformed data{
  int P = 8;
  real Cc = 2/1.7;

}

parameters {
  vector[P] gm;
  vector<lower=0>[P] tau_u;
  cholesky_factor_corr[P] L_u;    // Between participant cholesky decomposition
  matrix[S, P] param;    // Participant deviation from the group means


  vector[S] c0;
  vector[S] c11;

}

transformed parameters{


  vector[S]  beta = (param[,1]);
  vector[S]  sigma_e = (param[,2]);
  vector[S]  sigma_k = param[,3];  // confidence precision
  vector[S]  sigma_m = (param[,4]);  // meta uncertainty on correct trials
  vector[S]  meta_bias = param[,5];  // meta uncertainty on incorrect trials
  
  vector[S]  confprec = param[,6];  // meta uncertainty on incorrect trials
  vector[S]  sigmam_beta = param[,7];  // meta uncertainty on incorrect trials
  vector[S]  meta_bias_beta = param[,8];  // meta uncertainty on incorrect trials


  vector[S] sigma1 = sqrt(exp(sigma_e)^2 + exp(sigma_k)^2);  // Now squaring to get variances


 
  vector[N] theta;
  vector[N] theta_conf;

  for (n in 1:N) {
    
  real sigma2 = exp(sigma_e[S_id[n]])^2 + exp(sigma_m[S_id[n]]  +  sigmam_beta[S_id[n]] * interval[n])^2;
  
  real mu = (XD[n] - (inv_logit(beta[S_id[n]])-0.5)*2);

  theta[n] = Phi( mu / sigma1[S_id[n]]);

  if(a[n] == 1){
    theta_conf[n] = Phi( (((mu + exp(sigma_e[S_id[n]])^2/sigma1[S_id[n]] * lambda(mu/sigma1[S_id[n]])) * Cc*1) / sigma2) * 1/(sqrt(1 + ((Cc^2*1^2) / sigma2))));
    
  }else if(a[n] == 0){
    theta_conf[n] =  1-Phi( (((mu - exp(sigma_e[S_id[n]])^2/sigma1[S_id[n]] * lambda(-mu/sigma1[S_id[n]])) * Cc*1) / sigma2) * 1/(sqrt(1 + ((Cc^2*1^2) / sigma2))));
  }
  


  
  }
}


model {
  
  gm[1] ~ normal(0,0.5); //global mean of threshold 
  gm[2] ~ normal(-2,2); //global mean of sigma_e
  gm[3] ~ normal(-2,2); //global mean of sigma_k
  gm[4] ~ normal(-2,2); //global mean of sigma_m
  gm[5] ~ normal(0,0.5); //global mean of meta bias
  
  gm[6] ~ normal(2,2); //global mean confidence precision
  gm[7] ~ normal(0,0.5); //global mean of difference in sigma_m from PDW
  gm[8] ~ normal(0,0.5); //global mean of difference in meta_bias from PDW
  
for (s in 1:S){
  param[s,] ~ multi_normal_cholesky(gm, diag_pre_multiply(tau_u, L_u));
}
  // same as above (indicies) but between subject variances

  tau_u[1:6] ~ normal(1 , 1);
  tau_u[7] ~ normal(0 , 0.5);
  tau_u[8] ~ normal(0 , 0.5);
  
  // cholesky decomposition of correlation matrix
  L_u ~ lkj_corr_cholesky(2);

  target += bernoulli_lpmf(a | theta);   // likelihood for the first order choices 

  // likelihood for the second order confidence judgements 
  for (n in 1:N) {
    target += ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]) + meta_bias[S_id[n]]+  meta_bias_beta[S_id[n]] * interval[n], exp(confprec[S_id[n]]), c0[S_id[n]], c11[S_id[n]]);   // likelihood for confidence on correct trials 
  }

  // cutpoints for the ordered beta-likelihood
  for(s in 1:S){
    c0[s] ~ induced_dirichlet([1,10,1]', 0, 1, c0[s], c11[s]);
    c11[s] ~ induced_dirichlet([1,10,1]', 0, 2, c0[s], c11[s]);
  }


}


generated quantities {

  vector[S] c1 = c0 + exp(c11);
  matrix[P,P] correlation_matrix = L_u * L_u';

  
  vector[N] log_lik_bin;
  vector[N] log_lik_conf;
  vector[N] log_lik;

  for(n in 1:N){
    
    log_lik_bin[n] = binomial_lpmf(a[n] | 1, theta[n]);
  
    log_lik_conf[n] = ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]) + meta_bias[S_id[n]]+  meta_bias_beta[S_id[n]] * interval[n], exp(confprec[S_id[n]]), c0[S_id[n]], c11[S_id[n]]);
    log_lik[n] = log_lik_bin[n] + log_lik_conf[n];
  }



}
