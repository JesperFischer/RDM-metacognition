functions {

real lambda(real z){
  if (z > -5.0) {
    return exp(std_normal_lpdf(z)) / Phi(z);
  } else {
    // Asymptotic expansion for z << 0
    return -1.0 / z;
  }
}
  
  real psycho_ACC(real x, real beta, real alpha, real lapse){
    return (lapse + (1-2*lapse) *Phi(beta * (x-alpha)));
   }
   
  real psycho_conf(real x, real beta, real alpha, real lapse){
    return (lapse + (1-2*lapse) * Phi(beta * abs(x-alpha)));
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
  vector[N] ACC; // Vector of deltaBPM values that match the binary response
  array[N] int a;


}

transformed data{
  int P = 9;
  real Cc = 2/1.7;

}

parameters {
  vector[P] gm;
  vector<lower=0>[P] tau_u;
  cholesky_factor_corr[P] L_u;    // Between participant cholesky decomposition
  matrix[P, S] z_expo;    // Participant deviation from the group means


  vector[S] c0;
  vector[S] c11;

}

transformed parameters{

  matrix[S, P] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';

  matrix[S, P] param;

  for(p in 1:P){
    param[,p]= gm[p] + indi_dif[,p];
  }

  vector[S]  beta = (param[,1]);
  vector[S]  sigma_e = (param[,2]);
  vector[S]  sigma_k = param[,3];  // confidence precision
  vector[S]  sigma_m = (param[,4]);  // meta uncertainty on correct trials
  vector[S]  meta_bias = param[,5];  // meta uncertainty on incorrect trials
  vector[S]  lapse = param[,6];  // meta uncertainty on incorrect trials
  
  vector[S]  confprec = param[,7];  // meta uncertainty on incorrect trials
  vector[S]  sigmam_beta = param[,8];  // meta uncertainty on incorrect trials
  vector[S]  meta_bias_beta = param[,9];  // meta uncertainty on incorrect trials


  vector[S] sigma1 = exp(sigma_e)^2 + exp(sigma_k)^2;  // Now squaring to get variances


 
  vector[N] theta;
  vector[N] theta_conf;

  for (n in 1:N) {
    
  real sigma2 = exp(sigma_e[S_id[n]])^2 + exp(sigma_m[S_id[n]]  +  sigmam_beta[S_id[n]] * interval[n])^2;
  
  real mu = (XD[n] - (inv_logit(beta[S_id[n]])-0.5)*2);

    
  theta[n] = (inv_logit(lapse[S_id[n]])/2) + (1-2*inv_logit(lapse[S_id[n]])/2) * Phi( mu / sqrt(sigma1[S_id[n]]));

  if(a[n] == 1){
    theta_conf[n] = Phi( (((mu + exp(sigma_e[S_id[n]])^2/sqrt(sigma1[S_id[n]]) * lambda(mu/sqrt(sigma1[S_id[n]]))) * Cc*X[n]) / sigma2) * 1/(sqrt(1 + ((Cc^2*X[n]^2) / sigma2))));
    
  }else if(a[n] == 0){
    theta_conf[n] =  1-Phi( (((mu - exp(sigma_e[S_id[n]])^2/sqrt(sigma1[S_id[n]]) * lambda(-mu/sqrt(sigma1[S_id[n]]))) * Cc*X[n]) / sigma2) * 1/(sqrt(1 + ((Cc^2*X[n]^2) / sigma2))));
  }
  


  
  }
}


model {
  
  gm[1] ~ normal(0,0.5); //global mean of threshold 
  gm[2] ~ normal(-2,2); //global mean of slope
  gm[3] ~ normal(-2,2); //global mean of confidence precision
  gm[4] ~ normal(-2,2); //global mean of meta uncertainty
  gm[5] ~ normal(0,0.5); //global mean of meta bias
  gm[6] ~ normal(-4,2); //global mean of meta bias
  
  gm[7] ~ normal(2,2); //global mean of meta bias
  gm[8] ~ normal(0,0.5); //global mean of meta bias
  gm[9] ~ normal(0,0.5); 
  
  to_vector(z_expo) ~ std_normal();

  tau_u[1:7] ~ normal(1 , 1);
  tau_u[8] ~ normal(0 , 0.5);
  tau_u[9] ~ normal(0 , 0.5);
  
  L_u ~ lkj_corr_cholesky(2);

  target += bernoulli_lpmf(a | theta);   // likelihood for the outcomes


  for (n in 1:N) {

    target += ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]) + meta_bias[S_id[n]]+  meta_bias_beta[S_id[n]] * interval[n], exp(confprec[S_id[n]]), c0[S_id[n]], c11[S_id[n]]);   // likelihood for confidence on correct trials 

  }
  
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
