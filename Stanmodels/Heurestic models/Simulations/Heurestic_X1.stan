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

  vector[N] C;

  vector[N] X;
  vector[N] XD;

  vector[N] ACC; // Vector of deltaBPM values that match the binary response
  array[N] int a;


}

transformed data{
  int P = 5;
  real Cc = 2/1.7;

}

parameters {
  vector[P] gm;

  real c0;
  real c11;


}

transformed parameters{




  real beta = gm[1];
  real sigma_e = gm[2];
  real sigma_k = gm[3];
  real sigma_m = gm[4];

  real confprec_p = gm[5];

  real sigma1 = exp(sigma_e)^2 + exp(sigma_k)^2;  // Now squaring to get variances
  real sigma2 = exp(sigma_e)^2 + exp(sigma_m)^2;  // Now squaring to get variances


 
  vector[N] theta;
  vector[N] theta_conf;

  for (n in 1:N) {
  real mu = (XD[n] - (inv_logit(beta)-0.5)*2);

    
  theta[n] = Phi( mu / sqrt(sigma1));

  if(a[n] == 1){
    theta_conf[n] = Phi( (((mu + exp(sigma_e)^2/sqrt(sigma1) * lambda(mu/sqrt(sigma1))) * Cc*1) / sigma2) * 1/(sqrt(1 + ((Cc^2*1) / sigma2))));
    
  }else if(a[n] == 0){
    theta_conf[n] =  1-Phi( (((mu - exp(sigma_e)^2/sqrt(sigma1) * lambda(-mu/sqrt(sigma1))) * Cc*1) / sigma2) * 1/(sqrt(1 + ((Cc^2*1) / sigma2))));
  }
  


  
  }
}


model {
  
  

  gm[1] ~ normal(0,0.5); //global mean of threshold 

  gm[2] ~ normal(-1,1); //global mean of beta

  gm[3] ~ normal(-2,1); //global mean of confprec

  gm[4] ~ normal(-1,1); //global mean of meta uncertainty

  gm[5] ~ normal(2,2); //global mean of confprec



  c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
  c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);

  


  for (n in 1:N) {
    
    target += binomial_lpmf(a[n] | 1, theta[n]);   // likelihood for the outcomes

    target += ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]), exp(confprec_p), c0, c11);   // likelihood for confidence on correct trials

  }


}

generated quantities {

  
  vector[N] log_lik_bin;
  vector[N] log_lik_conf;
  vector[N] log_lik;

  for(n in 1:N){
    
    log_lik_bin[n] = binomial_lpmf(a[n] | 1, theta[n]);
    
    log_lik_conf[n] = ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]), exp(confprec_p), c0, c11);   // likelihood for confidence on correct trials

    log_lik[n] = log_lik_bin[n]+log_lik_conf[n];
  }



}
