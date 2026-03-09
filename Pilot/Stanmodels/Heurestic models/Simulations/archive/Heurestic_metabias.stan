functions {
  
  real psycho_ACC(real x, real beta, real alpha, real lapse){
    return (lapse + (1-2*lapse) *Phi(beta * (x-alpha)));
   }
   
   
  real psycho_conf(real x, real beta, real alpha){
    return (Phi(beta * abs(x-alpha)));
   }
   
  //real entropy(real p){
    //return(-p * log(p) - (1-p) * log(1-p));}       no reaction times for now


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
  
  int<lower=0> N;    // integer: amount of trials
  array[N] int a;    // array: answers
  vector[N] C;      // vector: confidence ratings
  vector[N] X;      // vector: stimulus values
  vector[N] ACC; // vector: outcome binary response

}

transformed data{
  int P = 7;    // Number of subject-level parameters
}

parameters {
  vector[P] gm;   // vector: subject-level parameters

  real c0;        // confidence cut-offs
  real c11;

}

transformed parameters{

  real alpha = gm[1]; // threshold
  real beta = gm[2]; // slope
  real lapse = gm[3]; // slope
  real confprec = gm[4];  // confidence precision
  real metauncor = gm[5];  // meta uncertainty on correct trials
  real metauninc = gm[6];  // meta uncertainty on incorrect trials
  real metabias = gm[7];


  vector[N] conf_mu;
  vector[N] theta;
  vector[N] theta_conf;

  for (n in 1:N) {
  theta[n] = psycho_ACC(X[n], exp(beta), (inv_logit(alpha)-0.5)*2, inv_logit(lapse) / 2) ;

  if(ACC[n] == 1){
    theta_conf[n] = psycho_conf(X[n], exp(beta - exp(metauncor)), (inv_logit(alpha)-0.5)*2);
  }else if(ACC[n] == 0){
    theta_conf[n] = psycho_conf(X[n], metauninc, (inv_logit(alpha)-0.5)*2);
  }
  
  }
  

}
model {
  gm[1] ~ normal(0,0.5); //global mean of threshold
  gm[2] ~ normal(2,2); //global mean of slope
  gm[3] ~ normal(-4,2); //global mean of confidence precision
  gm[4] ~ normal(3,2); //global mean of confidence precision
  gm[5] ~ normal(0,2); //global mean of meta uncertainty for correct trials 
  gm[6] ~ normal(0,2); //global mean of meta uncertainty for incorrect trials 
  gm[7] ~ normal(0,0.5); //global mean of meta uncertainty for incorrect trials 




  for (n in 1:N) {
    
    target += binomial_lpmf(a[n] | 1, theta[n]);   // likelihood for the outcomes

    target += ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]) + metabias, exp(confprec), c0, c11);   // likelihood for confidence on correct trials 


  }

    c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);   // likelihood for confidence cut-offs
    c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);

}

generated quantities {

  real c1 = c0 + exp(c11);    // actual high confidence cut-off

  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik_conf = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);


  for(n in 1:N){
    log_lik_bin[n] = binomial_lpmf(a[n] | 1, theta[n]);
    
    log_lik_conf[n] = ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]) + metabias, exp(confprec), c0, c11);
    
    // log_lik_conf[n] = ord_beta_reg_lpdf(C[n] | logit(theta_conf[n]) + meta_bias, exp(conf_prec), c0, c11);
    log_lik[n] = log_lik_conf[n] + log_lik_bin[n];
  }



}
