functions{
  
  
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
  int<lower=1> N;

  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=-1, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  array[N] real C;    // reported confidence
}


parameters {
  real alpha;              // decision bias
  real beta;
  real beta_m;
  // real beta_mean;          // average sensitivity (beta + beta_m)/2
  // real beta_diff;          // difference in sensitivity (beta - beta_m)
  real conf_prec;          // confidence noise
  real delta;              // confidence bias
  vector[N] evidence_std;
  vector[N] evidence_m_std;
  real c0;
  real c11;
}


transformed parameters{
  vector[N] evidence;
  vector[N] evidence_m;
  vector[N] confidence;
  real p_correct = inv_logit(delta);
  
  // Transform back to beta and beta_m
  // real beta = beta_mean + beta_diff / 2.0;
  // real beta_m = beta_mean - beta_diff / 2.0;
  // real beta_m = beta + beta_diff;
  
  for(i in 1:N){
    
    evidence[i] = exp(beta) * (D[i] * X[i] - alpha) + evidence_std[i];
    evidence_m[i] = exp(beta_m) * (D[i] * X[i] - alpha) + evidence_m_std[i];
    


    real loglike_m_correct = (normal_lpdf(evidence_m[i] | exp(beta_m) * (a[i] * X[i] - alpha), 1));

    real loglike_m_incorrect = (normal_lpdf(evidence_m[i] | exp(beta_m) * (-a[i] * X[i] - alpha), 1));
    
    real logp_cor =  loglike_m_correct + log(p_correct);
    
    real logp_incor = loglike_m_incorrect + log1m(p_correct);

    confidence[i] = inv_logit(logp_cor - logp_incor);
  }
}

model {
  alpha ~ normal(0, 1);
  beta  ~ normal(1, 1);
  beta_m  ~ normal(1, 1);
  // beta_mean ~ normal(1, 1);      // Prior on average sensitivity
  // beta_diff ~ normal(0, 1);    // Prior on difference (often small)
  delta  ~ normal(0, 0.5);
  conf_prec ~ normal(0,3);
  evidence_std ~ std_normal();
  evidence_m_std ~ std_normal();
  
  for (i in 1:N) {

    real p_choice = Phi(exp(beta) * (D[i] * X[i] - alpha));
    target += bernoulli_lpmf((a[i] + 1) / 2 | p_choice);

    target += ord_beta_reg_lpdf(C[i] | logit(confidence[i]), exp(conf_prec), c0, c11);

  }
  c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
  c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);


}

