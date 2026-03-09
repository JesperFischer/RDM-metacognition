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
  array[N] int<lower=0, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  vector[N] C;                         // stimulus strength
  vector[N] interval;
  
  
}


parameters {
  vector[N] evidence_std;
  vector[N] evidence_ind;
  real sigma_choice_log;              // decision bias
  real sigma_m_log;              // decision bias
  real sigma_m_log_beta;              // decision bias
  real prec_conf_log;              // decision bias
  real sigma_e_log;              // decision bias
  real mean_choice;
  real c11;
  real c0;
}


transformed parameters{
  vector[N] evidence;
  vector[N] p_action;
  vector[N] ehat;
  vector[N] mu_conf;
  
  for(i in 1:N){
    evidence[i] = mean_choice + X[i] * D[i] + exp(sigma_e_log) * evidence_std[i];
    ehat[i] = evidence[i] + exp(sigma_m_log + sigma_m_log_beta * interval[i]) * evidence_ind[i];
    
    p_action[i] = Phi((evidence[i])/exp(sigma_choice_log));
    if(a[i] == 1){
      mu_conf[i] = inv_logit(2*ehat[i]*X[i] / (exp(sigma_e_log)^2 + exp(sigma_m_log+ sigma_m_log_beta * interval[i])^2));
    }else if(a[i] == 0){
      mu_conf[i] = 1- inv_logit(2*ehat[i]*X[i] / (exp(sigma_e_log)^2 + exp(sigma_m_log+ sigma_m_log_beta * interval[i])^2));
    }
  }
}

model {
   sigma_choice_log ~ normal(-1, 1);
  sigma_e_log ~  normal(-1, 1);
  sigma_m_log ~  normal(-1, 1);
  prec_conf_log ~ normal(3, 3);
  mean_choice ~ normal(0, 0.3);
  sigma_m_log_beta ~ normal(0,0.1);
  
  evidence_std ~ std_normal();
  evidence_ind ~ std_normal();
  
  
  a ~ bernoulli(p_action);
  
  
  for(i in 1:N){
    
    target += ord_beta_reg_lpdf(C[i] | logit(mu_conf[i]), exp(prec_conf_log), c0, c11);

    
  }
  
  c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
  c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);
  

}

generated quantities{
  real c1 = c0 + exp(c11);

  real sigma_choice = exp(sigma_choice_log);
  real sigma_m = exp(sigma_m_log);
  real sigma_e = exp(sigma_e_log);
  real prec_conf = exp(prec_conf_log);


}
