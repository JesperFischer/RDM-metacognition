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


data{
  
 int<lower=1> N;

  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=0, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  vector[N] RT;                         // stimulus strength
  vector[N] C;                         // stimulus strength

  real minRT;
  vector[N] XD;                         // stimulus strength
  
  
}

parameters{
  real a_int;
  real a_beta;
  real v_int;
  real v_beta;

  real<lower=0, upper =1> b;
  real<lower=0, upper = minRT> ndt;
  
  vector[N] evidence_ind;

  
  real sigma_m_log;              // decision bias
  real prec_conf_log;              // decision bias
  real c11;
  real c0;
  
}

transformed parameters{
  vector[N] a_i;
  vector[N] v_i;
  vector[N] ehat;
  vector[N] mu_conf;
  
  a_i = exp(a_int + a_beta * RT);
  v_i = v_int + v_beta * XD;
  
  for(i in 1:N){
    ehat[i] = a_i[i] + exp(sigma_m_log) * evidence_ind[i];

    if(a[i] == 1){
      mu_conf[i] = inv_logit(2*ehat[i]*X[i] / (1^2 + exp(sigma_m_log)^2));
    }else if(a[i] == 0){
      mu_conf[i] = 1- inv_logit(2*ehat[i]*X[i] / (1^2 + exp(sigma_m_log)^2));
    }
    
    
  }
  
  
}

model{
  
  sigma_m_log ~ normal(-1, 1);
  
  a_int ~ normal(5,2);
  a_beta ~ normal(-1,2);
  
  v_int ~ normal(0,3);
  v_beta ~ normal(2,4);
  
  b ~ normal(0,0.5);
  ndt ~ normal(0.4,0.3);
  
  
  evidence_ind ~ std_normal();

  
  for(i in 1:N){
    
    target += ord_beta_reg_lpdf(C[i] | logit(mu_conf[i]), exp(prec_conf_log), c0, c11);
    if(a[i] == 1){
      target += wiener_lpdf(RT[i] | a_i[i], ndt, b,v_i[i]);
    }else if(a[i] == 0){
     target += wiener_lpdf(RT[i] | a_i[i], ndt, 1-b,-v_i[i]);
  }

  }
  
}

