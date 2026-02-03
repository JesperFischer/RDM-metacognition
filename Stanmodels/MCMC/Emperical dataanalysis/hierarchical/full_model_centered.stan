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
  int<lower=0> N;
  int<lower=0> S;
  array[N] int S_id;

  vector[N] C;

  vector[N] X;
  vector[N] XD;

  vector[N] ACC; // Vector of deltaBPM values that match the binary response
  array[N] int a;                      // stimulus strength
  
  
}

transformed data{
  int P = 5;
}

parameters {
  
  vector[N] evidence;
  vector[N] ehat;

  vector[P] gm;
  vector<lower=0>[P] tau_u;
  cholesky_factor_corr[P] L_u;    // Between participant cholesky decomposition
  matrix[P, S] z_expo;    // Participant deviation from the group means


  vector[S] c0;
  vector[S] c11;
}


transformed parameters{
  
  // Extracting individual deviations for each subject for each parameter
  matrix[S, P] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';

  matrix[S, P] param;

  for(p in 1:P){
    param[,p]= gm[p] + indi_dif[,p];
  }

  vector[S] mean_choice = (param[,1]);
  vector[S] sigma_e_log = (param[,2]);


  vector[S]  conf_prec1 = param[,3];  // confidence precision
  vector[S]  sigma_m_log = param[,4];  // meta uncertainty on correct trials
  vector[S]  sigma_choice_log = param[,5];  // meta uncertainty on incorrect trials


  // vector[N] evidence;
  vector[N] p_action;
  // vector[N] ehat;
  vector[N] mu_conf;
  
  for(i in 1:N){
    // evidence[i] = mean_choice[S_id[i]] + XD[i] + exp(sigma_e_log[S_id[i]]) * evidence_std[i];

    // ehat[i] = evidence[i] + exp(sigma_m_log[S_id[i]]) * evidence_ind[i];
    
    p_action[i] = Phi((evidence[i])/exp(sigma_choice_log[S_id[i]]));
    
    if(a[i] == 1){
      mu_conf[i] = inv_logit(2*ehat[i]*X[i] / (exp(sigma_e_log[S_id[i]])^2 + exp(sigma_m_log[S_id[i]])^2));
    }else if(a[i] == 0){
      mu_conf[i] = 1- inv_logit(2*ehat[i]*X[i] / (exp(sigma_e_log[S_id[i]])^2 + exp(sigma_m_log[S_id[i]])^2));
    }
  }
}

model {
  
  gm[1] ~ normal(0,0.3); //global mean of threshold 
  gm[2] ~ normal(-1,1); //global mean of slope
  gm[3] ~ normal(3,2); //global mean of confidence precision
  gm[4] ~ normal(-1,1); //global mean of meta uncertainty
  gm[5] ~ normal(-2,1); //global mean of meta bias


  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(0 , 0.3);
  tau_u[2] ~ normal(0 , 2);
  tau_u[3:5] ~ normal(0 , 2);
  L_u ~ lkj_corr_cholesky(2);

  for(i in 1:N){
    evidence[i] ~ normal(mean_choice[S_id[i]] + XD[i], exp(sigma_e_log[S_id[i]]));
    ehat[i] ~ normal(mean_choice[S_id[i]] + XD[i], sqrt(exp(sigma_e_log[S_id[i]])^2 + exp(sigma_m_log[S_id[i]])^2));
  }
  
  
  a ~ bernoulli(p_action);
  
  for(i in 1:N){
    
    target += ord_beta_reg_lpdf(C[i] | logit(mu_conf[i]), exp(conf_prec1[S_id[i]]), c0[S_id[i]], c11[S_id[i]]);

    
  }
  
  for(s in 1:S){
    c0[s] ~ induced_dirichlet([1,10,1]', 0, 1, c0[s], c11[s]);
    c11[s] ~ induced_dirichlet([1,10,1]', 0, 2, c0[s], c11[s]);
  }

}
