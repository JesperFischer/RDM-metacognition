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
  int<lower=1> S;
  array[N] int S_id; 
  
  array[N] int<lower=-1, upper=1> D;   // true stimulus
  array[N] int<lower=0, upper=1> a;   // observed choice
  vector[N] X;                         // stimulus strength
  vector[N] C;                         // stimulus strength
  
  
}

transformed data{
  int P2 = 3;
  int P_norm = 2;
}

parameters {
  vector[N] evidence_std;
  vector[N] evidence_ind;
  
  
  vector[P_norm] gm;
  vector<lower=0>[P_norm] tau_u;
  matrix[P_norm, S] z_expo;    // Participant deviation from the group means


  vector<lower=0>[P2] gm2;
  vector<lower=0>[P2] tau_u2;
  
  vector<lower=0>[S] sigma_choice;
  vector<lower=0>[S] sigma_m;
  vector<lower=0>[S] sigma_e;
  

  
  vector[S] c11;
  vector[S] c0;
}


transformed parameters{

  matrix[S, P_norm] param;

  for(p in 1:P_norm){
    param[,p]= to_vector(gm[p] + tau_u[p] * z_expo[p,]);
  }

  vector[S] mean_choice = param[,1];
  vector[S] prec_conf_log = param[,2];
  
  
  vector[N] evidence;
  vector[N] p_action;
  vector[N] ehat;
  vector[N] mu_conf;
  
  
  for(i in 1:N){
    evidence[i] = mean_choice[S_id[i]] + X[i] * D[i] + (sigma_e[S_id[i]]) * evidence_std[i];
  
    ehat[i] = evidence[i] + (sigma_m[S_id[i]]) * evidence_ind[i];

    p_action[i] = Phi((evidence[i])/(sigma_choice[S_id[i]]));
    
    if(a[i] == 1){
      mu_conf[i] = inv_logit(2*ehat[i]*X[i] / ((sigma_e[S_id[i]])^2 + (sigma_m[S_id[i]])^2));
    }else if(a[i] == 0){
      mu_conf[i] = 1- inv_logit(2*ehat[i]*X[i] / ((sigma_e[S_id[i]])^2 + (sigma_m[S_id[i]])^2));
    }
  }
}

model {

  gm[1] ~ normal(0,0.3);
  gm[2] ~ normal(3,2);
  
  // alpha parameters
  gm2[1] ~ gamma(10,3);
  gm2[2] ~ gamma(10,3);
  gm2[3] ~ gamma(10,3);


  tau_u[1] ~ normal(0,0.3);
  tau_u[2] ~ normal(0,2);
  
  
  tau_u2[1] ~ gamma(10,2);
  tau_u2[2] ~ gamma(10,2);
  tau_u2[3] ~ gamma(10,2);

  
  
  sigma_choice~ gamma(gm2[1],tau_u2[1]);
  sigma_m~ gamma(gm2[2],tau_u2[2]);
  sigma_e~ gamma(gm2[3],tau_u2[3]);
  
  
  
  to_vector(z_expo) ~ std_normal();
  


  evidence_std ~ std_normal();
  evidence_ind ~ std_normal();
  
  
  a ~ bernoulli(p_action);

  for(i in 1:N){
    target += ord_beta_reg_lpdf(C[i] | logit(mu_conf[i]), exp(prec_conf_log[S_id[i]]), c0[S_id[i]], c11[S_id[i]]);
  }
  
  for(s in 1:S){
      c0[s] ~ induced_dirichlet([1,10,1]', 0, 1, c0[s], c11[s]);
      c11[s] ~ induced_dirichlet([1,10,1]', 0, 2, c0[s], c11[s]);
  }


}

generated quantities{
  
}
