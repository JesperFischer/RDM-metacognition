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

// Stan model for binary choice + confidence
// Evidence marginalized out, psychometric function for choice
// Beta distribution for confidence based on P(correct)

data {
  int<lower=0> N;                    // Number of trials
  array[N] int<lower=-1, upper=1> D; // True stimulus state (-1 or 1)
  vector[N] X;                       // Stimulus strength (unsigned)
  array[N] int<lower=-1, upper=1> a; // Decision/action (-1 or 1)
  vector<lower=0, upper=1>[N] C;   // Confidence (bounded [0.5, 1])
}


parameters {
  real alpha;                        // Decision bias/threshold
  real beta;                // Decision sensitivity
  real delta; // Confidence bias
  real conf_prec;                // Decision sensitivity
  real c0;
  real c11;
  
}

model {
  // Priors
  alpha ~ normal(0, 1);
  beta ~ normal(2, 1);              // Weakly informative, mean = 4
  delta ~ normal(0, 1);            // Confidence bias
  conf_prec ~ normal(2, 1);              // Weakly informative, mean = 4
  
  // Likelihood
  for (i in 1:N) {
    // Choice: Psychometric function using Bernoulli distribution
    real prob_action_1 = Phi(exp(beta) * (D[i] * X[i] - alpha));
    target += bernoulli_lpmf((a[i] + 1) / 2 | prob_action_1);
    
    // Confidence: Based on P(correct)
    // P(correct) depends on the true stimulus direction D[i]
    real p_correct;
    if (D[i] == 1) {
      p_correct = prob_action_1;  // P(action=1) when stimulus is +1
    } else {
      p_correct = 1 - prob_action_1;  // P(action=-1) when stimulus is -1
    }


    target += ord_beta_reg_lpdf(C[i] | logit(p_correct) + delta, exp(conf_prec), c0, c11);

  }
  c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
  c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);
    
}

generated quantities{
  real beta1 = exp(beta);
}