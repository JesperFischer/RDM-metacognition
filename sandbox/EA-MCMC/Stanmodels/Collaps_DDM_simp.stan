
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
  

}

transformed parameters{
  vector[N] a_i;
  vector[N] v_i;
  
  
  a_i = exp(a_int + a_beta * RT);
  v_i = v_int + v_beta * XD;
  
}

model{
  

  a_int ~ normal(5,2);
  a_beta ~ normal(-1,2);
  
  v_int ~ normal(0,3);
  v_beta ~ normal(2,4);
  
  b ~ normal(0,0.5);
  ndt ~ normal(0.4,0.3);
  


  
  for(i in 1:N){
    
    if(a[i] == 1){
      target += wiener_lpdf(RT[i] | a_i[i], ndt, b,v_i[i]);
    }else if(a[i] == 0){
     target += wiener_lpdf(RT[i] | a_i[i], ndt, 1-b,-v_i[i]);
  }

  }
  
}

