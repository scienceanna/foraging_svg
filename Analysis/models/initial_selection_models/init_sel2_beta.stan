data {
  int<lower=0> N;
  int<lower=0> L; // number of observers
  int<lower=0,upper=L> Z[N]; // observer levels
  real<lower=0,upper=1> x[N];
  real<lower=0,upper=1> y[N];
}

parameters {
  real<lower=0> a_x[2];
  real<lower=0> a_y[2];
  
  real<lower=0> b_x[2];
  real<lower=0> b_y[2];

  real<lower=0, upper=1> lambda[L];
}

transformed parameters {

}

model {

  // priors for mean (x, y)
  a_x ~ normal(1, 0.25);
  a_y ~ normal(1, 0.25);

  b_x ~ normal(1, 1);
  b_y ~ normal(1, 1);


  for (n in 1:N) {
    
    int z = Z[n];
    
   target += log_mix(lambda[z],
                    beta_lpdf(x[n] | a_x[1], b_x[1]),
                    beta_lpdf(x[n] | a_x[2], b_x[2])); 

   target += log_mix(lambda[z],
                    beta_lpdf(y[n] | a_y[1], b_y[1]),
                    beta_lpdf(y[n] | a_y[2], b_y[2])); 

  }

 }
