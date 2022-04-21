data {
  int<lower=0> N;
  int<lower=0> L; // number of observers
  int<lower=0,upper=L> Z[N]; // observer levels
  real<lower=0,upper=1> x[N];
  real<lower=0,upper=1> y[N];
}

parameters {
  real<lower=0> mu_x[2];
  real<lower=0> mu_y[2];
  
  real sd_x[2];
  real sd_y[2];

  real<lower=0, upper=1> lambda[L];
}

transformed parameters {

}

model {

  // priors for mean (x, y)
  mu_x ~ normal(0.5, 0.1);
  mu_y ~ normal(0.0, 0.1);

  sd_x ~ normal(0, 1);
  sd_y ~ normal(0, 1);


  for (n in 1:N) {
    
    int z = Z[n];

   target += log_mix(lambda[z],
                    beta_lpdf(x[n] | mu_x[1], sd_x[1]),
                    beta_lpdf(x[n] | mu_x[2], sd_x[2])); 

  }

 }
