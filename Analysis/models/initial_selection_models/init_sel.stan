data {
  int<lower=0> N;
  int<lower=0> L; // number of observers
  int<lower=0,upper=L> Z[N]; // observer levels
  real<lower=0,upper=1> x[N];
  real<lower=0,upper=1> y[N];
}

parameters {
  real b[4];
  vector<lower=0>[4] sig_mu; //

  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[4] L_u;
  matrix[4, L] z_u;  // random effect matrix
}

transformed parameters {

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[4,L] u;
  u = diag_pre_multiply(sig_mu, L_u) * z_u;
}

model {

  // priors for mean (x, y)
  b[1] ~ normal(0.5, 0.1);
  b[2] ~ normal(0.5, 0.1);
  sig_mu ~ normal(0, 0.1);
  // priors for sd dev of (x, y)  
  b[3] ~ normal(-3, 0.1);
  b[4] ~ normal(-3, 0.1);
  // random effect structure
  L_u ~ lkj_corr_cholesky(1.5); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 1);  
  
  for (n in 1:N) {
    
    int z = Z[n];

    x[n] ~ normal(b[1] + u[1, z], exp(b[3] + u[3, z])) T[0, 1];
    y[n] ~ normal(b[2] + u[2, z], exp(b[4] + u[4, z])) T[0, 1];
  }

 }
