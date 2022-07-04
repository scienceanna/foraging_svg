data {
  int <lower = 1> N; // Sample size - number of rows (trials x n_targets_found x observers)
  int <lower = 1> L;  // number of levels (i.e, observers, chickens, etc)
  int <lower = 1> n_trials;   
  int <lower = 1> K; // number of conditions
  int <lower = 1> nT; // max number of targets per trial
  int <lower = 1, upper = N> trialID[N]; // trial ID for each row
  int <lower = 1> Y[N]; // which target is which
  int <lower = -1, upper = nT> selection_order[N]; // selection order
  matrix<lower = 0>[N, nT] D; // distance measures
  matrix<lower = 0>[N, nT] E; // direction measures
  int <lower = 1, upper = K> X[N]; // trial features
  matrix<lower = -1, upper = 1>[N/nT, nT] W; // target features per trial
  int <lower = 1, upper = L> Z[N]; // random effect levels 
  real <lower = 0> prior_sig_b;
  real <lower = 0> prior_sig_z;
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

  
  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {  

    // if selection_order[ii], we are at the start of a new trial!
    if (selection_order[ii] == 1) {

      trl = trl + 1; // update trial counter     
      cc = X[ii]; // get conditions of current target/trial
      ll = Z[ii]; // get observer level for this trial

      wx = log_mix(lambda[z],
                    beta_lpdf(x[n] | a_x[1], b_x[1]),
                    beta_lpdf(x[n] | a_x[2], b_x[2])); 

      wy = log_mix(lambda[z],
                    beta_lpdf(y[n] | a_y[1], b_y[1]),
                    beta_lpdf(y[n] | a_y[2], b_y[2])); 



      weights = inv_logit(bA * W[trl]);

    } 
    
  }
}

