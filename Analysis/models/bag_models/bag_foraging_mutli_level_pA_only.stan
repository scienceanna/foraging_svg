data {
  int <lower = 1> N;  // Sample size - total number of trials
  int <lower = 1> L;  // number of levels (i.e, observers, chickens, etc)
  int <lower = 1> K; // number of conditions
  int <lower = 1> nT; // max number of targets per trial
  int <lower = 1> nA[N]; // the number of targets of type A
  int <lower = 1> nB[N]; // the number of targets of type B
  int <lower = -1, upper = 1> Y[N, nT]; // matrix of target data to fit model to - each row is a trial
  int X[N]; // condition info, x = 1, 2
  int <lower = 1, upper = L> Z[N]; // observer ID info
  int <lower = 1, upper = nT> n_found[N]; // number of targets found on each trial
  real <lower = 0> prior_sig_b;
  real <lower = 0> prior_sig_z;// prior for our A and stick weights (i.e, ~ N(0, prior_sigma)

}

parameters {

  real b[K]; // we have two types of bias per condition [bA[1], bA[2], ..., bS[1], bS[2], ...]  
  vector<lower=0>[K] sig_b; // random effect sigma for biases

  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[K] L_u;
  matrix[K,L] z_u;  // random effect matrix
}

transformed parameters {

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[K,L] u;
  u = diag_pre_multiply(sig_b, L_u) * z_u;
}

model {

  real wa; // weight for picking ball a
  real wb; // weight for picking ball b

  real pA;
  real pS; 
  
  // we will now step through trial-by-trial, but first, make some variables for tracking things
  int ll; // for storing the current observer level
  int kk; // for storing the current condition level
  int balls[2]; // 2D array for counting how many targets (balls) of each type are left

  // priors
  b ~ normal(0, prior_sig_b);
  sig_b ~ normal(0, prior_sig_z);
  L_u ~ lkj_corr_cholesky(1.5); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 1);  
  
  for (n in 1:N) {

    // put balls in the bag!
    balls[1] = nA[N]; // the number of balls (targets) of type A = 1
    balls[2] = nB[N]; // the number of balls (targets) of type B = 0
   
    ll = Z[n]; // get observer level for this trial
    kk = X[n]; // get condition level for this trial

    /*
    compute likelihood for the first target...
    first, calculate weights for a and b
    */
    pA = inv_logit(b[kk] + u[kk, ll]);

    wa = pA * nA[N]; 
    wb =  (1-pA) * nB[N];

    Y[n, 1] ~ bernoulli(wa / (wa + wb));

    /* 
    remove the first ball from our counter
    if the current ball = 1, then 2-1 = 1
    if the current ball = 0, then 2-0 = 2
    */
    balls[2-Y[n,1]] = balls[2-Y[n,1]] - 1;

    // now iterate over remaining targets
    for (ii in 2:n_found[n]) {

      // probabilities will depend on prev target
      if (Y[n, ii-1] == 1)
      {
        wa = pA * balls[1];
        wb = (1-pA) * balls[2];
      } else  
      { 
        wa = pA  * balls[1];
        wb = (1-pA) * balls[2];
      }

      Y[n, ii] ~ bernoulli( wa / (wa + wb) );

      // update number of balls
      balls[2-Y[n,ii]] = balls[2-Y[n,ii]] - 1; 

    } // end counting over balls (targets) 
  } // end counting over trials
}

generated quantities {
  real b_prior = normal_rng(0, prior_sig_b);
}
