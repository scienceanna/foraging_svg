/*
Considers pS and pA
*/

data {
  int <lower = 1> N; // sample size - total number of trials
  int <lower = 1> K; // number of conditions
  int <lower = 1> nT; // number of targets in scene.
  int <lower = 1> nA[N]; // the number of targets of type A
  int <lower = 1> nB[N]; // the number of targets of type B
  int nF[N]; // number of targets found
  int <lower = 0, upper = 1> Y[N, nT]; // Matrix of data - each row is a trial
  int X[N]; // condition info, x = 1, 2
  real <lower = 0> prior_sig_bA;
  real <lower = 0> prior_sig_bS;
}

parameters {
  // I should change all these to use link functions
  real bA[K]; // 
  real bS[K]; //
}

transformed parameters {
  real pA[K] = inv_logit(bA);
  real pS[K] = inv_logit(bS);
}

model {

  real a; // weight for picking ball a
  real b; // weight for picking ball b
  int cc; // what condition is this trial?
  int balls[2]; // 2D array for counting how many targets (balls) of each type are left

  // define priors for bA and bS  
  bA ~ normal(0, prior_sig_bA);
  bS ~ normal(0, prior_sig_bS);

  // step through trial-by-trial
  for (n in 1:N) {

    cc = X[n]; // get condition
    balls = {nA[n], nB[n]}; // init our ball counter

    // give likelihood for first ball
    // first, calculate pA
    a = pS[cc] * pA[cc] * nA[n]; 
    b = (1-pS[cc]) * (1-pA[cc]) * nB[n];

    Y[n, 1] ~ bernoulli(a / (a + b));

    // remove the first ball from our counter
    balls[2-Y[n,1]] = balls[2-Y[n,1]] - 1;


    // now iterate over remaning targets
    for (ii in 2:nF[n]) {

      // probabilities will depend on prev target
      if (Y[n, ii-1] == 1) {
        a = pS[cc] * pA[cc] * balls[1];
        b = (1-pS[cc]) * (1-pA[cc]) * balls[2];
        } else { 
          a = pA[cc] * (1-pS[cc]) * balls[1];
          b = pS[cc] * (1-pA[cc]) * balls[2];
        }

        Y[n, ii] ~ bernoulli( a / (a + b) );

        // update number of balls
        balls[2-Y[n,ii]] = balls[2-Y[n,ii]] - 1; 
      } 

    }
  }

  generated quantities {
    real bA_prior = normal_rng(0, prior_sig_bA);
    real bS_prior = normal_rng(0, prior_sig_bS);
    real pA_prior = inv_logit(prior_sig_bA);
    real pS_prior = inv_logit(prior_sig_bS);
  }
