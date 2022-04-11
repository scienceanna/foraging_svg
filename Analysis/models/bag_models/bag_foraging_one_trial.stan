/*
This is the simple versin of our model.
Only considers the stick/switch bias pS
One trial, one codition, one participant. 
*/

data {  
  int <lower = 1> nT; // number of targets in scene.
  int <lower = 1> nA; // the number of targets of type A
  int <lower = 1> nB; // the number of targets of type B
  int <lower = 1, upper = nT> nF; // number of targets found 
  int <lower = -1, upper = 1> Y[nF]; // Matrix of data - each row is a trial
  real <lower = 0> prior_sig_bS;
}

parameters {
  real bS; 
}

transformed parameters {
  real pS = inv_logit(bS);
}

model {
  real a; // weight for picking ball a
  real b; // weight for picking ball b

  int balls[2]; // counting how many targets are left
  balls = {nA, nB};

  // define priors for bS
  bS ~ normal(0, prior_sig_bS);

  // give likelihood for first ball
  a = nA; 
  b = nB;
  Y[1] ~ bernoulli(a / (a + b));
  // remove the first ball from our counter
  balls[2-Y[1]] = balls[2-Y[1]] - 1;

  // now iterate over remaning targets
  for (ii in 2:nF) {
    // probabilities will depend on prev target
    if (Y[ii-1] == 1) {
      a = pS * balls[1];
      b = (1-pS) * balls[2];
    } else { 
      a = (1-pS) * balls[1];
      b = pS * balls[2];
    }

    Y[ii] ~ bernoulli( a / (a + b) );
    // update number of balls
    balls[2-Y[ii]] = balls[2-Y[ii]] - 1; 
  }
}

generated quantities {
  real b_prior = normal_rng(0, prior_sig_bS);
  real p_prior = inv_logit(b_prior);
}
