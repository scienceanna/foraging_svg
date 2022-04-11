data {
  int <lower = 1> N; // Sample size - number of rows (trials x n_targets_found x observers)
  int <lower = 1> n_trials;   
  int <lower = 1> K; // number of conditions
  int <lower = 1> nT; // max number of targets per trial
  int <lower = 1, upper = N> trialID[N]; // trial ID for each row
  int <lower = 1> Y[N]; // which target is which
  int <lower = -1, upper = nT> selection_order[N]; // selection order
  matrix<lower = 0>[nT*n_trials,nT] D; // distance measures
  matrix<lower = 0>[N, nT] E; // direction measures
  int <lower = 1, upper = K> X[N]; // trial features
  matrix<lower = -1, upper = 1>[n_trials, nT] W; // target features per trial
  // matrix<lower = -1, upper = 1>[N, nT] S; // stick/switch (does this targ match prev targ) 
  real <lower = 0> prior_sig_b;

}

parameters {

  // b contains our fixed effects
  real b[4*K]; // we have two types of bias per condition [bA[1], bA[2], ..., bS[1], bS[2], ...]  
  
}

transformed parameters {

}

model {

  row_vector[nT] remaining_items;
  row_vector[nT] weights;
  int trl = 0; // counter for trial number
  int cc; // what condition are we currently in?

  row_vector[nT] m; // does this target match the previous target?

  real bA;
  real bS;
  real bP;
  real bD;

  //////////////////////////////////////////////////
  // Define Priors
  //////////////////////////////////////////////////

  // priors for fixed effects
  for (kk in 1:K) {
      b[1+4*(kk-1)] ~ normal(0, prior_sig_b);
      b[2+4*(kk-1)] ~ normal(0, prior_sig_b);
      b[3+4*(kk-1)] ~ normal(20, 10);
      b[4+4*(kk-1)] ~ normal(0, 1);
  } 


  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {  

    // if selection_order[ii], we are at the start of a new trial!
    if (selection_order[ii] == 1) {

      trl = trl + 1; // update trial counter     
      cc = X[ii]; // get conditions of current target/trial 

      // reset the remaining_items tracker
      remaining_items = rep_row_vector(1, nT);
      
      // get biases based on curent condition
      bA = b[1+4*(cc-1)]; 
      bS = b[2+4*(cc-1)];
      bP = b[3+4*(cc-1)];
      bD = b[4+4*(cc-1)];

      weights = inv_logit(bA * W[trl]);

    } 
    else 
    {
      // get biases based on curent condition

      bA = b[1+4*(cc-1)];
      bS = b[2+4*(cc-1)];
      bP = b[3+4*(cc-1)];
      bD = b[4+4*(cc-1)];

     // print(weights .* remaining_items);

      // check what matches previous target
      m = W[trl] * W[trl, Y[ii-1]];
      if (selection_order[ii] == 2) {
         weights =   inv_logit(bA*W[trl] + bS*m) .* exp(-bP*D[ii-1, 1:nT]); //

        } else {
          weights = inv_logit(bA*W[trl] + bS*m) .* exp(-bP*D[ii-1, 1:nT]) .* exp(-bD*E[ii-1, 1:nT]);
        }

      // set weights for found targets to 0 
      weights = weights .* remaining_items;
    
    }

    // normalise
    weights = weights/sum(weights);

    // likelihood! 
    Y[ii] ~ categorical(to_vector(weights));

    // now remove found target from list of remaining remaining_items
    remaining_items[Y[ii]] = 0;
  }
}

generated quantities {
  //real b_prior = normal_rng(0, prior_spread);
  //real p_prior = inv_logit(b_prior);eal 
  real prox_prior =  lognormal_rng(4, 1);//real prox_prior = lognormal_rng(0, prior_prox_sigma);
}
