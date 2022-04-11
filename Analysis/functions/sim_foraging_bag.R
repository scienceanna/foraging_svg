library(tidyverse)

############################################
# Functions for data generation
############################################
simulate_observer <- function(ob_id = 1, observer_sig = 0.1) {
  
  # generate random weights
  obs_w_pref <- boot::inv.logit(boot::logit(prob_a) + rnorm(1, 0, observer_sig))
  obs_w_stay <- boot::inv.logit(boot::logit(prob_s) + rnorm(1, 0, observer_sig))
  
  
  d_ob <- map_df(1:n_conditions, 
                 simulate_data, 
                 n_trials = n_trials, n_targ = n_targ, prob_a =  obs_w_pref, prob_s =  obs_w_stay ) %>%
    mutate(
      observer = ob_id,
      condition = as.factor(condition), 
      obs_a = obs_w_pref[condition],
      obs_s = obs_w_stay[condition])
  
  return(d_ob)
  
}


simulate_data <- function(cc = 1, n_trials = 100, n_targ = 50, n_found =  NA, prob_a = 0.5, prob_s = 0.5, p_quit = 0.00) {
  
  # generate n_trials worth of data, all with the same parameters
  # include some basic run length statistics in the output
  
  if (is.na(n_found)) {
    n_found = n_targ
  }
  
  d <- map_df(1:n_trials, gen_foraging_seq, 
              ps = prob_s[cc], pa = prob_a[cc], n_tot = n_targ, n_found = n_found, p_q = p_quit) %>% 
    pivot_wider(names_from = "pick", values_from = "found") 
  
  ## add run-length and number-of-run stats
  d <- cbind(d, map_df(1:n_trials, get_run_info, d))
  
  d$condition = as.integer(cc)
  d$prob_a = prob_a[cc]
  d$prob_s = prob_s[cc]
  d$nA = as.integer(ceiling(n_targ/2))
  d$nB = as.integer(n_targ - d$nA)
  
  d <- select(d, condition, trial, nA, nB, prob_a, prob_s, n_found, max_run_length, all_run_lengths, n_runs, starts_with("t"))
  
  return(d)
}

gen_foraging_seq <- function(t = 1, ps = 0.5, pa = 0.5, p_q = 0, n_tot = 40, n_found = NA) {
  
  if (is.na(n_found)) {
    n_found = n_tot
  }
  
  # w_stay = weight on stick or twist
  # w_pref = bias for target type 1 over 0. 
  
  if (n_tot %% 2) {
    print("ERROR: n_tot should be even.")
  }
  
  nb <- c(n_tot/2, n_tot/2) # number of balls in bag
    
  found <- integer()
    
  # sample first ball
  p_1 = pa*nb[1] / (pa*nb[1] + (1-pa)*nb[2]) 
  pull = rbinom(1, 1, p_1)
  
  found <- append(found, pull)
  
  # remove ball from bag
  nb[2-pull] <- nb[2-pull]-1
  
  # has quit already
  given_up <- FALSE
  
  for (ii in 2:n_found) {
    # allow for early quitting
    if (given_up == FALSE & runif(1, 0, 1) < p_q) {
      given_up = TRUE
    }
    

    
    # change probability pA if we've ran out of a ball type
    if (nb[1] == 0) {pa = 0}
    if (nb[2] == 0) {pa = 1}
     
    # change probability of pS if we have ran out of ball type
    if (nb[2-found[ii-1]] == 0) {ps = 0}
    if (nb[2-(1-found[ii-1])] == 0) {ps = 1}
    
    if (given_up) {
      found <- append(found, -1)
    } else {
  
      # sample the ii-th ball
      p_1 = if_else(found[ii-1]==1, 
                    pa*ps*nb[1] / (pa*ps*nb[1] + (1-pa)*(1-ps)*nb[2]),
                    pa*(1-ps)*nb[1] / (pa*(1-ps)*nb[1] + (1-pa)*ps*nb[2]))
      
      pull = rbinom(1, 1, p_1) 
      found <- append(found, pull)
      # remove ball from bag
      nb[2-pull] <- nb[2-pull]-1
    }
  }
  
  return(tibble(trial = t, pick = paste("t", 1:n_found, sep = ""), found = found))
}

get_run_info <- function(t, d) {
  
  # calculate some simple run statistics for a trial
  
  trl <- filter(d, trial == t) 
  # remove trial index 
  trl <- trl[-1]
  # remove missing data values
  trl <- trl[trl != -1]
      
  rl <- rle(trl)
  return(list(
    n_found = length(trl),
    max_run_length = max(rl$lengths), 
    all_run_lengths =  list(rl$lengths), 
    n_runs = length(rl$lengths)))
}

get_run_info_multilevel <- function(d) {
  
  n_rows <- nrow(d) 
  
  n_runs <- integer(n_rows)
  max_run_length <- integer(n_rows)
  #all_run_lengths <- integer(n_rows)
  
  for (i in 1:n_rows){
    
    trl <- d %>%
      slice(i) %>%
      select(-trial, -condition, -observer, -n_found)
    
    rl <- rle(trl)
    
    n_runs[i] = length(rl$lengths)
    max_run_length[i] = max(rl$lengths)
    
  }
  
  d$n_runs <- n_runs
  d$max_run_length <- max_run_length
  
  return(d)
  
}

