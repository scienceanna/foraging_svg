library(tidyverse)


sim_spatial_trials <- function(n = 10, m = 10, tune_prox = 1, b_pref = 0, b_switch = 0, n_trials = 10, spat_config = FALSE) 
{
  # n is number of targets present
  # m is number of targets to find
  
  params <- list(
    n = n,
    m = m,
    tune_prox = tune_prox,
    b_pref = b_pref,
    b_switch = b_switch,
    spat_config = spat_config
  )
  
  dat <- map_df(1:n_trials, sim_spatial_trial, params = params)
  
  return(dat)
}

sim_spatial_trial <- function(params, trial = 1)
{
  
  # n is total number of targets
  # m is the number that we want to find

  d_trial <- tibble(
    trial = trial,
    id = 1:params$n,
    x  = runif(params$n),
    y = runif(params$n),
    type = rbinom(params$n, 1, 0.5),
    found = -1,
    dist = NA)
  
  if (params$spat_config==TRUE) 
  {
    
    d_trial %>% mutate(
      type = if_else(x < 0.5, 0, 1)
    ) -> d_trial
    
  }
  
  #d_trial <- mutate(d_trial, type = if_else(type ==0, -1, 1))
  
  # compute proximity info
  delta <- as.matrix(dist(select(d_trial, x, y), upper = TRUE, diag = TRUE))
  delta <- as_tibble(delta) 
  names(delta) = paste("dt", 1:params$n, sep = "_")
  d_trial <- bind_cols(d_trial, delta)
  
  #### TODO
  # use the precompute prox info in the simulation below
  
  d_remain <- d_trial %>% 
    mutate(
      dist = 0,
      prox = 0,
      b = boot::inv.logit(params$b_pref * type),
      b = b/sum(b))
  
  # pick a first point at random
  d_found <- sample_n(d_remain, 1, weight = b) %>%
    mutate(t = 1)
  
  # remove this point from the stimuli
  d_remain <- filter(d_remain, id != d_found$id)
  d_trial$found[d_found$id[1]] <- 1
  
 # varname <- paste("dt",d_found$id[1], sep = "_" )
 # d_trial[which(d_trial$found == -1),][varname] = 999
  
  for (t in 2:params$m) {
    
    prev_targ <- d_found$type[t-1]
    match_prev = if_else(d_remain$type == prev_targ, 1, -1)
    # compute proximity of remaining targets from current target
    d_remain %>% mutate(
      dist = (d_found$x[t-1] - x)^2 + (d_found$y[t-1] - y)^2,
      dist = sqrt(dist),
      prox = exp(-params$tune_prox * dist),
      b = params$b_pref * type + params$b_switch * match_prev,
      b = boot::inv.logit(b) * prox,
      b = b/sum(b)) -> d_remain
    
    # sample the next target
    
    d_found %>% add_row(
      sample_n(d_remain, 1, weight = b) %>% mutate(t = t)) -> d_found
    
    d_remain <- filter(d_remain, id != d_found$id[t]) 
    
    d_trial$found[d_found$id[t]] <- t
    d_trial$dist[d_found$id[t]] <- as.matrix(select(d_trial, starts_with("dt")))[d_found$id[t],d_found$id[t-1]]
    
  } 
  
  
  return(d_trial)
}

