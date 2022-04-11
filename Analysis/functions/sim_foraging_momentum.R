library(tidyverse)


sim_momentum_post <- function(df_params, n_trials = 10, n = 40)
{
  
  d <- tibble()
  
  for (ii in 1:nrow(df_params)) {
    
    print(ii)
    
    sim_momentum_trials(n = n, m = n, 
                        tune_prox = df_params$bP[ii], 
                        tune_angle = df_params$bM[ii], 
                        b_pref = df_params$bA[ii],
                        b_switch = df_params$bS[ii], 
                        n_trials = 10) %>%
      mutate(observer = df_params$observer[ii],
             condition = df_params$condition[ii],
             draw = df_params$.draw[ii]) %>%
      arrange(trial, found) %>%
      bind_rows(d) -> d
    
  }
  
  return(d)
  
}
  
get_runs <- function(d) {
  
  dout <- tibble()
  
  for (obs in unique(d$observer)) {
    for (cd in unique(d$condition)) {
      for (trl in unique(d$trial)) {
        
        targs <- filter(d, observer == obs, condition == cd, trial == trl)$targ_type
        
        if (length(targs)>0) {
          rl <- rle(targs)
          
          dout <- bind_rows(dout, 
                            tibble(observer = obs,
                                   trial = trl,
                                   condition = cd, 
                                   max_run = max(rl$lengths),
                                   num_run = length(rl$lengths),
                                   men_run = mean(rl$lengths)))
        }
      }
    }
  }
  
  return(dout)
  
}


sim_momentum_trials <- function(n = 10, m = 10, tune_prox = 1, tune_angle = 5, b_pref = 0, b_switch = 0, n_trials = 10, spat_config = FALSE) 
{
  # n is number of targets present
  # m is number of targets to find
  
  params <- list(
    n = n,
    m = m,
    tune_prox = tune_prox,
    tune_angle = tune_angle,
    b_pref = b_pref,
    b_switch = b_switch,
    spat_config = spat_config
  )
  
  #dat <- sim_spatial_trial(params)
  dat <- map_df(1:n_trials, sim_spatial_trial, params = params) %>%
    rename(targ_type = "type") %>% 
    mutate(observer = 1, condition = 1) %>%
    select(observer, condition, trial, id, targ_type, found, x, y, phi, dist)
  
  return(dat)
}

sim_spatial_trial <- function(params, trial = 1)
{
  
  # n is total number of targets
  # m is the number that we want to find

  d_trial <- 
    tibble(
      trial = trial,
          id = 1:params$n,
          x  = runif(params$n),
          y = runif(params$n),
          type = sample(rep(c(0, 1), params$n/2)),
          found = -1,
          dist = NA,
          phi = NA)
  
  d_trial$type[which(d_trial$type == 0)] <- -1
  
  if (params$spat_config==TRUE) 
  {
    
    d_trial %>% mutate(
      type = if_else(x < 0.5, -1, 1)
    ) -> d_trial
    
  }

  
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
  
 # second target select.. add in distance...
  t = 2
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
  d_trial$dist[d_found$id[t]] <- d_found$prox[t]
  
  for (t in 3:params$m) {
    
    prev_targ <- d_found$type[t-1]
    match_prev = if_else(d_remain$type == prev_targ, 1, -1)
    # compute proximity of remaining targets from current target
    d_remain %>% mutate(
      dist = (d_found$x[t-1] - x)^2 + (d_found$y[t-1] - y)^2,
      dist = sqrt(dist),
      phi = atan2(( d_found$y[t-1] - d_found$y[t-2]), (d_found$x[t-1] - d_found$x[t-2])) * 180/pi,
      phi = (atan2((y - d_found$y[t-1]), (x - d_found$x[t-1])) * 180/pi) - phi ,
      phi = pmin(abs((phi %% 360)), abs((-phi %% 360))),
      phi = phi/180,
      prox = exp((-params$tune_prox * dist) - (params$tune_angle * phi)),
      b = params$b_pref * type + params$b_switch * match_prev,
      b = boot::inv.logit(b) * prox,
      b = b/sum(b)) -> d_remain
    
    # sample the next target
    
    d_found %>% add_row(
      sample_n(d_remain, 1, weight = b) %>% mutate(t = t)) -> d_found
    
    d_remain <- filter(d_remain, id != d_found$id[t]) 
    
    d_trial$found[d_found$id[t]] <- t
    d_trial$dist[d_found$id[t]] <- d_found$prox[t]
    d_trial$phi[d_found$id[t]] <- d_found$phi[t]
    
  } 
  
  
  return(d_trial)
}


get_run_info_momentum <- function(d) {
  
  tot_trials <- max(d$id)
  
  n_rows <- nrow(d)/tot_trials 
  
  n_runs <- integer(n_rows)
  max_run_length <- integer(n_rows)
  run_info <- tibble(
    trial = 1:n_rows
  )
  
  for (i in 1:n_rows){
    
    trl <- d %>%
      slice((((i-1)*tot_trials)+1):(tot_trials*i)) %>%
      select(observer, condition, trial, targ_type, found) %>%
      arrange(observer, condition, trial, found)
    
    rl <- rle(trl$targ_type)
    
    n_runs[i] = length(rl$lengths)
    max_run_length[i] = max(rl$lengths)
    
  }
  
  run_info$n_runs <- n_runs
  run_info$max_run_length <- max_run_length
  
  return(run_info)
  
}
