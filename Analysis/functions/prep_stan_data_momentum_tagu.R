prep_stan_data <- function(d)
{
  
  # function to take a tibble/dataframe of "spatial foraging" data, and prepare
  # it for running through the stan model
  
  # d should be a tibble with the following columns:
  # - trial: 1-n_trials
  # - id: target id 1:n_targets. If this does not exist, it will be created from "found"
  # - targ_type: factor level of each target type. For now, = 0 or 1.  
  # - x: hori position of target. Usually normalized to be [0, 1]
  # - y: similar to x
  # - found: the order in which the target was found. Set missed targets to -1
  # - condition: 1, 2, 3, .etc
  
  if (!"condition" %in% names(d)) d$condition = 1
  if (!"id" %in% names(d)) d$id = d$found
  
  d %>% mutate(trial = paste(observer, condition, trial, sep = "-"),
                 trial = as_factor(trial),
                 trial = as.numeric(trial),
               targ_type = if_else(targ_type == 0, -1, targ_type)) -> d
  
  
  # arrange by trial and id before calc distance matrix
  d <- arrange(d, observer, condition, trial, id) 
  
  # for each trial, compute target-target distance matrix
  # delta is indexed by ID
  deltas <- tibble()
  for (obs in unique(d$observer))  {
    print(obs)
    for (cc in unique(d$condition))   {
      for (trl in unique(d$trial))      {
        d_trial <- filter(d, trial == trl, observer == obs, condition == cc) %>% select(x, y)
        delta <- as.matrix(dist(d_trial, upper = TRUE, diag = TRUE))
        delta <- as_tibble(delta) 
        names(delta) = paste("dt", 1:nrow(d_trial), sep = "_")
        deltas <- bind_rows(deltas, delta)
      }
    }
  }
  d <- bind_cols(d, deltas)
 
  # put in first row
 
  theta <- tibble()
 
  for (ii in 1:nrow(d)) {
    
    print(paste(ii, nrow(d)))
   
     d_trial <- filter(d, trial == d$trial[ii], 
                       observer == d$observer[ii], 
                       condition == d$condition[ii]) 
   
     if (d$found[ii] %in% c(-1, 1)) {
      phi <- rep(0, nrow(d_trial))
    } else
    {
      d_targ <- d[ii,]
      d_prev_targ <- filter(d_trial, found == d$found[ii] - 1)
      
      phi = atan2((d_targ$y -  d_prev_targ$y), (d_targ$x -  d_prev_targ$x)) * 180/pi
      phi = (atan2((d_trial$y - d_targ$y), (d_trial$x - d_targ$x)) * 180/pi) - phi 
      phi = pmin(abs((phi %% 360)), abs((-phi %% 360)))
      phi = phi/180
      phi[ii] = 1
    }
    #   
   
      
      names(phi) = paste("at", 1:nrow(d_trial), sep = "_")
    #   # order by ID rather than found
       theta <- bind_rows(theta, phi[d_trial$id])
    # }
  }
    

  d <- bind_cols(d, theta)
  
  
  # create a n_trial x n_targets array that stores targ_type 
  d %>% select(observer, condition, trial, id, targ_type) %>%
    arrange(observer, condition, trial, id) %>%
    pivot_wider(c(observer, trial, condition), names_from = id, values_from = targ_type) %>%
    select(-trial, -observer, -condition) -> features
  
  # remove rows for targets that were never found
  d <- filter(d, found != -1) %>% 
    arrange(observer, condition, trial, found)
  

  
  stan_data <- list(
    N = nrow(d),
    L = length(unique(d$observer)),
    n_trials = nrow(features),
    K = as.integer(max(d$condition)),
    nT = n_targets,
    trialID = as.integer(d$trial),
    Y = as.integer(d$id),
    D = select(d, starts_with("dt")),
    E = select(d, starts_with("at")),
    W = as.matrix(features),
    selection_order = d$found,
    # S = feature_match,
    X = as.integer(d$condition),
    Z = as.integer(d$observer),
    prior_sig_b = prior_b,
    prior_sig_z = 1
  )
  
  
  return(stan_data)
   }
