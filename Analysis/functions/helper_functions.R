prepare_stan_input <- function(d, stan_model_name, N = 1, nT, nF = NA,
                               L = NA, prior_b = 1, prox_prior = NA) 
{
  
  # define some default params
  if (is.na(nF)) nF = nT

  
  # now setup list of inputs for our stan model
  
  if (stan_model_name == "bag_foraging_one_trial") 
  {
    Y = as.numeric(select(d, starts_with("t"), -"trial"))
    
    stan_data <- list(
      nT = nT,
      nA = d$nA,
      nB = d$nB,
      nF = nF,
      Y = Y,
      prior_sig_bS = prior_b)
  }
  
  if (stan_model_name == "bag_foraging")
  {
    Y = as.matrix(select(d, starts_with("t"), -"trial"))
    
    stan_data <- list(
      N = N,
      K = length(unique(d$condition)),
      nT = nT,
      nA = d$nA,
      nB = d$nB,
      nF = d$n_found,
      Y = Y,
      X = as.numeric(d$condition),
      prior_sig_bS = prior_b,
      prior_sig_bA = prior_b)
  }
  
  if (stan_model_name == "bag_foraging_mutli_level") 
  {
    Y = as.matrix(select(d, starts_with("t"), -"trial"))
    
    stan_data <- list(
      N = N * n_conditions * L,
      L = L,
      K = length(unique(d$condition)),
      nT = nT,
      nA = d$nA,
      nB = d$nB,
      Y = Y,
      X = as.numeric(d$condition),
      Z = d$observer,
      n_found = d$n_found,
      prior_sig_b = prior_b,
      prior_sig_z = 1)
  }
  
  if (stan_model_name == "spatial_foraging")
  {
    stan_data <- prepare_prox_data(d, nT, prior_b, prox_prior)
  }
  
  if (stan_model_name == "spatial_foraging_multi_level")
  {
    stan_data <- prepare_prox_data(d, nT, prior_b, prox_prior, ml = TRUE)
  }
  
  
  return(stan_data)
  
}



prepare_prox_data <- function(d, nT, prior_b, prox_prior, ml = FALSE)
{
  
  # function to take a tibble/dataframe of "spatial foraging" data, and prepare
  # it for running through the stan model
  
  # d should be a tibble with the following columns:
  # - trial: 1-n_trials
  # - id: target id 1:n_targets. If this does not exist, it will be created from "found"
  # - targ_type: factor level of each target type. For now, = 0 or 1.  
  # - x: hori position of target. Usually normalised to be [0, 1]
  # - y: similar to x
  # - found: the order in which the target was found. Set missed targets to -1
  # - condition: 1, 2, 3, .etc
  
  if (!"condition" %in% names(d)) d$condition <- 1
  if (!"id" %in% names(d)) d$id = d$found
  if (!ml) d$observer <- 1
  
  # check that trial IDs are unique
  n_trials_by_condition <- nrow(d %>% group_by(condition, trial) %>% summarise(n(), .groups = "drop"))
  n_trials <- nrow(d %>% group_by(trial) %>% summarise(n(), .groups = "drop"))
  if (n_trials != n_trials_by_condition) {
    #recode so that trial IDs are unique
    d %>% mutate(trial = paste(condition, trial, sep = "-"),
                 trial = as_factor(trial),
                 trial = as.numeric(trial)) -> d
  }
  
    # arrange by trial and id before calc distance matrix
  d <- arrange(d, observer, trial, found)  %>%
    mutate(trl_id = paste(observer, condition, trial))
              
  
  # for each trial, compute target-target distance matrix
  deltas <- tibble()
  for (ii in unique(d$trl_id)) {
  
    
      d_trial <- filter(d, trl_id == ii) %>% select(x, y)
      delta <- as.matrix(dist(d_trial, upper = TRUE, diag = TRUE))
      delta <- as_tibble(delta) 
      delta <- delta[1:nrow(d_trial), ] # remove info for not found targets
      names(delta) = paste("dt", 1:nT, sep = "_")
      deltas <- bind_rows(deltas, delta)
  
  }
  print(dim(deltas))
  d <- bind_cols(d, deltas)
  rm(deltas, delta)
  
  # now reorder rows from first target found to last target found
  d <- arrange(d, observer, trial, condition, found)
  print(d)
  
  # create a n_trial x n_targets array that stores targ_type 

    d %>% select(observer, trial, found, targ_type, id) %>%
      arrange(observer, trial, found) %>%
      mutate(targ_type = as.integer(if_else(targ_type == 0, -1, 1))) %>%
      pivot_wider(c(observer, trial), names_from = id, values_from = targ_type) %>%
      select(-trial, -observer) -> features
 
  
  # remove rows for targets that were never found
  d <- filter(d, found != -1)
  
  # dprev <- 0
  # dprev[2:nrow(d)] <- d$targ_type[1:(nrow(d)-1)]
  # dprev <- ifelse(dprev == 0, -1, 1)
  # dprev[which(d$found == 1)] = NA 
  # 
  # feature_match = array(dim = c(nrow(d), nT))
  # for (ii in 1:nrow(d)) {
  #   feature_match[ii, ] = as.integer(if_else(
  #    dprev[ii] == features[d$trial[ii],], 1, -1))
  # }
  # feature_match[is.na(feature_match)] = 0

    stan_data <- list(
      N = nrow(d),
      L = length(unique(d$observer)),
      n_trials = as.integer(max(d$condition)*n_trials),
      K = as.integer(max(d$condition)),
      nT = as.integer(n_targets),
      trialID = as.integer(d$trial),
      Y = as.integer(d$id),
      selection_order = as.integer(d$found),
      D = select(d, starts_with("dt_")),
      W = as.matrix(features),
      X = as.integer(d$condition),
      Z = as.integer(d$observer),
      prior_sig_b = prior_b,
      prior_sig_z = 1
    )
 
  
  return(stan_data)
}

