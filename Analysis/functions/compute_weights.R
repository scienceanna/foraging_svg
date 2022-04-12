compute_weights <- function(trl, params, cond, obs)
{
  
  # function to step through empirical foraging data and compute the item weights
  # for every target selection
  
  # select the revelant trial
  d_trl <- filter(d, observer == obs, condition == cond, trial == trl)
  
  if (nrow(d_trl) == 0) {
    
    # if the data for this trial is missing, do nothing
    
  }
  else
  {
  
    # init column for the weight stats we want to select
    d_trl$b = NaN # weight of target next selected by participant
    d_trl$selected_max = NaN # did the participant select the most likely target?
    d_trl$max_b = NaN # weight of most likely target
    d_trl$model_pref = NaN # item ID for the most likely target
    
    # step through each target selection....
    for (t in 2:40) {
      
      # get info about the previously selected item
      prev_targ <- d_trl$targ_type[t-1]
      
      x0 <- d_trl$x[t-1]
      y0 <- d_trl$y[t-1]
      
      # and the item selected before that...
      x00 <- d_trl$x[t-2]
      y00 <- d_trl$y[t-2]
      
      # get list of items that haven't been selected yet
      d_remain <- filter(d_trl, found>=t)
      
      # check which items match the previously selected item
      match_prev = if_else(d_remain$targ_type == prev_targ, 1, -1)
      
      # compute proximity info
      d_remain %>% mutate(
        dist = (x0 - x)^2 + (y0 - y)^2,
        dist = sqrt(dist)) -> d_remain
      
      if (t == 2) {
        
        # if it's the second item selection, compute spatial weight based just on distance
        d_remain %>% mutate( 
          prox = exp(-params$bP * dist)) -> d_remain 
        
      } else {
        
        # but if t>2, we first compute phi and then compute spatial weight
        d_remain %>% mutate(
          phi = atan2(( y0 - y00), (x0 - y00)) * 180/pi,
          phi = (atan2((y - y0), (x - x0)) * 180/pi) - phi ,
          phi = pmin(abs((phi %% 360)), abs((-phi %% 360)))/180,
          prox = exp((-params$bP * dist) - (params$bM * phi))) -> d_remain
      }
      
      # now compute weights for each item
      d_remain %>% mutate(
        b = params$pA * targ_type + params$pS * match_prev,
        b = boot::inv.logit(b) * prox,
        b = b/sum(b)) -> d_remain
        
   
      # finally compute some stats about the item weights
      d_trl$b[t] <- d_remain$b[which(d_remain$found==t)]
      d_trl$selected_max[t] <- d_trl$b[t] == max(d_remain$b)
      d_trl$max_b[t] <-  max(d_remain$b) - d_trl$b[t]
      d_trl$model_pref[t] <- d_remain$found[which(d_remain$b == max(d_remain$b) )]
      
    }
  }
  
  return(d_trl)
}
