comp_weights <- function(trl, params, cond, obs)
{
  
  d_trl <- filter(d, observer == obs, condition == cond, trial == trl)
  
  if (nrow(d_trl) == 0) {
    
    
  }
  else
  {
  
  
  d_trl$b = NaN
  d_trl$bMax = NaN
  d_trl$err = NaN
  d_trl$model_pref = NaN
  
  
  for (t in 2:40) {
    
    prev_targ <- d_trl$targ_type[t-1]
    
    
    x0 <- d_trl$x[t-1]
    y0 <- d_trl$y[t-1]
    
    x00 <- d_trl$x[t-2]
    y00 <- d_trl$y[t-2]
    
    d_remain <- filter(d_trl, found>=t)
    
    match_prev = if_else(d_remain$targ_type == prev_targ, 1, -1)
    
    if (t == 2) {
      d_remain %>% mutate(
        dist = (x0 - x)^2 + (y0 - y)^2,
        dist = sqrt(dist),
        prox = exp((-params$bP * dist)),
        b = params$pA * targ_type + params$pS * match_prev,
        b = boot::inv.logit(b) * prox,
        b = b/sum(b)) -> d_remain
    } else {
      d_remain %>% mutate(
        dist = (x0 - x)^2 + (y0 - y)^2,
        dist = sqrt(dist),
        phi = atan2(( y0 - y00), (x0 - y00)) * 180/pi,
        phi = (atan2((y - y0), (x - x0)) * 180/pi) - phi ,
        phi = pmin(abs((phi %% 360)), abs((-phi %% 360))),
        phi = phi/180,
        prox = exp((-params$bP * dist) - (params$bM * phi)),
        b = params$pA * targ_type + params$pS * match_prev,
        b = boot::inv.logit(b) * prox,
        b = b/sum(b)) -> d_remain
      
    }
    
    d_trl$b[t] <- d_remain$b[which(d_remain$found==t)]
    d_trl$bMax[t] <- d_trl$b[t] == max(d_remain$b)
    d_trl$err[t] <-  max(d_remain$b) - d_trl$b[t]
    d_trl$model_pref[t] <- d_remain$found[which(d_remain$b == max(d_remain$b) )]
    
    
  }}
  
  return(d_trl)
}
