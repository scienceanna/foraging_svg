
library(tidyverse)


d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"),
         targ_type = as.numeric(targ_type))


calc_trl_path <- function(obs, trl, cond) {
  dtrl <- filter(d, observer == obs, trial == trl, condition == cond) %>%
    mutate(x0 = lag(x),
           y0 = lag(y),
           d = sqrt((x - x0)^2 + (y - y0)^2)) %>%
    select(observer, condition, trial,x, x0, d)
  
  path_length = sum(dtrl$d, na.rm = TRUE)
  return(path_length)
}

trial_list <- filter(d, found == 1) %>%
  select(obs = "observer", trl = "trial", cond = "condition")

paths <- pmap_df(trial_list, calc_trl_path)
