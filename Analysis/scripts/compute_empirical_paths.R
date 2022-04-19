library(tidyverse)
source("../functions/sim_foraging_momentum.R")


d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"),
         targ_type = as.numeric(targ_type))

fit <- readRDS("../scratch/qjep_model_fit.rda")

calc_trl_path <- function(obs, trl, cond) {
  
  # function to compute path length for a trial
  # trial is specified with obs(erver), trial and condition
  dtrl <- filter(d, 
                 observer == obs, trial == trl, condition == cond)  %>%
    mutate(x0 = lag(x),
           y0 = lag(y),
           d = sqrt((x - x0)^2 + (y - y0)^2)) %>%
    select(observer, condition, trial, id, targ_type,  x, y, d)
  
 
  
  # also compute the model's path length 
  # (allow it to start at the same item)
  
  d_trial = select(dtrl, trial, id, type = "targ_type", x,y) %>%
    mutate(found = -1, dist = NA, phi = NA )
  
  params = filter(fit, 
                  observer ==obs,
                  condition == cond) %>%
    mutate(spat_config = FALSE) %>%
    mutate(b_pref = boot::logit(pA),
           b_switch = boot::logit(pS),
           tune_prox = bP,
           tune_angle = bM,
           m = nrow(d_trial))
  
  sim_trial <- sim_momentum_trial(params, trial = trl, d_trial = d_trial,
                     init_sel = dtrl$id[1]) %>%
    arrange(found) %>%
    select(trial, x, y, found) %>% 
    mutate(x0 = lag(x),
           y0 = lag(y),
           d = sqrt((x - x0)^2 + (y - y0)^2))
  
  path_length = tibble(human = sum(dtrl$d, na.rm = TRUE),
                       simul = sum(sim_trial$d, na.rm = TRUE))
 
  return(path_length)
}

trial_list <- filter(d, found == 1) %>%
  select(obs = "observer", trl = "trial", cond = "condition")


trial_list$paths <- pmap_df(trial_list, calc_trl_path)

saveRDS(trial_list, "path_lengths.rda")
paths <- readRDS("path_lengths.rda")

trial_list %>% pivot_longer(c(human_path_length, model_path_length), 
                            names_to = "participant", values_to = "path_length") -> trial_list


ggplot(trial_list, aes(path_length, fill = participant)) + 
  geom_histogram(bins = 20, alpha = 0.5)


paths %>% mutate(diff = human_path_length - model_path_length) %>%
  ggplot(aes(y = diff, group = obs)) + geom_boxplot()

paths %>% mutate(diff = human_path_length - model_path_length) %>%
  group_by(obs) %>%
  summarise(median_diff = median(diff)) %>%
  ggplot(aes(median_diff)) + geom_histogram(bins = 20)

  