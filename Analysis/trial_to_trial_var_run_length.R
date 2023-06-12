# check trial-to-trial variability in sumamry stats

library(tidyverse)

d <- read_csv("../../foraging_spatial/data/clarke2020/clarke_2020_qjep.csv")

d %>% group_by(observer, condition, trial) %>%
  summarise(n =n()) %>%
  rename(pp = "observer", kk = "condition", tt = "trial") %>%
  select(-n) -> trials

get_run_info <- function(pp, kk, tt) {
  
  trl = filter(d, 
                 observer == pp,
                 condition == kk,
                 trial == tt)
    
    rl <- rle(trl$targ_type)
    
    n_runs = length(rl$lengths)
    max_run_length = max(rl$lengths)
    
dout <- tibble(
  observer = pp,
  condition = kk,
  trial = tt,
  n_runs = n_runs,
  max_rl = max_run_length
)
  
  return(dout)
  
}



druns <- pmap_df(trials, get_run_info)


# let's just look at conj
druns <- filter(druns, condition == 2)

druns %>% group_by(observer, condition) %>%
  summarise(meanruns = mean(n_runs), 
            sdruns = sd(n_runs)) -> dsd

dsd %>% ggplot(aes(meanruns, sdruns)) + geom_point()


druns %>% #filter(observer %in% hi_var$observer) %>%
  ggplot(aes(trial, n_runs, group = observer)) + 
  geom_hline(data = dsd, aes(yintercept = meanruns), linetype = 2) + 
  geom_path(alpha = 0.5, colour = "purple") + 
  facet_wrap(~observer) +
  theme(legend.position = "none") + 
  theme_minimal()

ggsave("n_runs_by_trial.pdf", width = 10, height = 10)

