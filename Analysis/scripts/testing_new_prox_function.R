library(tidyverse)
library(rstan)

source("../functions/prep_stan_data_momentum.R")

d <- read_csv("../reanalysis_data/clarke_2020_qjep.csv") %>%
 # mutate(condition = as_factor(condition),
 #        condition = fct_recode(condition, feature = "1", conjunction = "2")) %>% 
  filter(observer == 1) 


n_targets <- 40
n_trials = max(d$trial)
prior_b = 1

stan_data <- prep_stan_data(d)
stan_model_name <- "momentum_foraging"


fit_ids <- stan(paste("../models/momentum_models/", stan_model_name, ".stan", sep = ""),
                data = stan_data,
                chains = 5,
                iter = 1000,
                refresh = 10)
