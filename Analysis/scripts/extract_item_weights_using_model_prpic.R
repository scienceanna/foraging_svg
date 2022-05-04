library(tidyverse)
library(rstan)
library(tidybayes)

# this function takes a dataset and fitted model, then steps through the data 
# and computes the weights assigned to each item by the model

source("../functions/get_model_params.R")
source("../functions/compute_weights.R")

# read in data
d <- read_csv("../data/squirrel_spatial.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"),
         targ_type = as.numeric(targ_type))

# read in model and extract weights
fit <- get_model_params("../scratch/squirrel.model") %>% rename(observer = "obs")
saveRDS(fit, "../scratch/prpic_model_fit.rda")

# run for all participants...
a_feat <- map_dfr(unique(d$observer), compute_weights_trials, cond = "feature", n_trials = 5)
a_conj <- map_dfr(unique(d$observer), compute_weights_trials, cond = "conjunction", n_trials = 5)
a <- bind_rows(a_feat, a_conj)
rm(a_feat, a_conj)

saveRDS(a, "../scratch/prpic_model_weights.rda")