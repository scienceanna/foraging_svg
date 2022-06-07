library(tidyverse)
library(rstan)
library(tidybayes)

# this function takes a dataset and fitted model, then steps through the data 
# and computes the weights assigned to each item by the model

source("../functions/get_model_params.R")
source("../functions/compute_weights.R")

# read in data
d <- read_csv("../data/test_arni.csv") %>% # using test data (rather than full data)
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"),
         targ_type = as.numeric(targ_type))

# read in model and extract weights
fit <- get_model_params("../scratch/arni2014_train.model") %>% rename(observer = "obs")
saveRDS(fit, "../scratch/kristjansson_model_fit_train.rda")

# run for all participants...
a_feat <- map_dfr(unique(d$observer), compute_weights_trials, cond = "feature")
a_conj <- map_dfr(unique(d$observer), compute_weights_trials, cond = "conjunction")
a <- bind_rows(a_feat, a_conj)
rm(a_feat, a_conj)

saveRDS(a, "../scratch/kristjansson_model_weights_train.rda")
