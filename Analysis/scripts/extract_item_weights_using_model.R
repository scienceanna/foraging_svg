library(tidyverse)
library(rstan)
library(tidybayes)


# this function takes a dataset and fitted model, then steps through the data 
# and computes the weights assigned to each item by the model

source("../functions/get_model_params.R")
source("../functions/compute_weights.R")

# read in data
d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"))

# read in model and extract weights
fit <- get_model_params("../scratch/all_qjep_2020.rds") %>% rename(observer = "obs")
saveRDS(fit, "../scratch/qjep_model_fit.rda")

# run for all participants...
a_feat <- map_dfr(unique(d$observer), comp_trials, cond = "feature")
a_conj <- map_dfr(unique(d$observer), comp_trials, cond = "conjunction")
a <- bind_rows(a_feat, a_conj)
rm(a_feat, a_conj)

saveRDS(a, "../scratch/qjep_model_weights.rda")


