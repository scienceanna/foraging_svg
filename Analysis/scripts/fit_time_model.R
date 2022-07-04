library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)

options(mc.cores = 4)

d <- read_csv("../scratch/time_per_trial.csv") 

my_priors <- c(prior(normal(0.6, 0.1), class = "b", coef = "conditionfeature"),
               prior(normal(0.6, 0.1), class = "b", coef = "conditionconjunction"),
               prior(normal(0.0, 0.1), class = "b", coef = "conditionfeature:time_taken_s"),
               prior(normal(0.0, 0.1), class = "b", coef = "conditionconjunction:time_taken_s"))

m <- brm(acc ~ 0 + condition + time_taken_s:condition + 
           (0 + condition +  time_taken_s:condition|observer),
         data = d,
         family = "normal",
         prior = my_priors)

saveRDS(m, "../scratch/time_model_acc_model.rds")

