library(tidyverse)
library(rstan)
library(tidybayes)
library(truncnorm)

options(digits = 2, mc.cores = 10)

stan_list <- readRDS("../scratch/clarke2020qjep_stan_data.rds")


m <- stan("../models/initial_selection_models/momentum_foraging_multi_level.stan",
          data = stan_data,
          chains = 1,
          iter = 1000,
          refresh = 100,
          init = list(list(sd_x = c(1, 1))))

saveRDS(m, "../scratch/init_sel_model2.rds")

sample_beta <- function(c, a_y, b_y) {
  x <- seq(0.01, 0.99, 0.01)
  return(tibble(c = c, 
                x = x,
                y = dbeta(x, shape1 = a_y, shape2 = b_y)))
}

gather_draws(m, a_y[c], b_y[c]) %>%
  group_by(c, .variable) %>%
  summarise(.value = mean(.value)) %>%
  pivot_wider(names_from = .variable, values_from = ".value") %>%
  pmap_df(sample_beta) %>%
  mutate(c = as_factor(c)) %>%
  ggplot(aes(x, y)) +
  geom_path() + 
  facet_wrap(~c)


gather_draws(m, lambda[observer]) %>%
  ggplot(aes(.value, group = observer)) + 
  geom_density(alpha = 0.25, fill = "grey")
