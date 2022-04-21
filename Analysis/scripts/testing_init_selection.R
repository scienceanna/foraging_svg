library(tidyverse)
library(rstan)
library(tidybayes)
library(truncnorm)

options(digits = 2, mc.cores = 10)
d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  filter(found == 1)


stan_data <- list(N = nrow(d),
                  L = max(d$observer),
                  Z = d$observer,
                  x = d$x,
                  y = d$y)

m <- stan("../models/initial_selection_models/init_sel.stan",
                data = stan_data,
                chains = 4,
                iter = 1000,
                refresh = 100)



saveRDS(m, "../scratch/init_sel_model.rds")

m <- readRDS("../scratch/init_sel_model.rds")

gather_draws(m, u[dim, observer]) %>% 
  ungroup() %>%
  select(-.variable) %>%
  rename(obs_b = ".value") %>%
  full_join(gather_draws(m, b[dim])) %>%
  mutate(obs_b = .value + obs_b,
         .variable = if_else(dim %in% c(1,2), "mu", "sd"),
         dim = if_else(dim %in% c(1,3), "x", "y")) %>%
  select(observer, dim, .iteration, obs_b, param = ".variable") %>%
  mutate(obs_b = if_else(param == "sd", exp(obs_b), obs_b)) -> a 

ggplot(a, aes(obs_b, group = observer)) + 
  geom_density(alpha = 0.25) +
  facet_grid(param ~ dim, scales = "free")
 
a %>% group_by(observer, param, dim) %>%
  summarise(value = mean(obs_b) , .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = value) %>%
  rename(mean = "mu") -> am

obs_fits_x <- pmap_dfc(filter(am,  dim == "x") %>% select(mean, sd), dtruncnorm, a = 0, b= 1, x = z) %>%
  mutate(z = z) %>%
  pivot_longer(-z, names_to = "observer", values_to = "v") %>%
  ggplot(aes(x = z, y = v, group = observer)) + geom_path() +
  scale_x_continuous("x")

obs_fits_y <- pmap_dfc(filter(am,  dim == "y") %>% select(mean, sd), dtruncnorm, a = 0, b= 1, x = z) %>%
  mutate(z = z) %>%
  pivot_longer(-z, names_to = "observer", values_to = "v") %>%
  ggplot(aes(x = z, y = v, group = observer)) + geom_path() +
  scale_x_continuous("y")


obs_fits_x + obs_fits_y



m <- stan("../models/initial_selection_models/init_sel2.stan",
          data = stan_data,
          chains = 4,
          iter = 1000,
          refresh = 100)


