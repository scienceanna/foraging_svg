library(tidyverse)
library(rstan)
library(tidybayes)
library(truncnorm)

options(digits = 2, mc.cores = 10)
d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  filter(found == 1) %>%
  mutate(x = if_else(x == 0, 0.01, x),
         x = if_else(x == 1, 0.99, x),
         y = if_else(y == 0, 0.01, y),
         y = if_else(y == 1, 0.99, y)
         )


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
          chains = 1,
          iter = 1000,
          refresh = 100,
          init = list(list(sd_x = c(1, 1))))


gather_draws(m, lambda[observer]) %>%
  ggplot(aes(.value, group = observer)) + 
  geom_density(alpha = 0.1, fill = "cyan") -> plt_all

gather_draws(m, lambda[observer]) %>%
  group_by(observer) %>%
  summarise(lambda = mean(.value)) %>% 
  filter(lambda < 0.1) -> ll_peeps

gather_draws(m, lambda[observer]) %>%
  group_by(observer) %>%
  summarise(lambda = mean(.value)) %>% 
  filter(lambda > 0.4, lambda < 0.8) -> ml_peeps

gather_draws(m, lambda[observer]) %>%
  group_by(observer) %>%
  summarise(lambda = mean(.value)) %>% 
  filter(lambda > 0.9) -> hl_peeps


ggplot(filter(d, observer %in% ll_peeps$observer), aes(x, y)) + geom_hex(aes(colour = ..count..), bins = 12) +
  scale_fill_viridis_c() + scale_color_viridis_c() +
  my_theme -> plt_ll


ggplot(filter(d, observer %in% ml_peeps$observer), aes(x, y)) + geom_hex(aes(colour = ..count..), bins = 12) +
  scale_fill_viridis_c() + scale_color_viridis_c() +
  my_theme -> plt_ml

ggplot(filter(d, observer %in% hl_peeps$observer), aes(x, y)) + geom_hex(aes(colour = ..count..), bins = 12) +
  scale_fill_viridis_c() + scale_color_viridis_c() +
  my_theme -> plt_hl


plt_all / (plt_ll + plt_ml + plt_hl)

m <- stan("../models/initial_selection_models/init_sel2_beta.stan",
          data = stan_data,
          chains = 1,
          iter = 1000,
          refresh = 100,
          init = list(list(sd_x = c(1, 1))))

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
