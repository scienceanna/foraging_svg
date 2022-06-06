library(tidyverse)
library(rstan)
library(tidybayes)
library(truncnorm)

options(digits = 2, mc.cores = 10)
d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  filter(found == 1) %>%
  mutate(x = if_else(x < 0.01, 0.01, x),
         x = if_else(x > 0.99, 0.99, x),
         y = if_else(y < 0.01, 0.01, y),
         y = if_else(y > 0.99, 0.99, y)
         )


stan_data <- list(N = nrow(d),
                  L = max(d$observer),
                  Z = d$observer,
                  x = d$x,
                  y = d$y)

m <- stan("../models/initial_selection_models/init_sel_beta.stan",
                data = stan_data,
                chains = 1,
                iter = 1000,
                refresh = 100)



saveRDS(m, "../scratch/init_sel_model.rds")

m <- readRDS("../scratch/init_sel_model.rds")

gather_draws(m, u[dim, observer]) %>% 
  ungroup() %>%
  select(-.variable) %>%
  rename(obs_b = ".value") %>%
  left_join(gather_draws(m, b[dim]), 
            by = c("dim", ".chain", ".iteration", ".draw")) %>%
  mutate(obs_b = exp(.value + obs_b),
         .variable = if_else(dim %in% c(1,2), "a", "b"),
         dim = if_else(dim %in% c(1,3), "x", "y")) %>%
  select(observer, dim, .iteration, obs_b, param = ".variable") -> a 

ggplot(a, aes(obs_b, group = observer)) + 
  geom_density(alpha = 0.25) +
  facet_grid(param ~ dim, scales = "free")
 
a %>% group_by(observer, param, dim) %>%
  summarise(value = mean(obs_b) , .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = value) -> am

compute_beta_dist <- function(observer, a, b) {
  x <- seq(0.01,0.99,0.01)
  return(tibble(observer = observer, 
                x=  x,
                z = dbeta(x, a, b)))

}

obs_fits_x <- pmap_df(filter(am,  dim == "x") %>% select(observer, a, b), 
                      compute_beta_dist) %>%
  ggplot(aes(x = x, y = z, group = observer)) + geom_path() +
  scale_x_continuous("x")

obs_fits_x <- pmap_df(filter(am,  dim == "x") %>% select(observer, a, b), 
                      compute_beta_dist) %>%
  ggplot(aes(x = x, y = z, group = observer)) + geom_path() +
  scale_x_continuous("y")


obs_fits_x + obs_fits_y



m <- stan("../models/initial_selection_models/init_sel2_beta.stan",
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
