library(tidyverse)
library(rstan)
library(tidybayes)

options(digits = 2, mc.cores = 10)
d <- read_csv("data/clarke_2020_qjep.csv") %>%
  filter(found == 1)


stan_data <- list(N = nrow(d),
                  L = max(d$observer),
                  Z = d$observer,
                  x = d$x,
                  y = d$y)

m <- stan("models/initial_selection_models/init_sel.stan",
                data = stan_data,
                chains = 4,
                iter = 1000,
                refresh = 100)



saveRDS(m, "scratch/init_sel_model.rds")

gather_draws(m, u[dim, observer]) %>% 
  ungroup() %>%
  select(-.variable) %>%
  rename(obs_mu = ".value") %>%
  full_join(gather_draws(m, mu[dim])) %>%
  mutate(obs_mu = .value + obs_mu) -> a 


d %>% pivot_longer(c(x, y), names_to = "dim", values_to = "value") %>%
  mutate(dim = if_else(dim == "x", 1, 2)) -> d_plt

a %>% ggplot(aes(obs_mu, group = observer)) + geom_density() + 
  facet_grid(observer~dim, scales = "free") + 
  geom_density(data = d_plt, aes(value), colour = "purple")

%>%
  ggplot(aes(.value, group = observer)) + 
  geom_density(alpha = 0.1) + 
  facet_wrap(~.variable)
