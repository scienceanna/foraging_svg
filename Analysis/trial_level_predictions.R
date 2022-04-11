library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(ggrepel)

source("functions/get_model_params.R")
source("comp_weights.R")

d <- read_csv("data/clarke_2020_qjep.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"))

fit <- get_model_params() %>% rename(observer = "obs")

comp_trials <- function(obs, cond) {
  params <- filter(fit,  observer == obs, condition == cond )
  a <- map_dfr(1:20, comp_weights, params = params, cond = cond, obs = obs)
  return(a)
}

# a_feat <- map_dfr(unique(d$observer), comp_trials, cond = "feature")
# a_conj <- map_dfr(unique(d$observer), comp_trials, cond = "conjunction")
# a <- bind_rows(a_feat, a_conj)
# rm(a_feat, a_conj)

a <- readRDS("scratch/model_weights.rda")

a %>% group_by(condition, observer, found) %>%
  summarise(meanb = mean(b),
            prop_best = mean(bMax)) %>%
  filter(found > 1) %>%
  mutate(chance = 1/(41-found))-> a_agg

ggplot(a_agg, aes(x = found, y = meanb, colour = condition, fill = condition)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.25) + 
  stat_lineribbon(alpha = 0.25) + 
  geom_path(data = filter(a_agg, observer == 1, condition == "feature"), 
                          aes(y = chance), linetype = 2, colour = "black") + 
  theme_minimal() +
  ggthemes::scale_fill_fivethirtyeight() + 
  ggthemes::scale_colour_fivethirtyeight() +
  scale_x_continuous("target selection") + 
  scale_y_continuous("average weight from model") -> plt_b

ggplot(a_agg, aes(x = found, y = prop_best, colour = condition, fill = condition)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.25) + 
  stat_lineribbon(alpha = 0.25) + 
  geom_path(data = filter(a_agg, observer == 1, condition == "feature"), 
            aes(y = chance), linetype = 2, colour = "black") + 
  theme_minimal() +
  ggthemes::scale_fill_fivethirtyeight() + 
  ggthemes::scale_colour_fivethirtyeight() +
  scale_x_continuous("target selection") + 
  scale_y_continuous("proportion most likely was selected") -> plt_c

plt_b + plt_c + plot_layout(guides = "collect")
ggsave("../Figures/qjep_preds.png", width = 9, height = 4)

a %>% mutate(b_bin = cut(b, breaks = 100, labels = FALSE)) %>%
  group_by(b_bin) %>% 
  summarise(acc = mean(bMax)) %>%
  mutate(b_bin = as.numeric(b_bin)/100) %>%
  ggplot(aes(b_bin, acc)) + geom_path() +
  geom_abline(linetype = 2)



mistakes <- filter(a, observer == 1, trial == 3, condition == "conjunction", err> 0.5) %>%
  select(x, y, found, err, model_pref)

ggplot(filter(d, observer == 1, trial == 3, condition == "conjunction"), aes(x, y, colour = targ_type)) +
  geom_label(aes(label = found)) + geom_path() +
  geom_label_repel(data = mistakes, aes(x, y, label = paste(model_pref, ", p=", round(err,2))), colour = "red")



a %>% group_by(observer, condition) %>%
  summarise(acc = mean(bMax, na.rm = T)) %>%
  arrange((acc))


ggplot(filter(d, observer == 35, trial == 6, condition == "feature"), aes(x, y, colour = targ_type)) +
  geom_label(aes(label = found)) + geom_path() +
  geom_label_repel(data = mistakes, aes(x, y, label = paste(model_pref, ", p=", round(err,2))), colour = "red")


