library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(ggrepel)

source("functions/get_model_params.R")
source("functions/compute_weights.R")



# read in data
d <- read_csv("data/clarke_2020_qjep.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"))

# read in model and extract weights
fit <- get_model_params("scratch/all_qjep_2020.rds") %>% rename(observer = "obs")


# run for all participants...
a_feat <- map_dfr(unique(d$observer), comp_trials, cond = "feature")
a_conj <- map_dfr(unique(d$observer), comp_trials, cond = "conjunction")
a <- bind_rows(a_feat, a_conj)
rm(a_feat, a_conj)

saveRDS(a, "scratch/qjep_model_weights.rda")
a <- readRDS("scratch/qjep_model_weights.rda")



a %>% group_by(condition, observer, found) %>%
  summarise(meanb = mean(b),
            prop_best = mean(selected_max)) %>%
  filter(found > 1) %>%
  mutate(chance = 1/(41-found))-> a_agg

# plot target selected weights
ggplot(a_agg, aes(x = found, y = meanb, colour = condition, fill = condition)) + 
  geom_jitter(data = filter(a_agg, found<40), width = 0.1, height = 0, alpha = 0.2) + 
  stat_lineribbon(.width = 0.67, alpha = 0.50) +
  geom_path(data = filter(a_agg, observer == 1, condition == "feature"), 
                          aes(y = chance), linetype = 2, colour = "black") + 
  geom_point(data = tibble(x=40, y=1), aes(x, y), size = 1.5, colour = "black", fill = "grey") + 
  scale_x_continuous(breaks = c(2, 10, 20, 40), "target selection", expand = c(0.01, 0.01)) + 
  scale_y_continuous("average weight from model", expand = c(0.01, 0.01)) -> plt_b

ggplot(a_agg, aes(x = found, y = prop_best, colour = condition, fill = condition)) + 
  geom_jitter(data = filter(a_agg, found<40), width = 0.1, height = 0.00, alpha = 0.2) + 
  stat_lineribbon(.width = 0.67, alpha = 0.50) +
  geom_path(data = filter(a_agg, observer == 1, condition == "feature"), 
            aes(y = chance), linetype = 2, colour = "black") + 
  geom_point(data = tibble(x=40, y=1), aes(x, y), size = 1.5, colour = "black", fill = "grey") + 
  scale_x_continuous(breaks = c(2, 10, 20, 40), "target selection", expand = c(0.01, 0.01)) + 
  scale_y_continuous("proportion most likely was selected", expand = c(0.01, 0.01)) -> plt_c

plt_b + plt_c + plot_layout(guides = "collect")  &
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')
ggsave("../Figures/qjep_preds.png", width = 9, height = 4)

a %>% mutate(b_bin = cut(max_b, breaks = 100, labels = FALSE)) %>%
  group_by(condition, b_bin) %>% 
  summarise(acc = mean(selected_max)) %>%
  mutate(b_bin = as.numeric(b_bin)/100) %>%
  filter(is.finite(acc)) %>%
  ggplot(aes(b_bin, acc, colour = condition)) + 
  geom_point() + 
  geom_line(stat = "smooth", method = "loess", alpha = 0.5, size = 2, se = F) +
  geom_abline(linetype = 2)



mistakes <- filter(a, observer == 1, trial == 3, condition == "conjunction", err> 0.75) %>%
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


