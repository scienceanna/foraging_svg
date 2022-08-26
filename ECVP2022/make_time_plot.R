

a %>% select(observer, condition, trial, found, selected_max) %>%
  full_join(d, by = c("observer", "condition", "trial", "found")) %>%
  arrange(observer, condition, trial, found) %>%
  mutate(
    dist = sqrt((x-lag(x))^2 + (y-lag(y))^2),
    RTd = RT - lag(RT),
    RTd = if_else(found == 1, 0, RTd),
    logRTd = log2(RTd)) %>%
  filter(is.finite(logRTd)) %>%
  mutate(
    time_bin = cut(logRTd,breaks = 25)) %>%
  separate(time_bin, into = c("binS", "binE"), ",") %>%
  mutate(binS = parse_number(binS),
         binE = parse_number(binE),
         logRTd = (binS + binE)/2,
         RTd = 2^logRTd) -> qq

%>%
  group_by(observer, condition, RTd) %>%
  summarise(acc = mean(selected_max),
            n = n()) %>%
  filter(n > 5) -> b


ggplot(b, aes(x = RTd, y = acc, colour = condition)) + 
  geom_jitter(height = 0, width = 0.01, alpha = 0.5) +
  scale_x_log10("inter-item selection time") + 
  scale_y_continuous("model accuracy") + 
  stat_smooth(method = glm, method.args = list(family = "binomial"), se = FALSE, alpha = 0.5)


###########




a %>% select(observer, condition, trial, found, selected_max) %>%
  full_join(d, by = c("observer", "condition", "trial", "found")) %>%
  arrange(observer, condition, trial, found) %>%
  mutate(
    dist = sqrt((x-lag(x))^2 + (y-lag(y))^2),
    logdist = log(dist),
    RTd = RT - lag(RT),
    RTd = if_else(found == 1, 0, RTd),
    logRTd = log2(RTd)) %>%
  filter(is.finite(logRTd)) %>%
  select(observer, condition, trial, found, selected_max, logdist, logRTd) -> qq




library(brms)
options(mc.cores = 6)

my_priors = c(prior(normal(0, 1), class = "Intercept"),
              prior(normal(0, 0.5), class = "b") )


m <- brm(selected_max ~ logdist*logRTd + (1|observer), 
         data = filter(qq, found == 10), 
         family = "bernoulli",
         chains = 4, cores = 4)

qq %>% modelr::data_grid(logRTd = seq_range(logRTd, 5), logdist = seq_range(logdist, 100)) %>%
  add_epred_draws(m, re_formula = NA) %>%
  mutate(logRTd = as_factor(logRTd)) %>%
  group_by(logRTd, logdist) %>%
  median_hdci() %>%
  ggplot(aes(x = logdist, y = .epred, ymax = .upper, ymin = .lower, fill = logRTd)) +
  geom_path() + 
  geom_lineribbon(alpha = 0.33)


