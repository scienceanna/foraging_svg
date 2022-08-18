

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
  group_by(observer, condition, trial) %>%
  summarise(time_taken = max(RT),
            acc = mean(selected_max, na.rm = T), 
            .groups = "drop") %>%
  filter(time_taken < 60) %>%
  mutate(time_taken_s = as.numeric(scale(log(time_taken)))) -> timedat



a %>% select(observer, condition, trial, found, selected_max) %>%
  full_join(d, by = c("observer", "condition", "trial", "found")) %>%
  arrange(observer, condition, trial, found) %>%
  mutate(
    dist = sqrt((x-lag(x))^2 + (y-lag(y))^2),
    logdist = log(dist),
    RTd = RT - lag(RT),
    RTd = if_else(found == 1, 0, RTd),
    logRTd = log2(RTd)) %>%
  filter(is.finite(logRTd)) -> qq

%>%
  ggplot(aes(log(dist), logRTd)) + geom_point(alpha = 0.1)



library(brms)
options(mc.cores = 6)
m <- brm(selected_max ~ logdist*RTd*condition - logdist:RTd:condition + (1|observer), qq, 
         family = "binomial",
         chains = 4, cores = 4)



