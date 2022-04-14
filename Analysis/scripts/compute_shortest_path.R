library(igraph)
library(tidyverse)


d <- read_csv("../data/clarke_2020_qjep.csv") %>%
  mutate(condition = as_factor(condition),
         condition = fct_recode(condition, feature = "1", conjunction = "2"),
         targ_type = as.numeric(targ_type))


obs = 1
trl = 1
cond= 1


dtrl <- filter(d, observer == obs, trial == trl, condiition == cond)