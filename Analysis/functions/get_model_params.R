get_model_params <- function(filename) {
  
  
  cond_levels = c("feature", "conjunction")
  
  m <- readRDS(filename)
  
  m %>% spread_draws(b[cond], u[cond, obs]) %>%
    mutate(cond = as.factor(cond),
           cond = fct_recode(cond, bA_1 = "1", bS_1 = "2", bP_1 = "3", bM_1 = "4", 
                             bA_2 = "5", bS_2 = "6", bP_2 = "7", bM_2 = "8"), 
           bZ = b + u,
           pZ = boot::inv.logit(bZ)) %>%
    separate(cond, into = c("param", "condition")) %>%
    mutate(condition = factor(condition, labels = cond_levels)) -> fit1
  
  
  fit1 %>% select(obs, condition, param, pZ) %>%
    pivot_wider(names_from = "param", values_from = "pZ") %>%
    select(obs, condition, bA, bS) %>%
    pivot_longer(c('bA', 'bS'), names_to = "model", values_to = "bias") %>%
    unnest() -> qs
  
  fit1 %>% select(obs, condition, param, bZ) %>%
    pivot_wider(names_from = "param", values_from = "bZ") %>%
    select(obs, condition, bP, bM) %>%
    pivot_longer(c('bP', 'bM'), names_to = "model", values_to = "bias") %>%
    unnest() -> qs_p
  
  qs_all <- rbind(qs, qs_p) %>%
    mutate(model = fct_relevel(model, "bA", "bS", "bP", "bM"),
           model = plyr::revalue(model, c("bA" = "pA", "bS" = "pS"))) %>%
    group_by(obs, condition, model) %>%
    summarise(bias = median(bias)) %>%
    pivot_wider(names_from = "model", values_from = "bias")
  
  
  return(qs_all)
  
}
  
  