library(tidyverse)
library(tidybayes)
library(rstan)
library(patchwork)
library(wesanderson)


plot_multi_levels <- function(fit, cond_levels = c("value", "no value"))
{
  fit %>% spread_draws(b[cond], u[cond, obs]) %>%
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

  ggplot(qs_all, aes(x = pA,  y = pS, colour = condition)) +
    geom_point(alpha = 0.5) + 
    scale_x_continuous("pA (preference for target A)", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous("pS (stick preference)") + 
    scale_colour_manual(values = c("#DDAA33", "#004488")) + 
    facet_wrap(~condition) + 
    theme_bw() -> pApS
  
  ggplot(qs_all, aes(x = pA, y =bP, colour = condition)) +
    geom_point(alpha = 0.5) +
    scale_x_continuous("pA (preference for target A)", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous("bP (proximity weighting)") +
    scale_colour_manual(values = c("#DDAA33", "#004488")) + 
    facet_wrap(~condition) + 
    theme_bw() -> pAbP
  
  ggplot(qs_all, aes(x = pA,  y =bM,  colour = condition)) +
    geom_point(alpha = 0.5) +  
    scale_x_continuous("pA (preference for target A)", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous("bM (direction weighting)") +
    scale_colour_manual(values = c("#DDAA33", "#004488")) + 
    facet_wrap(~condition) + 
    theme_bw() -> pAbM
  
  ggplot(qs_all, aes(x = pS, y = bP, colour = condition)) +
    geom_point(alpha = 0.5) +  
    scale_x_continuous("pS (stick preference)", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous("bP (proximity weighting)") +
    scale_colour_manual(values = c("#DDAA33", "#004488")) + 
    facet_wrap(~condition) + 
    theme_bw() -> pSbP
  
  ggplot(qs_all, aes(x = pS, y =bM, colour = condition)) +
    geom_point(alpha = 0.5) + 
    scale_x_continuous("pS (stick preference)",  limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous("bM (direction weighting)") + 
    scale_colour_manual(values = c("#DDAA33", "#004488")) + 
    facet_wrap(~condition) + 
    theme_bw() -> pSbM
  
  ggplot(qs_all, aes(x = bP, y =bM,  colour = condition)) +
    geom_point(alpha = 0.5) + 
    facet_wrap(~condition) + 
    scale_x_continuous("bP (proximity weighting)") + 
    scale_y_continuous("bM (direction weighting)") + 
    scale_colour_manual(values = c("#DDAA33", "#004488")) + 
    geom_smooth(method = "lm", linetype = 3, colour = "black", fullrange = TRUE) +
    theme_bw() -> bPbM
  

  
  qs_all %>%
    pivot_longer(-c(obs, condition), names_to = "param", values_to = "bias") %>%
    pivot_wider(names_from = "condition", values_from = "bias") %>%
    mutate(param = as_factor(param), 
           param = fct_relevel(param, "pA", "pS"),
           param = fct_recode(param, 
                              `pA (preference for target A)` = "pA",
                              `pS (stick preference)` = "pS",
                              `bP (proximity weighting)` = "bP",
                              `bM (direction weighting)` = "bM")) -> d_plt
  
  ggplot(d_plt, aes(feature, conjunction)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~param, scales = "free", nrow = 1) +
    geom_smooth(data = filter(d_plt, param %in% c("bP (proximity weighting)", "bM (direction weighting)")),
                method = "lm", linetype = 3, colour = "black", fullrange = TRUE) +
    theme_bw()  -> plt_cond
  
  
  plt <- (pApS + pAbP + pAbM)  / (pSbP +pSbM+bPbM )/ (plt_cond) + 
    plot_layout(guides="collect", heights = c(0.7,0.7,1)) & theme(legend.position = "none") 
  ggsave( "qjep_scatter.pdf", width = 10, height = 7)
  
}


spat_pred <- function(ii, fit, x_max=1) {
  x = seq(0, x_max, 0.01)
  p = exp(-fit$b[ii] * x)
  out <- tibble(condition = fit$condition[ii], 
                .draw = fit$.draw[ii],
                x = x, p = p)
  return(out)
}


plot_fixed_effects <- function(fit, cond_levels = c("feature", "conjunction"), flip_ab = FALSE) {
  
  fit %>% spread_draws(b[cond]) %>%
    mutate(cond = as.factor(cond),
           cond = fct_recode(cond, bA_1 = "1", bS_1 = "2", bP_1 = "3", bM_1 = "4", 
                             bA_2 = "5", bS_2 = "6", bP_2 = "7", bM_2 = "8")       ) %>%
    separate(cond, into = c("param", "condition")) %>%
    mutate(condition = factor(condition, labels = cond_levels),
           p = boot::inv.logit(b)) -> fit
  
  
  if (flip_ab) {
    fit %>% mutate(p = if_else(param == "bA", 1 - p, p)) -> fit
  }

  # get p(diff)
  fit %>%  
    pivot_wider(names_from = "condition",  values_from = "b") %>%
    mutate(diff = !!sym(cond_levels[[2]]) - !!sym(cond_levels[[1]])) %>% # 
    group_by(param) %>%
    summarise(p_diff = max(mean(diff>0), mean(diff<0))) -> fit_p_diff
  
  #print(fit_p_diff)
  
  facet_names <- c(
    bA = "pA (preference for target A)",
    bS = "pS (stick preference)"
    
  )
  
  
 fit %>%
    filter(param %in% c("bA", "bS")) %>%
    ggplot(aes(x = p,  fill = condition)) + 
    geom_density(alpha = 0.5) +
   facet_grid(~param, scales="free", labeller = as_labeller(facet_names)) + 
    scale_fill_manual(values = c("#DDAA33", "#004488")) + 
   scale_x_continuous("weight", limits = c(0.4, 1), breaks = seq(0.5, 1, 0.25)) +
    geom_vline(xintercept = 0.5, linetype = 2) + 
    theme_bw() + 
   theme(legend.position = "right") -> plt_AS
  
  f <- filter(fit, param == "bP") 
  
  map_dfr(1:nrow(f), spat_pred, fit=f, 1) %>% 
    group_by(condition, x) %>%
    ggplot(aes(x, y = p, fill= condition)) +
    stat_lineribbon( .width = c(0.53, 0.89, 0.97), alpha = 0.4, linetype = 0) +
    theme_bw() +
    scale_fill_manual(values = c("#DDAA33", "#004488")) + 
    scale_x_continuous("distance (prop. width)", limits = c(0, 1), breaks = c(0, 0.5, 1)) + 
    geom_hline(yintercept = 1, linetype = 2) +
    scale_y_continuous("proximity weighting") +
    theme(legend.position = "none") -> plt_P
  
  f<- filter(fit, param == "bM") 
  
  map_dfr(1:nrow(f), spat_pred,  filter(f, param == "bM")) %>% 
    group_by(condition, x) %>%
    ggplot(aes(x, y = p, fill = condition)) +
    stat_lineribbon( .width = c(0.53, 0.89, 0.97), alpha = 0.4, linetype = 0) +
    theme_bw() +
    scale_fill_manual(values = c("#DDAA33", "#004488")) + 
    scale_x_continuous("direction", breaks = c(0, 0.5, 1), labels = c(0, "pi/2", expression("pi"))) + 
    scale_y_continuous("direction weighting") +
    geom_hline(yintercept = 1, linetype = 2) +
    theme(legend.position = "none") -> plt_D
  
  plt <- plt_AS / (plt_P + plt_D)# +  plot_layout(guides = 'collect') & theme(legend.position = 'top')
  return(plt)
}
  
