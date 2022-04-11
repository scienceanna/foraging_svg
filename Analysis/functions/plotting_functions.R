library(RColorBrewer)
library(wesanderson)
library(latex2exp)
library(ggridges)

plot_a_trial <- function(trial, n_col = 20, colour_scheme = c("#FC8D62", "#66C2A5"), shapes = 18)
  
{
  n <- length(trial)
  
  n_row = ceiling(n/n_col)
  
  # NA bad trial so that it is a multiple of n_col * n_row
  trial <- c(trial, rep(NA, n_col*n_row - n))
  
  labs <- paste("targets", n_col*(seq(n_row, 1, -1))-(n_col-1), "-", seq(n_row, 1, -1)*n_col)
  
  tibble(x = rep(1:n_col, n_row), 
         y = rep(seq(n_row, 1, -1), each = n_col),
         target = as.factor(trial)) %>%
    ggplot(aes(x, y, colour = target, fill = target)) +
    geom_hline(yintercept = 1:n_row, linetype = 2) +
    geom_point(size = 4, shape = shapes) +
    scale_y_continuous(breaks = 1:n_row, labels = labs) +
    scale_color_manual(values = colour_scheme) +
    scale_fill_manual(values = colour_scheme) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_blank()) -> plt
  return(plt)
}

plot_model_dists <- function(d_plt) {
  
  ggplot(d_plt, aes(x = bias, colour = model, fill = model)) +
    geom_density(alpha = 0.5) + 
    geom_vline(xintercept = 0.5, linetype = 2) + 
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    ylab("") + 
    scale_fill_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal1")) +
    scale_colour_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal1")) +
    theme(#legend.position = c(0.15, 0.65),
          #plot.margin = margin(2,20,20,1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) 
  
}

plot_simple_one_trial_model <- function(fit) {
  
  extract(fit, c("p_prior", "pS")) %>% 
    as_tibble() %>%
    rename(prior = "p_prior", posterior = "pS") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    mutate(model = as_factor(model), model = fct_relevel(model, "prior")) -> d_plt
  
  plot_model_dists(d_plt)
  
}

plot_simple_one_condition_model <- function(fit) {
  
  extract(fit, c("pA_prior", "pS_prior", "pS[1]", "pA[1]")) %>% 
    as_tibble() %>%
    rename(prior_a = "pA_prior", prior_s = "pS_prior", p_s = "pS[1]", p_a = "pA[1]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    mutate(model = as_factor(model), model = fct_relevel(model, "prior")) -> d_plt
  
  plot_model_dists(d_plt)
  
}


plot_simple_two_condition_model <- function(fit, cond_labels, gt) {
  
  extract(fit, c("pS[1]", "pS[2]", "pA[1]", "pA[2]")) %>% 
    as_tibble() %>%
    mutate(pa_diff = `pA[2]` - `pA[1]`, ps_diff = `pS[2]` - `pS[1]`) %>%
    rename( ps_1 = "pS[1]", ps_2  = "pS[2]", pa_1 = "pA[1]", pa_2  = "pA[2]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    separate(model, into = c("model", "condition"))   %>%
    mutate(
      model = as_factor(model),
      condition = as_factor(condition),
      condition = fct_relevel(condition, "diff")) -> d_plt
  
  levels(d_plt$condition) <- cond_labels
  
  ggplot(d_plt, aes(x = bias, y = model, fill = condition)) +
    geom_density_ridges(alpha = 0.5) + 
    geom_vline(xintercept = gt, linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 1) + 
    scale_x_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
    scale_fill_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal1")) +
    theme(plot.margin = margin(2,20,20,1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) 
} 

plot_simple_two_condition_model_nodiff <- function(fit, cond_labels, gt) {
  
  extract(fit, c("pS[1]", "pS[2]", "pA[1]", "pA[2]")) %>% 
    as_tibble() %>%
    rename( ps_1 = "pS[1]", ps_2  = "pS[2]", pa_1 = "pA[1]", pa_2  = "pA[2]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    separate(model, into = c("model", "condition"))   %>%
    mutate(
      model = as_factor(model),
      condition = as_factor(condition),
      condition = fct_relevel(condition, "diff")) -> d_plt
  
  levels(d_plt$condition) <- cond_labels
  
  ggplot(d_plt, aes(x = bias, y = model, fill = condition)) +
    geom_density_ridges(alpha = 0.5) + 
    geom_vline(xintercept = gt, linetype = 2) + 
    #scale_x_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
    scale_fill_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("GrandBudapest2")) +
    theme(plot.margin = margin(2,20,20,1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) 
} 

plot_simple_twocondition_model_facets <- function(fit, cond_labels, p_lines, palette_choice) {
  
  
  extract(fit, c("pS[1]", "pS[2]", "pA[1]", "pA[2]")) %>% 
    as_tibble() %>%
    rename( ps_1 = "pS[1]", ps_2  = "pS[2]", pa_1 = "pA[1]", pa_2  = "pA[2]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    separate(model, into = c("model", "condition"))   %>%
    mutate(
      model = as_factor(model),
      condition = as_factor(condition)) -> d_plt

  
  facet_names <- c(
    pa = "pA (preference for target A)",
    ps = "pS (stick preference)"

  )
  
  ggplot(d_plt, aes(x = bias, fill = condition)) +
    facet_wrap(~factor(model, levels = c('pa', 'ps')), labeller = as_labeller(facet_names)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = palette_choice, labels = cond_labels) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none") +
    ylab("") +
    xlab("parameter value") +
    geom_vline(data = p_lines, aes(xintercept = line_coords), linetype = 2) 


}


print_two_condition_summary_table <- function(fit) {
  extract(fit, c("pS[1]", "pS[2]", "pA[1]", "pA[2]")) %>%
    as_tibble() %>%
    mutate(pa_diff = `pA[2]` - `pA[1]`, ps_diff = `pS[2]` - `pS[1]`) %>%
    rename( ps_1 = "pS[1]", ps_2  = "pS[2]", pa_1 = "pA[1]", pa_2  = "pA[2]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    separate(model, into = c("param", "condition"))  -> d_fit
  
  d_fit %>% group_by(param, condition) %>%
    median_hdci(.width = c(0.53, 0.97)) %>% knitr::kable()
  
  
  d_fit %>% filter(condition == "diff") %>% group_by(param) %>%
    summarise(p_diff_greater_zero = mean(bias>0)) %>% knitr::kable()
}

plot_simple_three_condition_model <- function(fit, cond_labels, gt) {
  
  extract(fit, c("pS[1]", "pS[2]", "pS[3]", "pA[1]", "pA[2]", "pA[3]")) %>% 
    as_tibble() %>%
    rename( ps_1 = "pS[1]", ps_2  = "pS[2]", ps_3 = "pS[3]", pa_1 = "pA[1]", pa_2  = "pA[2]", pa_3 = "pA[3]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    separate(model, into = c("model", "condition"))   %>%
    mutate(
      model = as_factor(model),
      condition = as_factor(condition)) -> d_plt
  
  levels(d_plt$condition) <- cond_labels
  
  ggplot(d_plt, aes(x = bias, y = model, fill = condition)) +
    geom_density_ridges(alpha = 0.5) + 
    geom_vline(xintercept = gt, linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 1) + 
    scale_x_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
    scale_fill_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal2")) +
    theme(plot.margin = margin(2,20,20,1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) 
} 

print_three_condition_summary_table <- function(fit) {
  extract(fit, c("pS[1]", "pS[2]", "pS[3]", "pA[1]", "pA[2]", "pA[3]")) %>%
    as_tibble() %>%
    rename(ps_1 = "pS[1]", ps_2  = "pS[2]", ps_3 = "pS[3]", pa_1 = "pA[1]", pa_2  = "pA[2]", pa_3 = "pA[3]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    separate(model, into = c("param", "condition"))  -> d_fit
  
  d_fit %>% group_by(param, condition) %>%
    median_hdci(.width = c(0.53, 0.97)) %>% knitr::kable()
  
}

plot_multilevel <- function(fit, cond_labels){
  
  fit %>% spread_draws(bA[condition], 
                       bAz[observer, condition],
                       bS[condition], 
                       bSz[observer, condition]) %>%
    mutate(bAz = bA  + bAz,
           bSz = bS  + bSz,
           pAz = boot::inv.logit(bAz),
           pSz = boot::inv.logit(bSz)) %>%
    select(observer, condition, pAz, pSz) %>%
    pivot_longer(c('pAz', 'pSz'), names_to = "model", values_to = "bias") -> q
  
  ggplot(q, aes(x = bias, y = as.factor(observer), colour = as.factor(condition))) + 
    stat_pointinterval() +
    #geom_density_ridges(alpha = 0.7, colour = "white", scale = 5, panel_scaling = FALSE) + 
    facet_wrap(~model) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(colour = "grey90")) +
    #scale_fill_manual(values = wes_palette("GrandBudapest2"), name = "condition") +
    scale_colour_manual(values = wes_palette("GrandBudapest2"), name = "condition", labels = cond_labels) +  
    ylab('participant')
  
}

plot_simple_prox_model <- function(fit, bA, bS, tune_prox) {
  
  extract(fit, c("b[1]", "b[2]", "b[3]", "b[4]")) %>% 
    as_tibble() %>%
    rename(directionTune = "b[4]", proxTune = "b[3]", b_s = "b[2]", b_a = "b[1]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    mutate(model = as_factor(model), model = fct_relevel(model, "proxTune", after = Inf)) -> d_plt
  
  ggplot(d_plt, aes(x = bias, fill = model, colour = model)) +
    geom_density(alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    scale_x_continuous(limits = c(-5, 20), expand = c(0, 0)) + 
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    scale_colour_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal1")) +
    scale_fill_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal1")) +
    theme(legend.position = c(0.65, 0.65),
          plot.margin = margin(2,20,20,1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    geom_vline(xintercept = tune_prox, linetype = 2) +
    geom_vline(xintercept = bA, linetype = 2) +
    geom_vline(xintercept = bS, linetype = 2) +
    geom_vline(xintercept = 0, colour = "red")
  
}

plot_simple_one_condition_model_noprior <- function(fit) {
  
  extract(fit, c("bS[1]", "bA[1]")) %>% 
    as_tibble() %>%
    rename(b_s = "bS[1]", b_a = "bA[1]") %>%
    pivot_longer(everything(), names_to = "model", values_to = "bias") %>%
    mutate(model = as_factor(model))-> d_plt
  
  ggplot(d_plt, aes(x = bias, fill = model)) +
    geom_density(alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    scale_x_continuous(limits = c(-0.5, 3), expand = c(0, 0)) + 
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    scale_fill_manual(labels = TeX(levels(d_plt$model)), values = wes_palette("Royal1")) +
    theme(legend.position = c(0.5, 0.65),
          plot.margin = margin(2,20,20,1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) 
  
}
