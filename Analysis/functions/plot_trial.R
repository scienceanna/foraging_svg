mistake_x_y <- function(ii, trl, mistakes) {
  mistake <- mistakes[ii, ]
  
  t1 <- filter(trl, found == mistake$found-1)
  t2 <- filter(trl, found == mistake$model_pref)
  
  dout <- tibble(id = ii, 
                 x1 = t1$x, y1 = t1$y,
                 x2 = t2$x, y2 = t2$y)
  return(dout)
}



plot_trial <- function(id, trls, as, d) {
  trl_id <- trls[id,]
  
  trl <- filter(d, 
                observer == trl_id$observer, 
                condition == trl_id$condition, 
                trial == trl_id$trial)
  
  mistakes <- filter(as, 
                      observer == trl_id$observer, 
                      condition == trl_id$condition, 
                      trial == trl_id$trial) %>% 
    filter(selected_max == 0)
  
  mistakes_xy <- map_df(1:nrow(mistakes), mistake_x_y, trl, mistakes)
  
  p <- ggplot(trl, aes(x, y)) +
    geom_text_repel(aes(label = found)) + 
    geom_path(colour = "darkgrey") +
    geom_point(aes(colour = targ_type, shape = targ_type), size = 4) +
    geom_segment(data = mistakes_xy, size = 1,
                 aes(x = x1, y = y1, xend = x2, yend = y2), colour = "darkred",
                 arrow = arrow(length = unit(0.1, "inches"))) + 
    theme_void() + theme(legend.position = "none") + 
    coord_fixed() + 
    ggtitle(paste0("model accuracy = ", round(trl_id$prop_max, 2))) 
  
  return(p)
 
 
}
    
