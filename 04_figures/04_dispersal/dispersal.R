library(tidyverse)
library(cowplot)

# Import function to compute weighted cumulative density for dispersal distances
source("02_library/R/stat_ecdf_weighted.R")
# source("04_figures/02_posterior_distributions/summarise_posteriors.R")

# Import data on mating events
me <-read_csv(
  "03_analysis/04_mating_events/output/mating_events_over_chains.csv",
  col_types = cols()) %>% 
  filter( !is.na(father) ) # remove mating events with a missing father

# Histogram of disperal distances
linear_dispersal <- me %>%
  filter(offspring >=1) %>%
  mutate(iter = as.factor(iter)) %>% 
  ggplot( aes(x = distance, weights = prob) ) +
  geom_histogram(aes(
    y = stat(count / 1000) ),
    binwidth = 30, fill = 'white', colour = 'black',
    boundary=0
  ) + 
  labs(
    x = "Distance (m)",
    y = "Mating events"
  ) +
  theme_bw()

# Log cumulative distribution of disperal distances.
log_dispersal <- me %>% 
  filter(offspring >=1) %>%
  mutate(iter = as.factor(iter)) %>% 
  ggplot( aes(x = distance, weights = prob, groups = iter) ) +
  stat_ecdf(
    geom = "step",
    pad = FALSE
  ) +
  labs(
    x = element_blank(),#"Distance (m)",
    y = "Cumulative density"
  ) +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  theme_bw() +
  draw_plot_label("B")

# Inset cumulative distibution into histogram
inset_plot <- ggdraw() +
  draw_plot(linear_dispersal) + draw_plot_label("A") + 
  draw_plot(log_dispersal, x = 0.95, y = 0.95, width = 0.6, height = 0.6, hjust = 1, vjust = 1)

# Save to disk.
ggsave(filename = "05_manuscript/fig-dispersal.eps", 
       plot = inset_plot,
       width = 112, 
       height = 120,
       units = "mm",
       device = "eps")
  

# log_dispersal <- me %>% 
#   mutate(iter = as.factor(iter)) %>% 
#   ggplot( aes(x = distance, weights = prob) ) +
#   geom_histogram(aes(
#     y = stat(count / 1000) ),
#     bins = 20, fill = 'white', colour = 'black'
#     ) +
#   stat_ecdf(
#     aes(y=..y..*20),#, groups = iter),
#     geom = "step",
#     pad = FALSE
#   ) +
#   scale_y_continuous(
#     limits = c(0,20),
#     sec.axis=sec_axis(trans = ~./20 , name="Cumulative density", )
#     ) +
#   labs(
#     x = element_blank(),#"Distance (m)",
#     y = "Mating events"
#   ) +
#   scale_x_log10() +
#   annotation_logticks(sides = 'b') +
#   theme_bw()

# png(
#   filename = "04_figures/dispersal.png", 
#   units = "cm", height = 10, width = 11.2, 
#   res = 600
# )
# 
# plot_dispersal
# 
# dev.off()
