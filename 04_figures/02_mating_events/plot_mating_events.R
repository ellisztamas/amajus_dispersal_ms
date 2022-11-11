library(tidyverse)
library(ggpubr)

source("03_analysis/04_mating_events/n_mating_events.R")

p1 <- me %>% 
  # Keep only mating events with a sampled father, and with at least two offspring
  # (offspring number is weighted by the probability of having occurred, so 
  # this may be less than 1.)
  filter(!is.na(father), offspring >=0) %>% 
  # Average over iterations
  group_by(mother, array_size) %>% 
  summarise(
    full_sibships = sum(prob) / 1000
  ) %>% 
ggplot( aes( x= array_size, y = full_sibships)) + 
  geom_point() + 
  geom_line(data = pframe, aes(x = array_size, y = n)) +
  labs(
    x = "Maternal family size",
    y = "Number of full sibships"
  ) +
  theme_bw()

p2 <- me %>% 
  # Keep only mating events with a sampled father
  filter( !is.na(father) ) %>% 
  group_by(mother, father) %>% 
  summarise(
    sibship_size = sum(offspring) / 1000,
    p = sum(prob) / 1000,
    p_lower = quantile(prob, 0.02),
    p_upper = quantile(prob, 0.98)
  ) %>% 
  ggplot(aes(y = p, x = sibship_size)) + 
  geom_point() +
  geom_errorbar( aes( ymin = p_lower, ymax = p_upper )) +
  labs(
    x = "Full-sibship size",
    y = "Mean posterior probability"
  ) +
  theme_bw()

png(width = 16.9, height = 8, units = "cm", res = 600,
  filename = "05/plot_mating_events.png"
)
ggarrange(p1, p2, ncol = 2, labels = 'AUTO')

dev.off()


me %>% 
  # Keep only mating events with a sampled father
  filter( !is.na(father) ) %>% 
  group_by(mother, father) %>% 
  summarise(
    sibship_size = sum(offspring) / 1000,
    p = sum(prob) / 1000,
    p_lower = quantile(prob, 0.02),
    p_upper = quantile(prob, 0.98)
  ) %>% 
  ggplot(aes(y = p, x = sibship_size)) + 
  geom_point() +
  geom_errorbar( aes( ymin = p_lower, ymax = p_upper )) +
  labs(
    x = "Full-sibship size",
    y = "Mean posterior probability"
  ) +
  theme_bw()
