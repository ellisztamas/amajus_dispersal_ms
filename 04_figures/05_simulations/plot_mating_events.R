#' Plot patterns of simulated mating events.

library(tidyverse)
library(ggpubr)

sim_mating_events <- read_csv(
  "03_analysis/05_simulations/output/sims_mating.csv",
  col_types = 'ffidifddddddddd'
  )

plot_all_me <- sim_mating_events %>% 
  group_by(nloci, prop_purged,  data_type) %>% 
  summarise(
    `Less than one`  = mean(true_if_offs_ge1, na.rm = TRUE),
    `One or more`  = mean(arith_prob_correct, na.rm = TRUE)
  ) %>% 
  # Rejig the variables to make them plot nicely.
  mutate(
    # scale = paste0(scale, "m"),
    prop_purged = paste0(prop_purged*100, "% missing fathers"),
    # family_size = paste0(family_size, " siblings"),
    data_type = case_when(
      data_type == "paternity" ~ "Paternity only",
      data_type == "sibships" ~ "Pat., sibships",
      data_type == "full" ~ "Pat., sibs., dispersal"
    )
  ) %>% 
  rename(
    Inference = data_type
  ) %>% 
  pivot_longer(`Less than one` : `One or more`, names_to = "Family size", values_to = 'true_pos') %>% 
  ggplot(aes(x=nloci, y=true_pos, colour=Inference, linetype=`Family size`)) +
  geom_line() +
  geom_point() +
  labs(
    # title = "All mating events",
    x = "Number of markers",
    y = "Probability correct"
  ) +
  theme_bw() +
  guides(colour="none")+
  facet_grid(~prop_purged)

  


sim_mating_events <- sim_mating_events %>% 
  group_by(scale, nloci, prop_purged, family_size, data_type) %>% 
  summarise(
    true_pos  = mean(true_if_offs_ge1, na.rm = TRUE),
    false_neg = mean(true_if_offs_lt1, na.rm = TRUE),
    missing = mean(missing, na.rm = TRUE)
  ) %>% 
  # Rejig the variables to make them plot nicely.
  mutate(
    scale = paste0(scale, "m"),
    prop_purged = paste0(prop_purged*100, "% missing fathers"),
    family_size = paste0(family_size, " siblings"),
    data_type = case_when(
      data_type == "paternity" ~ "Paternity only",
      data_type == "sibships" ~ "Pat., sibships",
      data_type == "full" ~ "Pat., sibs., dispersal"
    )
  ) %>% 
  rename(
    `Mean dispersal` = scale,
    Inference = data_type
  )

plot_ge1 <- sim_mating_events %>% 
  ggplot(aes(x=nloci, y=true_pos, colour=Inference, linetype=`Mean dispersal`)) +
  geom_line() +
  geom_point() +
  labs(
    # title = "One or more offspring",
    x = "Number of markers",
    y = "Probability correct"
  ) +
  theme_bw() +
  facet_grid(family_size~prop_purged)

plot_lt1 <- sim_mating_events %>% 
  ggplot(aes(x=nloci, y=false_neg, colour=Inference, linetype=`Mean dispersal`)) +
  geom_line() +
  geom_point() +
  labs(
    # title = "Less than one offspring",
    x = "Number of markers",
    y = "Probability correct"
  ) +
  theme_bw() +
  facet_grid(family_size~prop_purged)

plot_mating_1 <- ggarrange(
  plot_lt1, plot_ge1,
  nrow=2, common.legend = TRUE, legend="right", labels = c("B", "C")
  )
ggarrange(
  plot_all_me, plot_mating_1, 
  nrow=2, heights = c(1,3), labels = c("A", "")
)

ggsave(
  filename = "05_manuscript/fig-sim_mating.eps",
  plot = plot_mating,
  units = "cm", width = 16.9, height = 20
  )
