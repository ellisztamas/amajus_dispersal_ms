library(tidyverse)
library(ggpubr)

sim_dispersal <- read_csv(
  "03_analysis/05_simulations/output/sims_dispersal.csv",
  col_types = 'ffidifd'
)

sim_dispersal <- sim_dispersal %>% 
  group_by(scale, nloci, prop_purged, family_size, data_type) %>% 
  summarise(
    deviation  = mean(deviation, na.rm = TRUE)
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

plot_dispersal <- sim_dispersal %>% 
  ggplot(aes(x=nloci, y=abs(deviation), colour=Inference, linetype=`Mean dispersal`)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Number of markers",
    y = "Error in median dispersal (m)"
  ) +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  facet_grid(family_size~prop_purged)

ggsave(
  filename = "05_manuscript/fig-sim_dipsersal.eps",
  plot = plot_dispersal,
  device = "eps",
  width = 16.9,
  height = 12,
  units = "cm"
)

