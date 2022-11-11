library(tidyverse)

sims <- read_csv("03_analysis/05_simulations/simulation_output.csv")

sim_plot <- sims %>% 
  group_by(iter) %>% 
  summarise(
    Paternity = mean(exp(paternity) > 0.99),
    Sibships  = mean(exp(sibships) > 0.99),
    `Sibships +\ndispersal` = mean(exp(sibs_covs) > 0.99)
  ) %>% 
  pivot_longer(Paternity: `Sibships +\ndispersal` ) %>% 
  mutate(
    name=fct_relevel(name,c("Paternity","Sibships","Sibships +\ndispersal"))
    ) %>% 
  ggplot(aes(x = name, y = value))+
  geom_boxplot() +
  labs(
    x = 'Information',
    y = "Proportion of simulations"
  ) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "05_manuscript/simulations.eps",
  plot = sim_plot,
  width = 80, 
  height = 100,
  units = "mm",
  device = "eps"
)