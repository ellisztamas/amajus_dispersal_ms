library(tidyverse)

# sims <- read_csv("03_analysis/05_simulations/compare_paternity_accuracy.csv",
#                  col_types = 'fccnnnnc')
sims <- read_csv(
  "03_analysis/05_simulations/simulation_output.csv",
  col_types = 'fccnnnnn'
)

mean(sims$false_negative > sims$sibs_covs)

mean( exp(sims$true_negative) >0.5)





# sim_plot <- 
sims %>% 
  group_by(iter) %>% 
  summarise(
    Paternity = mean(exp(paternity) > 0.95),
    Sibships  = mean(exp(sibships) > 0.95),
    `Sibships +\ndispersal` = mean(exp(sibs_covs) > 0.95)
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