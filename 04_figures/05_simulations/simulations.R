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


# sims %>% 
#   mutate(
#     paternity = exp(paternity),
#     sibships  = exp(sibships),
#     sibs_covs = exp(sibs_covs)
#   )
# 
# 
# 
# sims %>% 
#   pivot_longer(paternity:sibs_covs) %>% 
#   mutate(
#     value = exp(value),
#     bin = cut(value,  breaks =c(0,0.1, 0.9, 0.95, 1, 1.1), right=FALSE)) %>% 
#   group_by(iter, name, bin) %>% 
#   summarise(
#     n = n()
#     ) %>%
#   # group_by(iter, name) %>% 
#   # summarise(
#   #   n = sum(n)
#   # ) %>%  pull(n) %>%  table
#   group_by(name, bin) %>% 
#   summarise(
#     mean = mean(n),
#     n = n()
#   ) %>%
#   ggplot( aes(x = bin, y = mean, group = name, colour = name) ) +
#   geom_point() +
#   geom_line()
# 
#          