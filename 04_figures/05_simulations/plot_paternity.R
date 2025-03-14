library(tidyverse)
library(ggpubr)
# Import and format paternity accuracies

patsim <- read_csv(
  "03_analysis/05_simulations/output/simulate_fixed_family_sizes_paternity.csv",
  col_types = 'ffidifdddd'
) %>% 
  # Calculate rates
  mutate(
    n_pos = true_pos + false_pos,
    n_neg = true_neg + false_neg,
    p_neg = n_neg / (n_pos + n_neg),
    true_pos  = true_pos /  n_pos,
    false_pos = false_pos / n_pos,
    true_neg  = true_neg /  n_neg,
    false_neg = false_neg / n_neg
  ) %>% 
  # Average over replicates
  group_by(scale, nloci, prop_purged, family_size, data_type) %>% 
  summarise(
    true_pos  = mean(true_pos),
    false_pos = mean(false_pos),
    true_neg  = mean(true_neg),
    false_neg = mean(false_neg),
    p_neg     = mean(p_neg)
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

# Import and format mating results
sim_fixed <- read_csv(
  "03_analysis/05_simulations/output/simulate_fixed_family_sizes_mating.csv",
  col_types = 'ffdfffddddddddd')

sim_fixed <- sim_fixed %>% 
  group_by(scale, nloci, prop_purged, family_size, data_type) %>% 
  summarise(
    arith_prob_correct = mean(arith_prob_correct),
    weighted_prob_correct = mean(weighted_prob_correct),
    true_if_pp_is1 = mean(true_if_pp_is1),
    true_if_pp_lt1 = mean(true_if_pp_lt1),
    true_if_pp_is99 = mean(true_if_pp_is99),
    true_if_pp_lt99 = mean(true_if_pp_lt99),
    true_if_offs_ge1 = mean(true_if_offs_ge1),
    true_if_offs_lt1 = mean(true_if_offs_lt1),
    missing = mean(missing)
  ) 

plot_correct_paternity <- patsim %>% 
  ggplot(aes(x=nloci, y=true_pos, colour=Inference, linetype=`Mean dispersal`)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Number of markers",
    y = "Probability correct paternity"
  ) +
  lims(
    y = c(0.4,1)
  ) +
  theme_bw() +
  facet_grid(family_size~prop_purged)

plot_assigned_missing <- patsim %>%
  ggplot(aes(x=nloci, y=p_neg, colour=Inference, linetype=`Mean dispersal`)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Number of markers",
    y = "Offspring assigned to missing father"
  ) +
  lims(
    y = c(0,1)
  ) +
  theme_bw() +
  facet_grid(family_size~prop_purged)

ggarrange(
  plot_correct_paternity, plot_assigned_missing,
  nrow=2, ncol=1,
  labels="AUTO",
  common.legend = TRUE, legend = 'right'
)


ggsave(
  filename = "05_manuscript/fig-sim_paternity.eps",
  device = "eps",
  width = 16.9,
  height = 20,
  units = "cm"
)
