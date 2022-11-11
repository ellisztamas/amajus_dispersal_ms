library(tidyverse)
library(ggpubr)

files <- Sys.glob("03_analysis/04_mating_events/output/*/missing_sires.csv")

# Import mating events from multiple scenarios
mating_summary <- lapply(files, function(dir) {
  read_csv(dir, col_types = 'dddddd') %>%
    mutate(scenario = strsplit(dir, split = "/")[[1]][[4]])
}) %>% 
  do.call(what = "rbind") %>% 
  # Calculate a proportion of missing sires
  mutate(
    sim_missing_sires = sim_missing_sires / (sim_missing_sires + n_mating_events)
  ) %>% 
  # Change the names for each analysis to something human readable
  mutate(
    scenario = recode_factor(
      scenario,
      `01_mcmc_restrict_kurtosis` = "Restricted kurtosis",
      `02_mcmc_short_range_kurtosis` = "Short-range dispersal",
      `03_mcmc_strong_prior_on_kurtosis` = "Penalise kurtosis",
      `04_mcmc_022_missing_fathers` = "Fewer missing fathers",
      `05_mcmc_042_missing_fathers` = "More missing fathers",
      `06_no_mixture` = "No mixture")
  ) %>% 
  # Get means and 96% credible intervals
  pivot_longer(n_mating_events: sim_coef, names_to = "parameter") %>% 
  group_by(parameter, scenario) %>% 
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.02),
    upper = quantile(value, 0.98),
    .groups = "drop_last"
  )

plot_n_mating_events <- mating_summary %>% 
  filter(parameter == "n_mating_events") %>% 
  ggplot(aes( y = mean, x = scenario)) +
  geom_point() +
  geom_errorbar(aes( ymin=lower, ymax=upper), width=0.1 ) +
  labs(
    x = element_blank(),
    y = "Number of mating events"
  ) +
  theme_bw() +
  theme( axis.text.x = element_blank() )

plot_missing_sires <- mating_summary %>% 
  filter(parameter == "sim_missing_sires") %>% 
  ggplot(aes( y = mean, x = scenario)) +
  geom_point() +
  geom_errorbar(aes( ymin=lower, ymax=upper), width=0.1 ) +
  labs(
    x = element_blank(),
    y = "Missing fathers"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
  )

plot_mating_events <- ggarrange(plot_n_mating_events, plot_missing_sires,
          nrow = 2, heights = c(1,1.1), labels = "AUTO")

ggsave(
  filename = "05_manuscript/mating_events.eps",
  device = "eps",
  width = 8, height = 15,
  plot = plot_mating_events
)


# source("03_analysis/04_mating_events/n_mating_events.R")

plot_parameter("n_mating_events", "N. mating events") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
  )

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
