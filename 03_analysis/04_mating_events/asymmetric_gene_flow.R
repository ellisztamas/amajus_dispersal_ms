#' Check whether fathers are further East or West than mothers.
#' 
#' East-West distances between mother and father for each mating event, weighted
#' by posterior probability of each mating event having occurred.
#' Positive values indicate the father is further West of the mother.
#' 
#' Tom Ellis, 6th October 2023

library("tidyverse")
source("03_analysis/01_data_formatting/split_samples_by_road.R")

# Positions of plants on the lower road.
lower_road <- gps %>% 
  filter(lower_road) %>%
  select(PlantID_final, Easting, Northing)

# Import data table listing observed mating events
me = read_csv(
  "03_analysis/04_mating_events/output/01_mcmc_main/mating_events_over_chains.csv",
  col_types = 'iccdddddd') %>% 
  left_join( 
    read_csv("01_data/maternal_family_sizes.csv", col_types = 'ci'), by = "mother"
  ) %>% 
  rename(array_size = n_offspring)

# Import simulated mating events.
# Add columns for the East-West position, and the difference between mother and father
sim <-read_csv("03_analysis/05_simulations/simulated_mating_events.csv") %>% 
  filter(is.finite(offspring_real), q == 0.5) %>% 
  select(iter, mother, father) %>% 
  right_join(lower_road, by = c("mother" = 'PlantID_final')) %>% 
  right_join(lower_road, by = c("father" = 'PlantID_final'), suffix = c("_m", "_f")) %>% 
  mutate(diff = Easting_f - Easting_m)

# Histogram of relative positions of pollen donors to mothers in simulations.
sim %>% 
  group_by(iter) %>% 
  summarise(
    mean = mean(diff, na.rm = TRUE),
    goe = mean(diff > 0)
) %>% 
  ggplot(aes(x=mean)) +
  geom_histogram()
# Mean and 96% CIs of the relative positions in simulations
sim %>% 
  group_by(iter) %>% 
  summarise(
    means = mean(diff, na.rm = TRUE)
  ) %>% 
  summarise(
    mean=  mean(means, na.rm = TRUE),
    lower = quantile(means, 0.02, na.rm=TRUE),
    upper = quantile(means, 0.98, na.rm=TRUE)
  )

# Histogram of relative position of the father to mother in *real* mating events
me %>% 
  filter(!is.na(father), offspring >= 1) %>% 
  mutate(
    eastwest = Easting_father - Easting_mother
  ) %>%
  ggplot( aes(x = eastwest, weights = prob) ) +
  geom_histogram( aes(
    y = stat(count / sum(count))
  )) +
  labs(
    x = 'East-West distance from mother (m)',
    y = "Density"
  ) +
  theme_bw()
# Mean difference for real mating events
me %>% 
  filter(!is.na(father), offspring >= 1) %>% 
  mutate(
    eastwest = Easting_father - Easting_mother
  ) %>% 
  pull(eastwest) %>%  mean

# Count the number of mating events where the father was to the East or West of
# the mother.
me %>%
  filter(!is.na(father), offspring >= 1) %>%
  mutate(
    eastwest = Easting_father - Easting_mother
  ) %>%
  mutate( direction = ifelse(eastwest > 0, "East", "West")) %>%
  group_by(direction) %>%
  summarise(
    n = sum(prob) / max(me$iter),
  )

# Bootstrap resamples of observed mating distances.
# 97% of resamples have mean distance to the West.
east_west_distances <- me %>%
  filter(!is.na(father), offspring > 0) %>% 
  mutate(
    eastwest = Easting_father - Easting_mother
  ) 
set.seed(113)
perms <- sapply(unique(east_west_distances$iter), function(i) {
  d <- east_west_distances %>%
    filter(iter == i)
  
  mean(
    d[ sample(1:nrow(d), size=nrow(d), replace = TRUE), ]$eastwest
  )
})
mean(perms > 0)

hist(perms)
