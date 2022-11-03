#' Script to predict the number of mating events detected as a function of how
#' many progeny were included in the maternal array using a Poisson GLM, using
#' progeny with an identifiable father only.
#' 
#' Second, the script uses that GLM to predict how many separate full sibships
#' there are among progeny with an unidentified (unsampled) father.
#' 
#' Tom Ellis, 27th October 2022.

library(tidyverse)

# Import data table listing mating events
me = read_csv(
  "03_analysis/04_mating_events/output/mating_events_over_chains.csv",
  show_col_types = FALSE
  # col_types = 'ccddddddffffff'
  ) %>% 
  left_join( 
    read_csv("01_data/maternal_family_sizes.csv", show_col_types = FALSE), by = "mother"
  ) %>% 
  rename(array_size = n_offspring)

# Average the probability of mating events having occurred over each iteration 
# of the MCMC.
mean_me <- me %>% 
  # Keep only mating events with a sampled father, and with at least two offspring
  # (offspring number is weighted by the probability of having occurred, so 
  # this may be less than 1.)
  filter(!is.na(father)) %>% 
  # Average over iterations
  group_by(mother, array_size) %>% 
  summarise(
    n = sum(prob) / 1000
  )

me %>% 
  mutate(father_found = !is.na(father)) %>% 
  group_by(father_found, iter) %>% 
  summarise(
    n = sum( offspring * prob)
  ) %>% 
  group_by(father_found) %>% 
  summarise(
    mean = mean(n),
    lower = quantile(n, 0.02),
    upper = quantile(n, 0.98)
    ) %>% 
  View()

me %>% 
  filter(is.na(father)) %>% 
  group_by(iter) %>% 
  summarise(
    n = n()
  ) %>%  pull(n)

ms <- read_csv("03_analysis/04_mating_events/output/01_mcmc_restrict_kurtosis/summarise_mating.csv")

mean(ms$n_mating_events)
quantile(ms$n_mating_events, c(0.02, 0.98))
  
  

# Poisson GLM for the number of mating events as a function of the number of 
# progeny in the maternal family.
mod <- glm(n ~ array_size, data = mean_me, family = poisson)
# Dataframe giving predicted number of mating events for maternal families of
# 1 to 30 from the GLM
pframe <- data.frame(
  array_size = 1:30
  ) %>% 
  mutate(
    n = predict(mod, newdata = ., type = 'response')
  )

# A table of groups of progeny in each maternal family for whom no father was
# found.
# Then add the expected number of full sibships expected within each group from
# the GLM above.
missing <- me %>%
  # Keep only mating events with an unsampled father, and with at least one offspring
  # (offspring number is weighted by the probability of having occurred, so 
  # this may be less than 1.)
  filter(is.na(father), offspring >= 1) %>%
  # Average over iterations
  group_by(mother) %>% 
  summarise(
    array_size = sum(offspring*prob) / 1000
  ) %>% 
  # Get predicted number of families from the GLM
  mutate(
    pred = predict(object = mod, newdata = ., type = "response"),
    n = 1,
    p = pred / array_size
  )

# Draw a binomially distributed number of full sibships for offspring with 
# unsampled father within each maternal family
sim <- sapply(1:1000, function(x){
  rbinom(n = missing$n, size = round(missing$array_size), prob = missing$p) 
  }) %>% 
  colSums()
# Mean and 96% credible intervals on the number of full sibships with a missing
# father
mean(sim)
quantile(sim, c(0.02, 0.98))

132 / (212+132)
