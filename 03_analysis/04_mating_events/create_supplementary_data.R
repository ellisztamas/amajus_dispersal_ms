#' Save a table of mating events as a supplementary data table
#' 
#' Takes all ~220 mating events over the MCMC chains and averages the posterior
#' probability and offspring number over the chains.

library('tidyverse')

mating_events <- read_csv(
  "03_analysis/04_mating_events/output/mating_events_over_chains.csv"
  )

max_i <- max(mating_events$iter)

mating_events %>% 
  group_by(mother, father, distance) %>% 
  summarise(
    prob = sum(prob) / max_i,
    offspring = sum(offspring) / max_i
  ) %>% 
  write_csv(
    "05_manuscript/mating_events.csv"
  )


mating_events %>% 
  group_by(mother, father, distance) %>% 
  summarise(
    prob = sum(prob) / max_i,
    offspring = sum(offspring) / max_i
  ) %>% 
  filter(!is.na(father))

mating_events %>% 
  filter(!is.na(father)) %>% 
  group_by(iter) %>% 
  summarise(
    n = n()
  ) %>% 
  pull(n) %>%  quantile(0.96)
  
mating_events %>% 
  filter(! is.na(father) , prob >= 1) %>% 
  group_by(iter) %>% 
  summarise(
    n = n()
  ) %>% 
  pull(n) %>%  
  mean()
  quantile(0.02)
  
  
mating_events %>% 
  filter(is.na(father)) %>% 
  group_by(iter) %>% 
  summarise(
    mean = sum(offspring)
  ) %>% 
  pull(mean) %>% quantile(c(0.02, 0.96))
