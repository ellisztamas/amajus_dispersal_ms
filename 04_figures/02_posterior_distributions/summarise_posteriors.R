#' Script to prepare MCMC output for plotting summaries of means and CIs.
#' 
#' This goes through folders of results from multiple MCMC analyses to pull out
#' the full chain and summary files, and joins these. It then calculates means
#' and 96% credible intervals for parameters of interest ready to be plotted.

library(tidyverse)
library(ggpubr)

source("02_library/R/import_MCMC_results.R")

# Import all the MCMC results and concatenate into a single data.frame
mcmc_results <- list.dirs("03_analysis/03_mcmc", recursive = FALSE) %>%
  map(import_mcmc_results) %>%
  do.call(what = "rbind") %>%
  filter(iter > 500) %>% # Remove the first 500 iterations as burnin
  select(mixture, scale, shape, scenario) %>% 
  pivot_longer(mixture : shape, names_to = "parameter")

# Import summaries of mating events
sm <- Sys.glob("03_analysis/04_mating_events/output/*/summarise_mating.csv")
mating_summary <- lapply(sm, function(dir) {
  read_csv(dir, col_types = 'dddddd') %>% 
    mutate(scenario = strsplit(dir, split = "/")[[1]][[4]])
})
mating_summary <- do.call('rbind', mating_summary) %>% 
  select(scenario, n_mating_events, orphans, median_dispersal, `dispersal>500m`) %>% 
  rename(long_range_dispersal = `dispersal>500m`) %>% 
  pivot_longer(n_mating_events : long_range_dispersal, names_to = "parameter")

# Bind the two tables of results, and summarise means and CIs for each parameter
mcmc_summary <- mcmc_results %>% 
  rbind(mating_summary) %>% 
  group_by(scenario, parameter) %>% 
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.02),
    upper = quantile(value, 0.98)
  ) %>%  
  # Change the names for each analysis to something human readable
  mutate(
    scenario = recode_factor(
      scenario,
      `01_mcmc_restrict_kurtosis` = "Restricted\nkurtosis",
      `02_mcmc_short_range_kurtosis` = "Short-range\ndispersal",
      `03_mcmc_strong_prior_on_kurtosis` = "Penalise\nkurtosis",
      `04_mcmc_022_missing_fathers` = "Fewer missing\nfathers",
      `05_mcmc_042_missing_fathers` = "More missing\nfathers",
      `06_no_mixture` = "No mixture")
  )

# Function to plot means and credible intervals from the object created above.
plot_parameter <- function(p, ylab){
  mcmc_summary %>% 
    filter(parameter == p) %>% 
    ggplot(aes( y = mean, x = scenario)) +
    geom_point() +
    geom_errorbar(aes( ymin=lower, ymax=upper), width=0.1 ) +
    labs(
      x = element_blank(),
      y = ylab
    ) +
    theme_bw()
}
