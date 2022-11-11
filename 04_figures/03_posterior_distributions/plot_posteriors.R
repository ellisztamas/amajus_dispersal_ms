#' Script to calculte and plot posterior means and 96% credible intervals for 
#' MCMC parameters

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


# Summarise means and CIs for each parameter
mcmc_summary <- mcmc_results %>% 
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


plist <- list(
  plot_parameter("shape", "Shape") + theme( axis.text.x = element_blank() ),
  plot_parameter("mixture", "Mixture parameter")  + theme( axis.text.x = element_blank() ),
  plot_parameter("scale", "Scale") + 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
      )
)

plot_posteriors <- ggarrange(plotlist = plist, ncol=1, heights = c(1,1,1.5), labels = 'AUTO')

ggsave(
  filename = "05_manuscript/posterior_distributions.eps",
  plot = plot_posteriors,
  device = "eps",
  width = 8,
  height = 22,
  units = "cm"
)
