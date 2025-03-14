#' Script to calculte and plot posterior means and 96% credible intervals for 
#' MCMC parameters

library(tidyverse)
library(ggpubr)

source("02_library/R/import_MCMC_results.R")

#x Import all the MCMC results and concatenate into a single data.frame
mcmc_results <- import_mcmc_results("03_analysis/03_mcmc/01_mcmc_main") %>% 
  filter(iter > 500) # Remove the first 500 iterations as burnin

scale <- mcmc_results %>% 
  ggplot(aes(x = scale, colour = chain, groups=chain, fill = chain)) +
    geom_histogram() +
    labs(
      x = "Scale",
      y = "Posterior density"
    ) +
    theme_bw() +
    theme(
      strip.text.y = element_blank()
    )

shape <- mcmc_results %>% 
  ggplot(aes(x = shape, colour = chain, groups=chain, fill = chain)) +
  geom_histogram() +
  labs(
    x = "Shape",
    y = "Posterior density"
  ) +
  theme_bw() +
  theme(
    strip.text.y = element_blank()
  )

mixture <- mcmc_results %>% 
  ggplot(aes(x = mixture, colour = chain, groups=chain, fill = chain)) +
  geom_histogram() +
  labs(
    x = "Mixture parameter",
    y = "Posterior density"
  ) +
  theme_bw() +
  theme(
    strip.text.y = element_blank()
  )


plot_posteriors <- ggarrange(
  shape, scale, mixture,
  ncol=3,
  common.legend = TRUE, legend = 'bottom',
  labels = "AUTO"
  )

ggsave(
  filename = "05_manuscript/fig-posterior_distributions.eps",
  plot = plot_posteriors,
  device = "eps",
  width = 16.9,
  height = 8,
  units = "cm"
)
