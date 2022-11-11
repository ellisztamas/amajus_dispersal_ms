source("04_figures/02_posterior_distributions/summarise_posteriors.R")

plist <- list(
  plot_parameter("scale", "Scale") + theme( axis.text.x = element_blank() ),
  plot_parameter("shape", "Shape") + theme( axis.text.x = element_blank() ),
  plot_parameter("mixture", "Mixture parameter")  + theme( axis.text.x = element_blank() ),
  # plot_parameter("orphans", "Offsping w. unsampled father") + theme( axis.text.x = element_blank() ),
  plot_parameter("n_mating_events", "N. mating events") + 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
      )
)

plot_posteriors <- ggarrange(plotlist = plist, ncol=1, heights = c(1,1,1,1.5), labels = 'AUTO')

ggsave(
  filename = "05_manuscript/posterior_distributions.eps",
  plot = plot_posteriors,
  device = "eps",
  width = 16.9,
  height = 22,
  units = "cm"
)
