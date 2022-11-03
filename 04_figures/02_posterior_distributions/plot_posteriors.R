source("04_figures/02_posterior_distributions/summarise_posteriors.R")

plist <- list(
  plot_parameter("scale", "Scale") + theme( axis.text.x = element_blank() ),
  plot_parameter("shape", "Shape") + theme( axis.text.x = element_blank() ),
  plot_parameter("mixture", "Mixture parameter")  + theme( axis.text.x = element_blank() ),
  plot_parameter("n_mating_events", "N. mating events") + theme( axis.text.x = element_blank() ) + 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
      )
)

png(
  filename = "04_figures/02_posterior_distributions/posterior_distributions.png",
  height = 22, width = 8, units = 'cm',
  res = 600
)

ggarrange(plotlist = plist, ncol=1, heights = c(1,1,1,1.3), labels = 'AUTO')

dev.off()
