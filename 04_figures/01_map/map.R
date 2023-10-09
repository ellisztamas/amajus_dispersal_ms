# Script to draw the map of the sampled plants, highlighting maternal plants.

library('tidyverse')

gps <- read_csv("01_data/processed_GPS_positions.csv", col_types = 'cdd')

fathers = read_csv(
  "03_analysis/04_mating_events/output/01_mcmc_main/mating_events_over_chains.csv",
  col_types = 'iccdddddd'
  ) %>%
  filter(offspring >= 1, !is.na(father)) %>% 
  pull(father)

os <- 20
map_data <- gps %>%
  arrange(is_mother) %>% 
  mutate(
    sample_type = ifelse(PlantID_final %in% fathers, "Pollen donor", "Other plants"),
    sample_type = ifelse(is_mother, "Maternal plant", sample_type),
    sample_type = fct_relevel(sample_type, c("Maternal plant", "Pollen donor", "Other plants")),
    Easting = Easting - min(Easting),
    Northing = Northing - min(Northing),
    # Offset mothers and fathers a little
    Northing = ifelse(is_mother, Northing + os, Northing),
    Northing = ifelse(sample_type == "Pollen donor", Northing - os, Northing)
  )

map_plot <- map_data %>% 
  ggplot(aes(x = Easting, y = Northing, colour = sample_type)) +
  geom_point() +
  geom_point(
    data = map_data %>% filter(sample_type != "Other plants")
  ) +
  scale_colour_manual(values= c( "#E69F00", "#009E73", 'grey')) +
  coord_fixed() +
  labs(
    x = "Easting (m)",
    y = "Northing (m)"
  ) + 
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.title=element_blank()
  )

ggsave(
  filename = "05_manuscript/map.eps",
  plot = map_plot,
  device = "eps",
  width = 16.9,
  height = 10,
  units = "cm"
  )
