# Script to draw the map of the sampled plants, highlighting maternal plants.

library('tidyverse')

gps <- read_csv("01_data/processed_GPS_positions.csv", col_types = 'cdd')


map <- gps %>% 
  arrange(is_mother) %>% 
  mutate(
    is_mother = ifelse(is_mother, "Maternal plant", "Candidate father"),
    Easting = Easting - min(Easting),
    Northing = Northing - min(Northing)
  ) %>% 
  ggplot(aes(x = Easting, y = Northing, colour = is_mother)) +
  geom_point() +
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

# setEPS()
# postscript(
#   file = "04_figures/01_map/map.eps",
#   width = 16.9/2.54,
#   height = 10/2.54
# )
png(
  filename = "04_figures/01_map/map.png",
  units = "cm",
  width = 16.9,
  height = 10,
  res = 300
)

map

dev.off()
