# Add a column to the GPS data to indicate whether samples are from the upper
# or lower roads.

# Tom Ellis, 6th October 2023

library(tidyverse)

gps <- read_csv("01_data/processed_GPS_positions.csv")

# Find three transects that divide the two roads
sv <- rbind(
  c(-3000, -500),
  c(-1700,  25),
  c(-100, -330),
  c(2000, -330)
)
plot(gps$Easting, gps$Northing)
points(sv[,1], sv[,2], type = 'l')

# Treat each transect a a regression line
mod <- list(
  lm(sv[1:2,2] ~ sv[1:2,1]),
  lm(sv[2:3,2] ~ sv[2:3,1]),
  lm(sv[3:4,2] ~ sv[3:4,1])
)

# Coordinates for plants in each transect along the East-West axis.
ix <- list(
  (gps$Easting >= -3000) & (gps$Easting < -1700),
  (gps$Easting >= -1700) & (gps$Easting < -100),
  (gps$Easting >= -100) & (gps$Easting < 2000)
)

# Create a column indicating whether a sample is from the lower road or not.
gps$lower_road <- NA
for(i in 1:length(ix)){
  a <- mod[[i]]$coefficients[1]
  b <- mod[[i]]$coefficients[2]
  x <- gps$Easting[ix[[i]]]
  
  pred <- a+b*x
  gps$lower_road[ix[[i]]] <- ifelse(gps$Northing[ix[[i]]] < pred, TRUE, FALSE)
}
#Plot to make sure it works.
plot(gps$Easting, gps$Northing, col=gps$lower_road+1)

write_csv(gps, "01_data/GPS_with_road.csv")

