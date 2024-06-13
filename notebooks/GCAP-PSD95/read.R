library(spatstat)
library(lgcp)
library(tidyverse)


original_psd95 = read.csv('Wouter data/LTP_Ch1_PSD95.csv')
original_GKAP = read.csv('Wouter data/LTP_Ch2_GKAP.csv')

Cell_ID = 3
PSD_ID = 1

subset_psd95 <- original_psd95 %>%
  filter(Cell == Cell_ID, PSDnumber == PSD_ID)

subset_gkap <- original_GKAP %>%
  filter(Cell == Cell_ID, PSDnumber == PSD_ID)

x_coordinates_subset <- subset_psd95$x
y_coordinates_subset <- subset_psd95$y
W <- owin(c(min(x_coordinates_subset), max(x_coordinates_subset)), c(min(y_coordinates_subset), max(y_coordinates_subset)))
ppPSD <- ppp(x_coordinates_subset, y_coordinates_subset, W)

x_coordinates_subset <- subset_gkap$x
y_coordinates_subset <- subset_gkap$y
W <- owin(c(min(x_coordinates_subset), max(x_coordinates_subset)), c(min(y_coordinates_subset), max(y_coordinates_subset)))
ppGKAP <- ppp(x_coordinates_subset, y_coordinates_subset, W)


# Plot an empty window as the base layer
plot(W, main="Combined point patterns")

# Add points from the first ppp object
points(ppGKAP, pch=19, col="red",cex = 0.5)

# Add points from the second ppp object
points(ppPSD, pch=19, col="blue",cex = 0.5)


