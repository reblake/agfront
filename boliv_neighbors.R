### Script to get neighboring pixel info for Bolivia

### Load libraries
library(tidyverse)
library(raster)

### Open data
boliv <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bol_dry_2008_102033.tif")

### Function to calculate % neighbors of different types
percent_different_neighbors <- function(cell) {
  # pull neighboring values
  # see ?adjacent for the arguments used here
  x <- getValues(boliv)[adjacent(boliv, 
                                 cells = cell, 
                                 directions = 8, 
                                 pairs = F, 
                                 include = T)]
  # calculate the percent of them that are not equal to the focal cell (i.e., not equal to the first value of x)
  perc_diff <- sum(x!=x[1])/(length(x)-1)*100
  return(perc_diff)
}

### Run on sample matrix: apply function to each cell
perc_diff_values <- map_dbl(1:ncell(boliv),
                            ~percent_different_neighbors(.))

### Make a new raster with the cells
perc_diff_raster <- boliv %>% setValues(perc_diff_values)

writeRaster(perc_diff_raster, 
            filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/bol_perc_diff.tiff",
            format = "GTiff",
            overwrite = TRUE)
