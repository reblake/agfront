### Script to get neighboring pixel info for Bolivia
# Optimized version by QDR, 19 Feb 2021

### Load libraries
library(raster)

### Open data
# boliv <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bol_dry_2008_102033.tif")
# brazi <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bra_dry_2008_102033.tif")
# peru <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_per_dry_2008_102033.tif")
# peru18 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_per_dry_2018_102033.tif")
# boliv18 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bol_dry_2018_102033.tif")
brazi18 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bra_dry_2018_102033.tif")

### Define function taking the raster as input.
# All mapping across cells is internal in the function so the function only needs to be run once.
# There is no parallelization anymore because the whole thing runs in a few minutes!
percent_different_neighbors <- function(r) {
  
  cells <- 1:ncell(r) # Get the number of cells in the raster.
  col_index <- colFromCell(r, cells) # Get the column index of each cell in the raster.
  r_values <- getValues(r) # Get the values of the raster.
  nc <- ncol(r) # Get the number of columns in the raster.
  
  # For each cell, the same adjustment is used to find the index of its neighbors.
  nbrs <- c(0, -nc-1, -nc, -nc+1, -1, +1, +nc-1, +nc, +nc+1) 
  
  ad <- sapply(nbrs, function(x) cells + x) # Add that adjustment value to every cell index.
  colnames(ad) <- c("center", "NW", "N", "NE", "W", "E", "SW", "S", "SE") # Just FYI to show you which neighbor each column represents.
  
  # these three lines deal with edge issues, replacing the appropriate neighbor index values for cells that are on the edge of the raster with NA
  ad[col_index == 1, c(2, 5, 7)] <- NA
  ad[col_index == nc, c(4, 6, 9)] <- NA
  ad[ad<1] <- NA
  
  # For each cell, get the values corresponding to the indices of its neighbors.
  advalues <- t(apply(ad, 1, function(idx) r_values[idx]))
  
  # For each of the sets of neighbor values, calculate the percent different from the focal cell and return the result.
  apply(advalues, 1, function(x) {
    x <- x[!is.na(x)]
    sum(x!=x[1])/(length(x)-1)*100
  })
  
}

### Call the function on the boliv raster
# perc_diff_values <- percent_different_neighbors(boliv)
perc_diff_values <- percent_different_neighbors(brazi18)

### Make a new raster with the cells
# perc_diff_raster <- setValues(boliv, perc_diff_values)
perc_diff_raster <- setValues(brazi18, perc_diff_values)

# writeRaster(perc_diff_raster, 
#             filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/bol_perc_diff.tiff",
#             format = "GTiff",
#             overwrite = TRUE)

# writeRaster(perc_diff_raster, 
#             filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/bra_perc_diff.tiff",
#             format = "GTiff",
#             overwrite = TRUE)

# writeRaster(perc_diff_raster, 
#             filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/per_perc_diff.tiff",
#             format = "GTiff",
#             overwrite = TRUE)

# writeRaster(perc_diff_raster, 
#             filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/per18_perc_diff.tiff",
#             format = "GTiff",
#             overwrite = TRUE)

# writeRaster(perc_diff_raster,
#             filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/bol18_perc_diff.tiff",
#             format = "GTiff",
#             overwrite = TRUE)

writeRaster(perc_diff_raster,
            filename = "/nfs/agfrontiers-data/Remote Sensing/KS files/bra18_perc_diff.tiff",
            format = "GTiff",
            overwrite = TRUE)