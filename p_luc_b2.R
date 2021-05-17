### Code to make LUC models and apply to 2018 LUC maps

### Open packages
library(tidyverse)
library(sf)
library(raster)
# library(spatialEco)

####################################################################################
####################################################################################

### Loop through simulated maps and make raster that adds values of all map layers
### This code is for the Bolivia case study, Models 3 & 4
### Making 5 stacks per model because the code didn't work when I tried to stack all 1000 rasters at once.


###################
########## Model 3
###################
### Bolivia Model 3 layer subset
test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/bolivia_3_project",
                          pattern = "*.tif",
                          full.names = TRUE)

test_layers_1 <- test_layers[1:200]
test_layers_2 <- test_layers[201:400]
test_layers_3 <- test_layers[401:600]
test_layers_4 <- test_layers[601:800]
test_layers_5 <- test_layers[801:1000]

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop for stack 1
for (i in 1:length(test_layers_1)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_1[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}

### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_1.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop for stack 2
for (i in 1:length(test_layers_2)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_2[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}
### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_2.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop for stack 3
for (i in 1:length(test_layers_3)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_3[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}
### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_3.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop for stack 4
for (i in 1:length(test_layers_4)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_4[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}
### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_4.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop for stack 5
for (i in 1:length(test_layers_5)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_5[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}
### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_5.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)


###################
########## Model 4-- batches 1-3 already ran
###################
test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/bolivia_4_project",
                          pattern = "*.tif",
                          full.names = TRUE)

# test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
test_layers_4 <- test_layers[601:800]
test_layers_5 <- test_layers[801:1000]

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

# ### Write loop for stack 1
# for (i in 1:length(test_layers_1)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_1[i])
#   
#   ###################################
#   ### Add raster values
#   ###################################
#   ### Stack with raster_fill
#   raster_fill <- stack(raster_fill, ch_rast)
#   
#   ### Add raster values
#   raster_fill <- calc(raster_fill, sum)
# }
# 
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_1.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# rm(raster_fill)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop for stack 2
# for (i in 1:length(test_layers_2)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_2[i])
#   
#   ###################################
#   ### Add raster values
#   ###################################
#   ### Stack with raster_fill
#   raster_fill <- stack(raster_fill, ch_rast)
#   
#   ### Add raster values
#   raster_fill <- calc(raster_fill, sum)
# }
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_2.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# rm(raster_fill)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop for stack 3
# for (i in 1:length(test_layers_3)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_3[i])
#   
#   ###################################
#   ### Add raster values
#   ###################################
#   ### Stack with raster_fill
#   raster_fill <- stack(raster_fill, ch_rast)
#   
#   ### Add raster values
#   raster_fill <- calc(raster_fill, sum)
# }
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_3.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# rm(raster_fill)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0

### Write loop for stack 4
for (i in 1:length(test_layers_4)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_4[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}
### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_4.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop for stack 5
for (i in 1:length(test_layers_5)) {
  
  ### Open raster
  ch_rast <- raster(test_layers_5[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
}
### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_5.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
rm(raster_fill)