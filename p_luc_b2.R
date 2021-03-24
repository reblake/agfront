### Code to make LUC models and apply to 2018 LUC maps

### Open packages
library(tidyverse)
library(sf)
library(raster)
library(spatialEco)

####################################################################################
####################################################################################

### Loop through simulated maps and: 
### 1. make raster that adds values of all map layers
### 2. calculate forest loss within Jamanxim NF and put that value in a vector
### 3. calculate FRAGSTATS and put in dataframe -- can't allocate vector...

# ### Add library for FRAGSTATS
# library(landscapemetrics)

### Open Jamanxim without buffer
# jaman <- st_read("/nfs/agfrontiers-data/Case Study Info/Brazil/JamanBuffer/jamanxim_wdpa_march21_102033.shp")

###################
########## Model 1
###################
# ### Try on subset of layers in Brazil Model 1
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/brazil_1_project", 
#                           pattern = "*.tif", 
#                           full.names = TRUE)
# 
# ### First set
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# # test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Vector for area converted from forest to ag within Jaman
# jam_area_vec <- rep("", times = 200)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
# for (i in 1:length(test_layers_5)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_5[i])
#   
#   ###################################
#   ### Add raster values
#   ###################################
#   ### Stack with raster_fill
#   raster_fill <- stack(raster_fill, ch_rast)
#   
#   ### Add raster values
#   raster_fill <- calc(raster_fill, sum)
#   
#   ###################################
#   ### Calculate forest loss within JNP
#   ###################################  
#   ### Clip raster to shapefile
#   rast_crop <- crop(ch_rast, jaman)
#   rast_crop <- mask(rast_crop, jaman)
#   
#   ### Calculate # pixels that convert
#   forest_loss <- freq(rast_crop, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss <- forest_loss * res(rast_crop)[1] * res(rast_crop)[2]
#   
#   ### Add to vector
#   jam_area_vec[i] <- forest_loss
#   
# }
# 
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_m1_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Save vector of area lost
# write.csv(jam_area_vec, 
#           file = "/nfs/agfrontiers-data/luc_model/brazil_m1_jamanxim_forestloss_5.csv", 
#           row.names=FALSE)

###################
########## Model 2
###################
# ### Try on subset of layers in Brazil Model 2
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/brazil_2_project", 
#                           pattern = "*.tif", 
#                           full.names = TRUE)
# 
# ### First set
# # test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# # test_layers_4 <- test_layers[601:800]
# # test_layers_5 <- test_layers[801:1000]
# 
# ### Vector for area converted from forest to ag within Jaman
# jam_area_vec <- rep("", times = 200)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
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
#   
#   ###################################
#   ### Calculate forest loss within JNP
#   ###################################  
#   ### Clip raster to shapefile
#   rast_crop <- crop(ch_rast, jaman)
#   rast_crop <- mask(rast_crop, jaman)
#   
#   ### Calculate # pixels that convert
#   forest_loss <- freq(rast_crop, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss <- forest_loss * res(rast_crop)[1] * res(rast_crop)[2]
#   
#   ### Add to vector
#   jam_area_vec[i] <- forest_loss
#   
# }
# 
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_m2_stacked_2.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Save vector of area lost
# write.csv(jam_area_vec, 
#           file = "/nfs/agfrontiers-data/luc_model/brazil_m2_jamanxim_forestloss_2.csv", 
#           row.names=FALSE)

###################
########## Model 3
###################
# ### Try on subset of layers in Brazil Model 3
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/brazil_3_project", 
#                           pattern = "*.tif", 
#                           full.names = TRUE)
# 
# ### First set
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# # test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Vector for area converted from forest to ag within Jaman
# jam_area_vec <- rep("", times = 200)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
# for (i in 1:length(test_layers_5)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_5[i])
#   
#   ###################################
#   ### Add raster values
#   ###################################
#   ### Stack with raster_fill
#   raster_fill <- stack(raster_fill, ch_rast)
#   
#   ### Add raster values
#   raster_fill <- calc(raster_fill, sum)
#   
#   ###################################
#   ### Calculate forest loss within JNP
#   ###################################  
#   ### Clip raster to shapefile
#   rast_crop <- crop(ch_rast, jaman)
#   rast_crop <- mask(rast_crop, jaman)
#   
#   ### Calculate # pixels that convert
#   forest_loss <- freq(rast_crop, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss <- forest_loss * res(rast_crop)[1] * res(rast_crop)[2]
#   
#   ### Add to vector
#   jam_area_vec[i] <- forest_loss
#   
# }
# 
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_m3_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Save vector of area lost
# write.csv(jam_area_vec, 
#           file = "/nfs/agfrontiers-data/luc_model/brazil_m3_jamanxim_forestloss_5.csv", 
#           row.names=FALSE)

###################
########## Model 4
###################
### Try on subset of layers in Brazil Model 4
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/brazil_4_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# ### First set
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# # test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Vector for area converted from forest to ag within Jaman
# jam_area_vec <- rep("", times = 200)
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
# for (i in 1:length(test_layers_5)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_5[i])
#   
#   ###################################
#   ### Add raster values
#   ###################################
#   ### Stack with raster_fill
#   raster_fill <- stack(raster_fill, ch_rast)
#   
#   ### Add raster values
#   raster_fill <- calc(raster_fill, sum)
#   
#   ###################################
#   ### Calculate forest loss within JNP
#   ###################################
#   ### Clip raster to shapefile
#   rast_crop <- crop(ch_rast, jaman)
#   rast_crop <- mask(rast_crop, jaman)
#   
#   ### Calculate # pixels that convert
#   forest_loss <- freq(rast_crop, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss <- forest_loss * res(rast_crop)[1] * res(rast_crop)[2]
#   
#   ### Add to vector
#   jam_area_vec[i] <- forest_loss
#   
# }
# 
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_m4_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Save vector of area lost
# write.csv(jam_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/brazil_m4_jamanxim_forestloss_5.csv",
#           row.names=FALSE)

### Peru Models
### model 1
# ### Try on subset of layers in Peru Model 1
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_1_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
# # test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
# for (i in 1:length(test_layers_4)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_4[i])
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
#             filename = "/nfs/agfrontiers-data/luc_model/peru_m1_stacked_4.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

# #### model 2
# ### Try on subset of layers in Peru Model 2
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_2_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
# # test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
# for (i in 1:length(test_layers_4)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_4[i])
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
#             filename = "/nfs/agfrontiers-data/luc_model/peru_m2_stacked_4.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###################
########## Model 3
###################
# ### Try on subset of layers in Peru Model 3
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_3_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# # test_layers_3 <- test_layers[401:600]
# # test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
# for (i in 1:length(test_layers_5)) {
#   
#   ### Open raster
#   ch_rast <- raster(test_layers_5[i])
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
#             filename = "/nfs/agfrontiers-data/luc_model/peru_m3_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###################
########## Model 4
###################
### Try on subset of layers in Peru Model 4
test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_4_project",
                          pattern = "*.tif",
                          full.names = TRUE)

# test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
test_layers_5 <- test_layers[801:1000]

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop
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
            filename = "/nfs/agfrontiers-data/luc_model/peru_m4_stacked_5.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))