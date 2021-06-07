### Simulation histogram

### Open brazil 2008 map
bra_2008 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bra_dry_2008_102033.tif")

### Open brazil 2018 map
bra_2018 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bra_dry_2018_102033.tif")

### Open buffer
jaman <- st_read("/nfs/agfrontiers-data/Case Study Info/Brazil/JamanBuffer/JamanBuffer.shp") %>%
  st_transform(., crs = "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0
+ellps=aust_SA +units=m +no_defs")

### Crop LUC maps
jaman_08 <- crop(bra_2008, jaman)
jaman_08 <- mask(jaman_08, jaman)
jaman_18 <- crop(bra_2018, jaman)
jaman_18 <- mask(jaman_18, jaman)

### Reclassification matrix
reclass_df <- c(0, 1.5, 0.05,       ## ag
                1.6, 2.2, 0.1,      ## forest
                2.6, 3.2, 0.15,     ## bare soil
                3.6, 4.2, 0.2,      ## urban
                4.6, 5.2, 0.25,     ## wetland
                5.6, 7.2, 0.3)      ## water
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

### Reclassify 2018
bra_18_rc <- reclassify(jaman_18,
                        reclass_m)

### Reclassification for subtracted rasters
reclass_sub_df <- c(0, 1.94, 0,       ## other
                    1.945, 1.955, 1,  ## missed conversion to ag
                    1.96, 7.93, 0,    ## other
                    7.94, 7.96, 2,    ## stable ag from 2008
                    7.97, 9.88, 0,    ## other
                    9.89, 9.91, 3,    ## incorrectly predicted ag conversion
                    9.92, 9.94, 0,    ## other
                    9.945, 9.96, 4,   ## correctly predicted ag conversion
                    9.97, 15, 0)      
reclass_sub_m <- matrix(reclass_sub_df,
                        ncol = 3,
                        byrow = TRUE)

### Rasters for model 1
m1_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/brazil_1_project_2018",
                          pattern = "*.tif",
                          full.names = TRUE)

### Data frame
### do I need to make a blank data frame to fill here?

### Write loop
for (i in 1:length(m1_layers)) {
  
  ### Open first raster
  j_m1 <- raster(m1_layers[i])

  ### Restrict to extent of national forest
  j_m1 <- crop(j_m1, jaman)
  j_m1 <- mask(j_m1, jaman)
  

  ### Mask cells that were ag in 2008
  j1_mask_08 <- mask(j1, 
                     jaman_08,
                     inverse = FALSE,
                     maskvalue = 1,
                     updatevalue = 8)
  
  ### Subtract rasters
  j1_mask_08_sub <- j1_mask_08 - bra_18_rc
  
  ### Reclassify subtracted raster
  j1_mask_08_sub_rc <- reclassify(j1_mask_08_sub,
                                  reclass_sub_m)
  
  ### Record % of pixels with each cell value in dataframe (one row per raster)

}

### Save dataframe

