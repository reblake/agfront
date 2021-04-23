### Code to make LUC models and apply to 2018 LUC maps

### Open packages
library(tidyverse)
library(sf)
library(raster)
library(spatialEco)

##################
### peru model 1 
##################

### Data prep
#############
# 
# ### Open sample points
# per_dat_all <- read.csv("/nfs/agfrontiers-data/luc_model/peru_data_for_model.csv")
# 
# 
# ### Make column for transition/none
# per_dat_all <- per_dat_all %>%
#   mutate(transition = ifelse(lc_2018 == "2",
#                              "no_change",
#                              ifelse(lc_2018 == "1",
#                                     "f_to_ag",
#                                     "other_change")),
#          trans_rc = ifelse(transition == "no_change",
#                            "0", ifelse(transition == "f_to_ag",
#                                        "1", "2")))
# 
# ### Make columns factor
# fact_cols <- c("uid", "prot_status",
#                "minconc",
#                "nonforconc", "npchg",
#                "permprod", "protforest",
#                "refor",
#                "lc_2018", "trans_rc")
# per_dat_all <- per_dat_all %>%
#   mutate_each_(funs(factor(.)), fact_cols)
# 
# ### Restrict to transitions
# per_dat_use <- per_dat_all %>%
#   filter(trans_rc == "1" | trans_rc == "0")

### Model 1
###########

# ### Run regression
# fit_per <- glm(trans_rc ~ aspect + cities + 
#                  slope + 
#                  crop_suit + prot_status + popdens + 
#                  precip + rivers + dist_ag + perc_diff,
#                data = per_dat_use,
#                family = binomial(link = "logit"))
# 
# ### 2018 data
# #############
# 
# ### Open rasterbrick
# per_layers_all <- brick("/nfs/agfrontiers-data/luc_model/peru_all_layers.tif")
# 
# ### Rename columns
# names(per_layers_all) <- c("aspect", "cities",
#                            "crop_suit", "elev",
#                            "mines", "prot_status",
#                            "popdens", "poverty",
#                            "precip", "rivers",
#                            "roads", "slope",
#                            "soil", "dist_ag",
#                            "control", "commun",
#                            "fire", "forest",
#                            "illegal_mines", "minconc",
#                            "nonforconc", "npchg",
#                            "permprod", "protforest",
#                            "refor", "tourism")
# 
# ### Replace precip layer
# per_layers_all <- dropLayer(per_layers_all, 9)
# 
# ### Open new precip layer
# per_precip <- raster("/nfs/agfrontiers-data/luc_model/per_precip_new.tif")
# 
# ### Add new precip layer
# per_layers_all <- addLayer(per_layers_all, per_precip)
# names(per_layers_all)[26] <- "precip"
# 
# ### Open 2018 percent of neighboring pixels raster
# neigh_2018 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/per18_perc_diff.tif")
# 
# ### Add to other layers
# per_18_all <- addLayer(per_layers_all, neigh_2018)
# names(per_18_all)[27] <- "perc_diff"


# ### Open distance to ag in 2018 layer
# d_ag_18 <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Peru/Processed Data/peru_dist_ag_18.tif")
# 
# ### Drop 2008 distance to ag, replace with 2018 distance to ag
# per_18_all <- addLayer(per_18_all, d_ag_18)
# per_18_all <- dropLayer(per_18_all, 13)
# names(per_18_all)[27] <- "dist_ag"

### Run prediction
##################

# ### Run prediction
# pp_peru <- raster::predict(per_18_all, fit_per,
#                              na.rm = TRUE,
#                              type = "response")
# 
# ### Open 2018 LUC map
# peru18 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_per_dry_2018_102033.tif")
# 
# ### Mask cells that were not forested in 2018
# pp_peru_mask <- mask(pp_peru, peru18,
#                      # filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m1_mask.tif",
#                      inverse= TRUE,
#                      maskvalue = 2)
# 
# writeRaster(pp_peru_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/peru_pp_masked.tif",
#                         format = "GTiff",
#                         overwrite = TRUE,
#                         options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

##################
### peru model 2
##################
# ### Run regression
# fit_per_2 <- glm(trans_rc ~ dist_ag +
#                    control + fire + illegal_mines +
#                    nonforconc + npchg + permprod + protforest +
#                    refor + tourism,
#                  data = per_dat_use,
#                  family = binomial(link = "logit"))
# 
# ### Run prediction
# ##################
# pp_peru_2 <- raster::predict(per_18_all, fit_per_2,
#                              na.rm = TRUE,
#                              type = "response")
#                              # filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m2_logit.tif",
#                              # overwrite = TRUE)
# 
# ### Mask cells that were not forested in 2018
# pp_peru_2_mask <- mask(pp_peru_2, peru18,
#                        inverse= TRUE,
#                        maskvalue = 2)
# 
# writeRaster(pp_peru_2_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/peru_pp_2_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

##################
### peru model 3 
##################
# ### Run regression
# fit_per_3 <- glm(trans_rc ~ aspect + slope + 
#                    crop_suit + prot_status + popdens + 
#                    precip + rivers + dist_ag + perc_diff +
#                    control + fire + illegal_mines +
#                    nonforconc + npchg + permprod + protforest +
#                    refor + tourism,
#                  data = per_dat_use,
#                  family = binomial(link = "logit"))
# 
# ### Run prediction
# ##################
# pp_peru_3 <- raster::predict(per_18_all, fit_per_3,
#                                na.rm = TRUE,
#                                type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_peru_3_mask <- mask(pp_peru_3, peru18,
#                          inverse= TRUE,
#                          maskvalue = 2)
# 
# writeRaster(pp_peru_3_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/peru_pp_3_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

##################
### peru model 4 
##################
# ### Run regression
# fit_per_4 <- glm(trans_rc ~ slope + rivers + 
#                    dist_ag + crop_suit + popdens + 
#                    perc_diff + prot_status +
#                    control + fire + illegal_mines + 
#                    nonforconc + permprod + protforest +
#                    refor + tourism,
#                  data = per_dat_use,
#                  family = binomial(link = "logit"))

# ### Run prediction
# ##################
# pp_peru_4 <- raster::predict(per_18_all, fit_per_4,
#                                na.rm = TRUE,
#                                type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_peru_4_mask <- mask(pp_peru_4, peru18,
#                          inverse= TRUE,
#                          maskvalue = 2)
# 
# writeRaster(pp_peru_4_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/peru_pp_4_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# 
# ######################
# ### Save model outputs 
# ######################
# 
# # ### Model 1
# fit_per_summ <- summary(fit_per)
# fit_per_coeff <- as.data.frame(fit_per_summ$coefficients)
# write.csv(fit_per_coeff,
#           file = "/nfs/agfrontiers-data/luc_model/fit_per_model1.csv")
# 
# ### Model 2
# fit_per_summ2 <- summary(fit_per_2)
# fit_per_coeff2 <- as.data.frame(fit_per_summ2$coefficients)
# write.csv(fit_per_coeff2,
#           file = "/nfs/agfrontiers-data/luc_model/fit_per_model2.csv")
# 
# ### Model 3
# fit_per_summ3 <- summary(fit_per_3)
# fit_per_coeff3 <- as.data.frame(fit_per_summ3$coefficients)
# write.csv(fit_per_coeff3,
#           file = "/nfs/agfrontiers-data/luc_model/fit_per_model3.csv")
# 
# ### Model 4
# fit_per_summ4 <- summary(fit_per_4)
# fit_per_coeff4 <- as.data.frame(fit_per_summ4$coefficients)
# write.csv(fit_per_coeff4,
#           file = "/nfs/agfrontiers-data/luc_model/fit_per_model4.csv")

##################
### Compare models 
##################

# ### anova
# anova(fit_per, fit_per_2, fit_per_3, fit_per_4, test ="Chisq")
# anova(fit_per_3, fit_per_4, test ="Chisq")
# 
# ### psuedo r2
# liperry(pscl)
# pR2(fit_per)
# pR2(fit_per_2)
# pR2(fit_per_3)
# pR2(fit_per_4)

##############################################################################
##############################################################################
##############################################################################

##################
### MC simulations 
##################

######################
### Model 1 projection
######################

### Open PA shps
# tambo <- st_read("/nfs/agfrontiers-data/Case Study Info/Peru/wdpa_mar2021_tambo_102033.shp")
# bahua <- st_read("/nfs/agfrontiers-data/Case Study Info/Peru/wdpa_mar2021_bahuaja_102033.shp")

# ### Open PP1 map
# pp_peru_mask <- raster("/nfs/agfrontiers-data/luc_model/peru_pp_masked.tif")

# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Tambo
# tambo_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Bahuaja
# bahua_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
# 
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp_peru_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
# 
#   ### Set CRS
#   j_crs <- crs(pp_peru_mask)
#   crs(j_blank) <- j_crs
# 
#   ### Set extent
#   extent(j_blank) <- extent(pp_peru_mask)
# 
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_peru_mask)
#   j_stack <- stack(j_blank, pp_peru_mask)
# 
#   ### Reclassification function (1 = forest converts, 2 = forest stays forest)
#   rc <- function(x1, x2) {
#     ifelse(x1 < x2, 1, 0)
#   }
# 
#   ### Make output raster
#   j_classified <- overlay(j_stack, fun = rc)
# 
#   ### Save raster
#   writeRaster(j_classified,
#               filename = paste0("/nfs/agfrontiers-data/luc_model/peru_1_project/peru_projected_",
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
#   ### Calculate area that converts
#   area_change <- freq(j_classified, value = 1)
#   forest_loss <- area_change * res(j_classified)[1] * res(j_classified)[2]
# 
#   ### Add to vector
#   area_vec[i] <- forest_loss
#   
#   #######################
#   ### Tambopata
#   #######################
#   ### Clip raster to shapefile
#   rast_tambo <- crop(j_classified, tambo)
#   rast_tambo <- mask(rast_tambo, tambo)
#   
#   ### Calculate # pixels that convert
#   forest_loss_tambo <- freq(rast_tambo, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_tambo <- forest_loss_tambo * res(rast_tambo)[1] * res(rast_tambo)[2]
#   
#   ### Add to vector
#   tambo_area_vec[i] <- forest_loss_tambo
#   
#   #######################
#   ### Bahuaja
#   #######################
#   ### Clip raster to shapefile
#   rast_bahua <- crop(j_classified, bahua)
#   rast_bahua <- mask(rast_bahua, bahua)
#   
#   ### Calculate # pixels that convert
#   forest_loss_bahua <- freq(rast_bahua, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_bahua <- forest_loss_bahua * res(rast_bahua)[1] * res(rast_bahua)[2]
#   
#   ### Add to vector
#   bahua_area_vec[i] <- forest_loss_bahua
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_1_project/peru_projected_loss.csv",
#           row.names=FALSE)
# write.csv(tambo_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_1_project/peru_projected_loss_tambopata.csv",
#           row.names=FALSE)
# write.csv(bahua_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_1_project/peru_projected_loss_bahuaja.csv",
#           row.names=FALSE)

######################
### Model 2 projection
######################
### Open PP2 map
# pp2_peru_mask <- raster("/nfs/agfrontiers-data/luc_model/peru_pp_2_masked.tif")

# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Tambo
# tambo_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Bahuaja
# bahua_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp2_peru_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp2_peru_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp2_peru_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp2_peru_mask)
#   j_stack <- stack(j_blank, pp2_peru_mask)
#   
#   ### Reclassification function (1 = forest converts, 2 = forest stays forest)
#   rc <- function(x1, x2) {
#     ifelse(x1 < x2, 1, 0)
#   }
#   
#   ### Make output raster
#   j_classified <- overlay(j_stack, fun = rc)
#   
#   ### Save raster
#   writeRaster(j_classified,
#               filename = paste0("/nfs/agfrontiers-data/luc_model/peru_2_project/peru_2_projected_",
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
#   
#   ### Calculate area that converts
#   area_change <- freq(j_classified, value = 1)
#   forest_loss <- area_change * res(j_classified)[1] * res(j_classified)[2]
#   
#   ### Add to vector
#   area_vec[i] <- forest_loss
#   
#   #######################
#   ### Tambopata
#   #######################
#   ### Clip raster to shapefile
#   rast_tambo <- crop(j_classified, tambo)
#   rast_tambo <- mask(rast_tambo, tambo)
#   
#   ### Calculate # pixels that convert
#   forest_loss_tambo <- freq(rast_tambo, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_tambo <- forest_loss_tambo * res(rast_tambo)[1] * res(rast_tambo)[2]
#   
#   ### Add to vector
#   tambo_area_vec[i] <- forest_loss_tambo
#   
#   #######################
#   ### Bahuaja
#   #######################
#   ### Clip raster to shapefile
#   rast_bahua <- crop(j_classified, bahua)
#   rast_bahua <- mask(rast_bahua, bahua)
#   
#   ### Calculate # pixels that convert
#   forest_loss_bahua <- freq(rast_bahua, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_bahua <- forest_loss_bahua * res(rast_bahua)[1] * res(rast_bahua)[2]
#   
#   ### Add to vector
#   bahua_area_vec[i] <- forest_loss_bahua
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_2_project/peru_2_projected_loss.csv",
#           row.names=FALSE)
# write.csv(tambo_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_2_project/peru_2_projected_loss_tambopata.csv",
#           row.names=FALSE)
# write.csv(bahua_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_2_project/peru_2_projected_loss_bahuaja.csv",
#           row.names=FALSE)

######################
### Model 3 projection
######################
### Open PP3 map
# pp3_peru_mask <- raster("/nfs/agfrontiers-data/luc_model/peru_pp_3_masked.tif")

# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Tambo
# tambo_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Bahuaja
# bahua_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp3_peru_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp3_peru_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp3_peru_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp3_peru_mask)
#   j_stack <- stack(j_blank, pp3_peru_mask)
#   
#   ### Reclassification function (1 = forest converts, 2 = forest stays forest)
#   rc <- function(x1, x2) {
#     ifelse(x1 < x2, 1, 0)
#   }
#   
#   ### Make output raster
#   j_classified <- overlay(j_stack, fun = rc)
#   
#   ### Save raster
#   writeRaster(j_classified,
#               filename = paste0("/nfs/agfrontiers-data/luc_model/peru_3_project/peru_3_projected_",
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
#   
#   ### Calculate area that converts
#   area_change <- freq(j_classified, value = 1)
#   forest_loss <- area_change * res(j_classified)[1] * res(j_classified)[2]
#   
#   ### Add to vector
#   area_vec[i] <- forest_loss
#   
#   #######################
#   ### Tambopata
#   #######################
#   ### Clip raster to shapefile
#   rast_tambo <- crop(j_classified, tambo)
#   rast_tambo <- mask(rast_tambo, tambo)
#   
#   ### Calculate # pixels that convert
#   forest_loss_tambo <- freq(rast_tambo, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_tambo <- forest_loss_tambo * res(rast_tambo)[1] * res(rast_tambo)[2]
#   
#   ### Add to vector
#   tambo_area_vec[i] <- forest_loss_tambo
#   
#   #######################
#   ### Bahuaja
#   #######################
#   ### Clip raster to shapefile
#   rast_bahua <- crop(j_classified, bahua)
#   rast_bahua <- mask(rast_bahua, bahua)
#   
#   ### Calculate # pixels that convert
#   forest_loss_bahua <- freq(rast_bahua, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_bahua <- forest_loss_bahua * res(rast_bahua)[1] * res(rast_bahua)[2]
#   
#   ### Add to vector
#   bahua_area_vec[i] <- forest_loss_bahua
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_3_project/peru_3_projected_loss.csv",
#           row.names=FALSE)
# write.csv(tambo_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_3_project/peru_3_projected_loss_tambopata.csv",
#           row.names=FALSE)
# write.csv(bahua_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_3_project/peru_3_projected_loss_bahuaja.csv",
#           row.names=FALSE)

######################
### Model 4 projection
######################
# ### Open PP4 map
# pp4_peru_mask <- raster("/nfs/agfrontiers-data/luc_model/peru_pp_4_masked.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Tambo
# tambo_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Bahuaja
# bahua_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp4_peru_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp4_peru_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp4_peru_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp4_peru_mask)
#   j_stack <- stack(j_blank, pp4_peru_mask)
#   
#   ### Reclassification function (1 = forest converts, 2 = forest stays forest)
#   rc <- function(x1, x2) {
#     ifelse(x1 < x2, 1, 0)
#   }
#   
#   ### Make output raster
#   j_classified <- overlay(j_stack, fun = rc)
#   
#   ### Save raster
#   writeRaster(j_classified,
#               filename = paste0("/nfs/agfrontiers-data/luc_model/peru_4_project/peru_4_projected_",
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
#   
#   ### Calculate area that converts
#   area_change <- freq(j_classified, value = 1)
#   forest_loss <- area_change * res(j_classified)[1] * res(j_classified)[2]
#   
#   ### Add to vector
#   area_vec[i] <- forest_loss
#   
#   #######################
#   ### Tambopata
#   #######################
#   ### Clip raster to shapefile
#   rast_tambo <- crop(j_classified, tambo)
#   rast_tambo <- mask(rast_tambo, tambo)
#   
#   ### Calculate # pixels that convert
#   forest_loss_tambo <- freq(rast_tambo, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_tambo <- forest_loss_tambo * res(rast_tambo)[1] * res(rast_tambo)[2]
#   
#   ### Add to vector
#   tambo_area_vec[i] <- forest_loss_tambo
#   
#   #######################
#   ### Bahuaja
#   #######################
#   ### Clip raster to shapefile
#   rast_bahua <- crop(j_classified, bahua)
#   rast_bahua <- mask(rast_bahua, bahua)
#   
#   ### Calculate # pixels that convert
#   forest_loss_bahua <- freq(rast_bahua, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_bahua <- forest_loss_bahua * res(rast_bahua)[1] * res(rast_bahua)[2]
#   
#   ### Add to vector
#   bahua_area_vec[i] <- forest_loss_bahua
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_4_project/peru_4_projected_loss.csv",
#           row.names=FALSE)
# write.csv(tambo_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_4_project/peru_4_projected_loss_tambopata.csv",
#           row.names=FALSE)
# write.csv(bahua_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/peru_4_project/peru_4_projected_loss_bahuaja.csv",
#           row.names=FALSE)

##################################################################
### Stack layers

###################
########## Model 1
###################
# ### Try on subset of layers in Peru Model 1
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_1_project",
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
#             filename = "/nfs/agfrontiers-data/luc_model/peru_m1_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

# ###################
# ########## Model 2
# ###################
# ### Try on subset of layers in Peru Model 2
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_2_project",
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
#             filename = "/nfs/agfrontiers-data/luc_model/peru_m2_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

# ###################
# ########## Model 3
# ###################
# ### Try on subset of layers in Peru Model 3
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/peru_3_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# # test_layers_1 <- test_layers[1:200]
# # test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
# # test_layers_4 <- test_layers[601:800]
# # test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
# ### Write loop
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
# 
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/peru_m3_stacked_3.tif",
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
test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop
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
            filename = "/nfs/agfrontiers-data/luc_model/peru_m4_stacked_4.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

#############################
### Run predictions 2008-2018
#############################

##################
### Peru model 1 
##################

### Data prep
#############

### Open sample points
per_dat_all <- read.csv("/nfs/agfrontiers-data/luc_model/peru_data/peru_data_for_model.csv")

### Make column for transition/none
per_dat_all <- per_dat_all %>%
  mutate(transition = ifelse(lc_2018 == "2",
                             "no_change",
                             ifelse(lc_2018 == "1",
                                    "f_to_ag",
                                    "other_change")),
         trans_rc = ifelse(transition == "no_change",
                           "0", ifelse(transition == "f_to_ag",
                                       "1", "2")))

### Make columns factor
fact_cols <- c("uid", "prot_status",
               "minconc",
               "nonforconc", "npchg",
               "permprod", "protforest",
               "refor",
               "lc_2018", "trans_rc")
per_dat_all <- per_dat_all %>%
  mutate_each_(funs(factor(.)), fact_cols)

### Restrict to transitions
per_dat_use <- per_dat_all %>%
  filter(trans_rc == "1" | trans_rc == "0")

### Model 1
###########

### Run regression
fit_per <- glm(trans_rc ~ aspect + cities +
                 slope +
                 crop_suit + prot_status + popdens +
                 precip + rivers + dist_ag + perc_diff,
               data = per_dat_use,
               family = binomial(link = "logit"))

### 2008 data
#############

### Open rasterbrick
per_layers_all <- brick("/nfs/agfrontiers-data/luc_model/peru_all_layers.tif")

### Rename columns
names(per_layers_all) <- c("aspect", "cities",
                           "crop_suit", "elev",
                           "mines", "prot_status",
                           "popdens", "poverty",
                           "precip", "rivers",
                           "roads", "slope",
                           "soil", "dist_ag",
                           "control", "commun",
                           "fire", "forest",
                           "illegal_mines", "minconc",
                           "nonforconc", "npchg",
                           "permprod", "protforest",
                           "refor", "tourism")

### Replace precip layer
per_layers_all <- dropLayer(per_layers_all, 9)

### Open new precip layer
per_precip <- raster("/nfs/agfrontiers-data/luc_model/per_precip_new.tif")

### Add new precip layer
per_layers_all <- addLayer(per_layers_all, per_precip)
names(per_layers_all)[26] <- "precip"

### Open 2008 percent of neighboring pixels raster
neigh_2008 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/per_perc_diff.tif")

### Add to other layers
per_layers_all <- addLayer(per_layers_all, neigh_2008)
names(per_layers_all)[27] <- "perc_diff"

### Run prediction
##################

### Run prediction
pp_peru <- raster::predict(per_layers_all, fit_per,
                              na.rm = TRUE,
                              type = "response",
                              overwrite = TRUE)

### Open 2008 LUC map
peru08 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_per_dry_2008_102033.tif")

### Mask cells that were not forested in 2008
pp_peru_mask_08 <- mask(pp_peru, peru08,
                           inverse= TRUE,
                           maskvalue = 2)

writeRaster(pp_peru_mask_08,
            filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m1_masked_2018.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###########
### Model 2
###########
### Run regression
fit_per_2 <- glm(trans_rc ~ dist_ag +
                   control + fire + illegal_mines +
                   nonforconc + npchg + permprod + protforest +
                   refor + tourism,
                 data = per_dat_use,
                 family = binomial(link = "logit"))

### Run prediction
pp_peru_2 <- raster::predict(per_layers_all, fit_per_2,
                                na.rm = TRUE,
                                type = "response",
                                # filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m2_2018.tif",
                                overwrite = TRUE)

### Mask cells that were not forested in 2008
pp_peru_2_mask_08 <- mask(pp_peru_2, peru08,
                             inverse= TRUE,
                             maskvalue = 2)

writeRaster(pp_peru_2_mask_08,
            filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m2_masked_2018.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###########
### Model 3
###########
### Run regression
fit_per_3 <- glm(trans_rc ~ aspect + slope +
                   crop_suit + prot_status + popdens +
                   precip + rivers + dist_ag + perc_diff +
                   control + fire + illegal_mines +
                   nonforconc + npchg + permprod + protforest +
                   refor + tourism,
                 data = per_dat_use,
                 family = binomial(link = "logit"))

### Run prediction
pp_peru_3 <- raster::predict(per_layers_all, fit_per_3,
                             na.rm = TRUE,
                             type = "response",
                             # filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m3_2018.tif",
                             overwrite = TRUE)

### Mask cells that were not forested in 2008
pp_peru_3_mask_08 <- mask(pp_peru_3, peru08,
                          inverse= TRUE,
                          maskvalue = 2)

writeRaster(pp_peru_3_mask_08,
            filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m3_masked_2018.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###########
### Model 4
###########
### Run regression
fit_per_4 <- glm(trans_rc ~ slope + rivers +
                   dist_ag + crop_suit + popdens +
                   perc_diff + prot_status +
                   control + fire + illegal_mines +
                   nonforconc + permprod + protforest +
                   refor + tourism,
                 data = per_dat_use,
                 family = binomial(link = "logit"))


### Run prediction
pp_peru_4 <- raster::predict(per_layers_all, fit_per_4,
                             na.rm = TRUE,
                             type = "response",
                             # filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m4_2018.tif",
                             overwrite = TRUE)

### Mask cells that were not forested in 2008
pp_peru_4_mask_08 <- mask(pp_peru_4, peru08,
                          inverse= TRUE,
                          maskvalue = 2)

writeRaster(pp_peru_4_mask_08,
            filename = "/nfs/agfrontiers-data/luc_model/peru_pp_m4_masked_2018.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
