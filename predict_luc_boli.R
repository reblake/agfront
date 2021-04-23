### Code to make LUC models and apply to 2018 LUC maps

### Open packages
library(tidyverse)
library(sf)
library(raster)
# library(spatialEco)

# ####################
# ### Bolivia model 1 
# ####################
# 
# ### Open sample points
# bol_dat_all <- read.csv("/nfs/agfrontiers-data/luc_model/boli_data_for_model.csv")
# 
# ### Rename columns
# colnames(bol_dat_all) <- c("uid", "lon", "lat", "lc_2008", "lc_2018", 
#                            "id", "mines", "roads", "rivers", "cities",     
#                            "popdens", "poverty", "prot_status", "crop_suit", 
#                            "dist_ag", "dist_fire", "fire_dens", 
#                            "ill_mine", "field_res", "tourism", 
#                            "pes", "precip", "soil", "south_north", 
#                            "landtenure", "elev", "slope", "aspect", "perc_diff")
# 
# ### Make column for transition/none
# bol_dat_all <- bol_dat_all %>%
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
#                "pes", "south_north",
#                "landtenure",
#                "lc_2018", "trans_rc")
# bol_dat_all <- bol_dat_all %>%
#   mutate_each_(funs(factor(.)), fact_cols)
# 
# ### Model 1
# ###########
# 
# ### Run regression
# fit_bol <- glm(trans_rc ~ elev + aspect + slope + poverty +
#                  cities + mines + crop_suit + prot_status + popdens + 
#                  precip + rivers + dist_ag + perc_diff,
#                data = bol_dat_use,
#                family = binomial(link = "logit"))
# 
# ### Write out coefficients
# fit_bol_summ <- summary(fit_bol)
# fit_bol_coeff <- as.data.frame(fit_bol_summ$coefficients)
# write.csv(fit_bol_coeff, 
#           file = "/nfs/agfrontiers-data/luc_model/fit_bol_model1.csv")
# 
# ### 2018 data
# #############
# 
# ### Open rasterbrick
# bol_layers_all <- brick("/nfs/agfrontiers-data/luc_model/bolivia_all_vars.tif")
# 
# ### Rename columns
# names(bol_layers_all) <- c("aspect", "mines", "roads", "rivers",
#                            "cities", "elev", "popdens", "poverty", 
#                            "ppt", "prot_status", "slope", "soil", 
#                            "crop_suit", "dist_ag", "dist_fire", "fire_dens",
#                            "ill_mine", "field_res", "tourism", "pes")
# 
# ### Replace precip layer
# bol_layers_all <- dropLayer(bol_layers_all, 9)
# 
# ### Open new precip layer
# bol_precip <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_precip_use.tif")
# 
# ### Add new precip layer
# bol_layers_all <- addLayer(bol_layers_all, bol_precip)
# 
# ### Drop old soil moisture layer
# bol_layers_all <- dropLayer(bol_layers_all, 11)
# 
# ### Open new soil moisture layer
# bol_sm <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/sm_bol_use.tif")
# bol_sm <- raster::resample(bol_sm, bol_precip, "bilinear")
# 
# ### Add new precip layer
# bol_layers_all <- addLayer(bol_layers_all, bol_sm)
# 
# ### Add north-south layer
# bol_ns <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_south_north.tif")
# bol_layers_all <- addLayer(bol_layers_all, bol_ns)
# 
# ### Add land tenure layer
# bol_landten <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_landtenure.tif")
# bol_layers_all <- addLayer(bol_layers_all, bol_landten)
# 
# ### Drop DEM layers (wrong extent)
# bol_layers_all <- dropLayer(bol_layers_all, c(1, 6, 10))
# 
# ### Replace DEM layers
# bol_elev <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_elev_326.tif")
# bol_slope <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_slope_326.tif")
# bol_aspect <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_asp_326.tif")
# 
# ### Add DEM layers
# bol_layers_all <- addLayer(bol_layers_all, bol_elev)
# bol_layers_all <- addLayer(bol_layers_all, bol_slope)
# bol_layers_all <- addLayer(bol_layers_all, bol_aspect)
# 
# ### Add cattle
# bol_cattle <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/cattle_ambcar.tiff")
# bol_cattle_res <- raster::resample(bol_cattle, bol_precip, "ngb")
# bol_layers_all <- addLayer(bol_layers_all, bol_cattle_res)
# 
# ### Rename layers
# names(bol_layers_all) <- c("mines", "roads", "rivers", "cities", 
#                            "popdens", "poverty", "prot_status", "crop_suit", 
#                            "dist_ag", "dist_fire", "fire_dens", "ill_mine", 
#                            "field_res", "tourism", "pes", "precip", 
#                            "soil", "south_north", "landtenure", "elev",
#                            "slope", "aspect", "cattle")
# 
# ## Convert NAs to 0 in PES layer
# bol_layers_all[[15]][is.na(bol_layers_all[[15]][])] <- 0
# 
# ### Open 2018 percent of neighboring pixels raster
# neigh_2018 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/bol18_perc_diff.tif")
# 
# ### Add to other layers
# bol_layers_all <- addLayer(bol_layers_all, neigh_2018)
# names(bol_layers_all)[24] <- "perc_diff"
# 
# ### Run prediction
# ##################
# 
# ### Run prediction
# pp_boli <- raster::predict(bol_layers_all, fit_bol,
#                              na.rm = TRUE,
#                              type = "response")
# 
# ### Open 2018 LUC map
# boli18 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bol_dry_2018_102033.tif")
# 
# ### Mask cells that were not forested in 2018
# pp_boli_mask <- mask(pp_boli, boli18,
#                      # filename = "/nfs/agfrontiers-data/luc_model/boli_pp_m1_mask.tif",
#                      inverse= TRUE,
#                      maskvalue = 2)
# 
# writeRaster(pp_boli_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_pp_masked.tif",
#                         format = "GTiff",
#                         overwrite = TRUE,
#                         options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# writeRaster(pp_boli_mask,
#             filename = "/nfs/agfrontiers-data/boli_pp_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ###########
# ### Model 2
# ###########
# ### Run regression
# fit_bol_m2 <- glm(trans_rc ~ dist_ag + dist_fire +
#                     fire_dens + ill_mine + field_res +
#                     tourism + pes + south_north +
#                     landtenure + cattle,
#                   data = bol_dat_use,
#                   family = binomial(link = "logit"))
# 
# ### Write out coefficients
# fit_bol2_summ <- summary(fit_bol_m2)
# fit_bol2_coeff <- as.data.frame(fit_bol2_summ$coefficients)
# write.csv(fit_bol2_coeff, 
#           file = "/nfs/agfrontiers-data/luc_model/fit_bol_model2.csv")
# 
# ### Run prediction
# pp_boli_2 <- raster::predict(bol_layers_all, fit_bol_m2,
#                              na.rm = TRUE,
#                              type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_boli_2_mask <- mask(pp_boli_2, boli18,
#                      # filename = "/nfs/agfrontiers-data/luc_model/boli_pp_m2_mask.tif",
#                      inverse= TRUE,
#                      maskvalue = 2)
# 
# writeRaster(pp_boli_2_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_pp_m2_mask.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# 
# ###########
# ### Model 3
# ###########
# ### Run regression
# fit_bol_m3 <- glm(trans_rc ~ elev + aspect + slope + soil +
#                     cities + mines + crop_suit + prot_status + popdens + 
#                     precip + rivers + dist_ag + perc_diff +
#                     dist_fire + fire_dens + ill_mine + field_res +
#                     tourism + pes + south_north + landtenure + cattle,
#                   data = bol_dat_use,
#                   family = binomial(link = "logit"))
# 
# ### Write out coefficients
# fit_bol3_summ <- summary(fit_bol_m3)
# fit_bol3_coeff <- as.data.frame(fit_bol3_summ$coefficients)
# write.csv(fit_bol3_coeff, 
#           file = "/nfs/agfrontiers-data/luc_model/fit_bol_model3.csv")
# 
# ### Run prediction
# pp_boli_3 <- raster::predict(bol_layers_all, fit_bol_m3,
#                              na.rm = TRUE,
#                              type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_boli_3_mask <- mask(pp_boli_3, boli18,
#                        inverse= TRUE,
#                        maskvalue = 2)
# 
# writeRaster(pp_boli_3_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_pp_m3_mask.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ###########
# ### Model 4
# ###########
# ### Run regression
# fit_bol_m4 <- glm(trans_rc ~ elev + aspect + slope + soil +
#                     cities + mines + crop_suit + precip +
#                     prot_status + dist_ag + perc_diff +
#                     rivers + dist_fire + fire_dens + pes + 
#                     south_north + landtenure + tourism + cattle,
#                   data = bol_dat_use,
#                   family = binomial(link = "logit"))
# 
# ### Write out coefficients
# fit_bol4_summ <- summary(fit_bol_m4)
# fit_bol4_coeff <- as.data.frame(fit_bol4_summ$coefficients)
# write.csv(fit_bol4_coeff, 
#           file = "/nfs/agfrontiers-data/luc_model/fit_bol_model4.csv")
# 
# ### Run prediction
# pp_boli_4 <- raster::predict(bol_layers_all, fit_bol_m4,
#                              na.rm = TRUE,
#                              type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_boli_4_mask <- mask(pp_boli_4, boli18,
#                        inverse= TRUE,
#                        maskvalue = 2)
# 
# writeRaster(pp_boli_4_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_pp_m4_mask.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###############################
### Monte Carlo
##############################

# ### Open PA shps
# carra <- st_read("/nfs/agfrontiers-data/Case Study Info/Bolivia/wdpa_apr21_carrasco_102033.shp")
# amb_np <- st_read("/nfs/agfrontiers-data/Case Study Info/Bolivia/wdpa_apr21_amboro_np_102033.shp")
# amb_im <- st_read("/nfs/agfrontiers-data/Case Study Info/Bolivia/wdpa_apr21_amboro_im_102033.shp")

# #######
# ## M1
# #######
# ### Open PP1 map
# pp_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_masked.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro IMNA
# ambim_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
# 
#   ### Make blank raster, fill with random values
#   j_blank <- random.raster(pp_boli_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
# 
#   ### Set CRS
#   j_crs <- crs(pp_boli_mask)
#   crs(j_blank) <- j_crs
# 
#   ### Set extent
#   extent(j_blank) <- extent(pp_boli_mask)
# 
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_boli_mask)
#   j_stack <- stack(j_blank, pp_boli_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/boliv_1_project/boli_projected_",
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
#   ### Carrasco
#   #######################
#   ### Clip raster to shapefile
#   rast_carra <- crop(j_classified, carra)
#   rast_carra <- mask(rast_carra, carra)
# 
#   ### Calculate # pixels that convert
#   forest_loss_carra <- freq(rast_carra, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]
# 
#   ### Add to vector
#   carra_area_vec[i] <- forest_loss_carra
# 
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambnp <- crop(j_classified, amb_np)
#   rast_ambnp <- mask(rast_ambnp, amb_np)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambnp <- freq(rast_ambnp, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]
#   
#   ### Add to vector
#   ambnp_area_vec[i] <- forest_loss_ambnp
# 
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambim <- crop(j_classified, amb_im)
#   rast_ambim <- mask(rast_ambim, amb_im)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambim <- freq(rast_ambim, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]
#   
#   ### Add to vector
#   ambim_area_vec[i] <- forest_loss_ambim
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/boliv_1_project/boli_projected_loss.csv",
#           row.names=FALSE)
# write.csv(carra_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/boliv_1_project/boli_projected_loss_carrasco.csv",
#           row.names=FALSE)
# write.csv(ambnp_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/boliv_1_project/boli_projected_loss_amboro_np.csv",
#           row.names=FALSE)
# write.csv(ambim_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/boliv_1_project/boli_projected_loss_amboro_imna.csv",
#           row.names=FALSE)

# #######
# ## M2
# #######
# ### Open PP2 map
# pp_2_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_m2_mask.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro IMNA
# ambim_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
# 
#   ### Make blank raster, fill with random values
#   j_blank <- random.raster(pp_2_boli_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
# 
#   ### Set CRS
#   j_crs <- crs(pp_2_boli_mask)
#   crs(j_blank) <- j_crs
# 
#   ### Set extent
#   extent(j_blank) <- extent(pp_2_boli_mask)
# 
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_2_boli_mask)
#   j_stack <- stack(j_blank, pp_2_boli_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_",
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
#   ### Carrasco
#   #######################
#   ### Clip raster to shapefile
#   rast_carra <- crop(j_classified, carra)
#   rast_carra <- mask(rast_carra, carra)
# 
#   ### Calculate # pixels that convert
#   forest_loss_carra <- freq(rast_carra, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]
# 
#   ### Add to vector
#   carra_area_vec[i] <- forest_loss_carra
# 
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambnp <- crop(j_classified, amb_np)
#   rast_ambnp <- mask(rast_ambnp, amb_np)
# 
#   ### Calculate # pixels that convert
#   forest_loss_ambnp <- freq(rast_ambnp, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]
# 
#   ### Add to vector
#   ambnp_area_vec[i] <- forest_loss_ambnp
# 
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambim <- crop(j_classified, amb_im)
#   rast_ambim <- mask(rast_ambim, amb_im)
# 
#   ### Calculate # pixels that convert
#   forest_loss_ambim <- freq(rast_ambim, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]
# 
#   ### Add to vector
#   ambim_area_vec[i] <- forest_loss_ambim
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss.csv",
#           row.names=FALSE)
# write.csv(carra_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss_carrasco.csv",
#           row.names=FALSE)
# write.csv(ambnp_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss_amboro_np.csv",
#           row.names=FALSE)
# write.csv(ambim_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss_amboro_imna.csv",
#           row.names=FALSE)

#######
## M3
#######
# ### Open PP3 map
# pp_3_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_m3_mask.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro IMNA
# ambim_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank raster, fill with random values
#   j_blank <- random.raster(pp_3_boli_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp_3_boli_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp_3_boli_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_3_boli_mask)
#   j_stack <- stack(j_blank, pp_3_boli_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_",
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
#   ### Carrasco
#   #######################
#   ### Clip raster to shapefile
#   rast_carra <- crop(j_classified, carra)
#   rast_carra <- mask(rast_carra, carra)
#   
#   ### Calculate # pixels that convert
#   forest_loss_carra <- freq(rast_carra, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]
#   
#   ### Add to vector
#   carra_area_vec[i] <- forest_loss_carra
#   
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambnp <- crop(j_classified, amb_np)
#   rast_ambnp <- mask(rast_ambnp, amb_np)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambnp <- freq(rast_ambnp, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]
#   
#   ### Add to vector
#   ambnp_area_vec[i] <- forest_loss_ambnp
#   
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambim <- crop(j_classified, amb_im)
#   rast_ambim <- mask(rast_ambim, amb_im)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambim <- freq(rast_ambim, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]
#   
#   ### Add to vector
#   ambim_area_vec[i] <- forest_loss_ambim
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss.csv",
#           row.names=FALSE)
# write.csv(carra_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss_carrasco.csv",
#           row.names=FALSE)
# write.csv(ambnp_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss_amboro_np.csv",
#           row.names=FALSE)
# write.csv(ambim_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss_amboro_imna.csv",
#           row.names=FALSE)

#######
## M4
#######
# ### Open PP4 map
# pp_4_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_m4_mask.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro IMNA
# ambim_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank raster, fill with random values
#   j_blank <- random.raster(pp_4_boli_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp_4_boli_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp_4_boli_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_4_boli_mask)
#   j_stack <- stack(j_blank, pp_4_boli_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/bolivia_4_project/boli_4_projected_",
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
#   ### Carrasco
#   #######################
#   ### Clip raster to shapefile
#   rast_carra <- crop(j_classified, carra)
#   rast_carra <- mask(rast_carra, carra)
#   
#   ### Calculate # pixels that convert
#   forest_loss_carra <- freq(rast_carra, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]
#   
#   ### Add to vector
#   carra_area_vec[i] <- forest_loss_carra
#   
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambnp <- crop(j_classified, amb_np)
#   rast_ambnp <- mask(rast_ambnp, amb_np)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambnp <- freq(rast_ambnp, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]
#   
#   ### Add to vector
#   ambnp_area_vec[i] <- forest_loss_ambnp
#   
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambim <- crop(j_classified, amb_im)
#   rast_ambim <- mask(rast_ambim, amb_im)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambim <- freq(rast_ambim, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]
#   
#   ### Add to vector
#   ambim_area_vec[i] <- forest_loss_ambim
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_4_project/boli_4_projected_loss.csv",
#           row.names=FALSE)
# write.csv(carra_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_4_project/boli_4_projected_loss_carrasco.csv",
#           row.names=FALSE)
# write.csv(ambnp_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_4_project/boli_4_projected_loss_amboro_np.csv",
#           row.names=FALSE)
# write.csv(ambim_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_4_project/boli_4_projected_loss_amboro_imna.csv",
#           row.names=FALSE)

##################################################################
### Stack layers

###################
########## Model 1
###################
# ### Try on subset of layers in Bolivia Model 1
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/boliv_1_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m1_stacked_1.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m1_stacked_2.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m1_stacked_3.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 4
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m1_stacked_4.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 5
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m1_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###################
########## Model 2
###################
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/bolivia_2_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m2_stacked_1.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m2_stacked_2.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m2_stacked_3.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 4
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m2_stacked_4.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 5
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m2_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###################
########## Model 3
###################
# ### Try on subset of layers in Bolivia Model 3
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/bolivia_3_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_1.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_2.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
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
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_3.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 4
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_4.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 5
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m3_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ###################
# ########## Model 4
# ###################
# ### Try on subset of layers in Bolivia Model 4
# test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/bolivia_4_project",
#                           pattern = "*.tif",
#                           full.names = TRUE)
# 
# test_layers_1 <- test_layers[1:200]
# test_layers_2 <- test_layers[201:400]
# test_layers_3 <- test_layers[401:600]
# test_layers_4 <- test_layers[601:800]
# test_layers_5 <- test_layers[801:1000]
# 
# ### Blank raster for adding all rasters
# raster_fill <- raster(test_layers[1])
# raster_fill[raster_fill > 0] <- 0
# 
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
# 
# ### Write loop for stack 4
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_4.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ### Write loop for stack 5
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
# ### Write out sum raster
# writeRaster(raster_fill,
#             filename = "/nfs/agfrontiers-data/luc_model/boli_m4_stacked_5.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

#############################
### Run predictions 2008-2018
#############################

##################
### Bolivia model 1 
##################

### Data prep
#############

# ### Open sample points
# bol_dat_all <- read.csv("/nfs/agfrontiers-data/luc_model/boli_data/boli_data_for_model.csv")
# 
# ### Rename columns
# colnames(bol_dat_all) <- c("uid", "lon", "lat", "lc_2008", "lc_2018",
#                            "id", "mines", "roads", "rivers", "cities",
#                            "popdens", "poverty", "prot_status", "crop_suit",
#                            "dist_ag", "dist_fire", "fire_dens",
#                            "ill_mine", "field_res", "tourism",
#                            "pes", "precip", "soil", "south_north",
#                            "landtenure", "elev", "slope", "aspect", "cattle", "perc_diff")
# 
# ### Make column for transition/none
# bol_dat_all <- bol_dat_all %>%
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
#                "pes", "south_north",
#                "landtenure",
#                "lc_2018", "trans_rc")
# bol_dat_all <- bol_dat_all %>%
#   mutate_each_(funs(factor(.)), fact_cols)
# 
# ### Restrict to transitions
# bol_dat_use <- bol_dat_all %>%
#   filter(trans_rc == "1" | trans_rc == "0")
# 
# ### Model 1
# ###########
# 
# ### Run regression
# fit_bol <- glm(trans_rc ~ elev + aspect + slope + poverty +
#                  cities + mines + crop_suit + prot_status + popdens +
#                  precip + rivers + dist_ag + perc_diff,
#                data = bol_dat_use,
#                family = binomial(link = "logit"))
# 
# ### 2008 data
# #############
# 
# ### Open rasterbrick
# bol_layers_all <- brick("/nfs/agfrontiers-data/luc_model/bolivia_all_vars.tif")
# 
# ### Rename columns
# names(bol_layers_all) <- c("aspect", "mines", "roads", "rivers",
#                            "cities", "elev", "popdens", "poverty",
#                            "ppt", "prot_status", "slope", "soil",
#                            "crop_suit", "dist_ag", "dist_fire", "fire_dens",
#                            "ill_mine", "field_res", "tourism", "pes")
# 
# ### Replace precip layer
# bol_layers_all <- dropLayer(bol_layers_all, 9)
# 
# ### Open new precip layer
# bol_precip <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_precip_use.tif")
# 
# ### Add new precip layer
# bol_layers_all <- addLayer(bol_layers_all, bol_precip)
# 
# ### Drop old soil moisture layer
# bol_layers_all <- dropLayer(bol_layers_all, 11)
# 
# ### Open new soil moisture layer
# bol_sm <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/sm_bol_use.tif")
# bol_sm <- raster::resample(bol_sm, bol_precip, "bilinear")
# 
# ### Add new precip layer
# bol_layers_all <- addLayer(bol_layers_all, bol_sm)
# 
# ### Add north-south layer
# bol_ns <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_south_north.tif")
# bol_layers_all <- addLayer(bol_layers_all, bol_ns)
# 
# ### Add land tenure layer
# bol_landten <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_landtenure.tif")
# bol_layers_all <- addLayer(bol_layers_all, bol_landten)
# 
# ### Drop DEM layers (wrong extent)
# bol_layers_all <- dropLayer(bol_layers_all, c(1, 6, 10))
# 
# ### Replace DEM layers
# bol_elev <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_elev_326.tif")
# bol_slope <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_slope_326.tif")
# bol_aspect <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/bol_asp_326.tif")
# 
# ### Add DEM layers
# bol_layers_all <- addLayer(bol_layers_all, bol_elev)
# bol_layers_all <- addLayer(bol_layers_all, bol_slope)
# bol_layers_all <- addLayer(bol_layers_all, bol_aspect)
# 
# ### Add cattle
# bol_cattle <- raster("/nfs/agfrontiers-data/DINAMICA/Simple Model/Bolivia/Processed Data/cattle_ambcar.tiff")
# bol_cattle_res <- raster::resample(bol_cattle, bol_precip, "ngb")
# bol_layers_all <- addLayer(bol_layers_all, bol_cattle_res)
# 
# ### Rename layers
# names(bol_layers_all) <- c("mines", "roads", "rivers", "cities",
#                            "popdens", "poverty", "prot_status", "crop_suit",
#                            "dist_ag", "dist_fire", "fire_dens", "ill_mine",
#                            "field_res", "tourism", "pes", "precip",
#                            "soil", "south_north", "landtenure", "elev",
#                            "slope", "aspect", "cattle")
# 
# ## Convert NAs to 0 in PES layer
# bol_layers_all[[15]][is.na(bol_layers_all[[15]][])] <- 0
# 
# ### Open 2008 percent of neighboring pixels raster
# neigh_2008 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/bol_perc_diff.tif")
# 
# ### Add to other layers
# bol_layers_all <- addLayer(bol_layers_all, neigh_2008)
# names(bol_layers_all)[24] <- "perc_diff"
# 
# ### Run prediction
# ##################
# 
# ### Run prediction
# pp_bolivia <- raster::predict(bol_layers_all, fit_bol,
#                              na.rm = TRUE,
#                              type = "response",
#                              filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m1_2018.tif",
#                              overwrite = TRUE)
# 
# ### Open 2008 LUC map
# boli08 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bol_dry_2008_102033.tif")
# 
# ### Mask cells that were not forested in 2008
# pp_bolivia_mask_08 <- mask(pp_bolivia, boli08,
#                           inverse= TRUE,
#                           maskvalue = 2)
# 
# writeRaster(pp_bolivia_mask_08,
#             filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m1_masked_2018.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ###########
# ### Model 2
# ###########
# ### Run regression
# fit_bol_m2 <- glm(trans_rc ~ dist_ag + dist_fire +
#                     fire_dens + ill_mine + field_res +
#                     tourism + pes + south_north +
#                     landtenure + cattle,
#                   data = bol_dat_use,
#                   family = binomial(link = "logit"))
# 
# ### Run prediction
# pp_bolivia_2 <- raster::predict(bol_layers_all, fit_bol_m2,
#                               na.rm = TRUE,
#                               type = "response",
#                               filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m1_2018.tif",
#                               overwrite = TRUE)
# 
# ### Mask cells that were not forested in 2008
# pp_bolivia_2_mask_08 <- mask(pp_bolivia_2, boli08,
#                            inverse= TRUE,
#                            maskvalue = 2)
# 
# writeRaster(pp_bolivia_2_mask_08,
#             filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m2_masked_2018.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ###########
# ### Model 3
# ###########
# ### Run regression
# fit_bol_m3 <- glm(trans_rc ~ elev + aspect + slope + soil +
#                     cities + mines + crop_suit + prot_status + popdens +
#                     precip + rivers + dist_ag + perc_diff +
#                     dist_fire + fire_dens + ill_mine + field_res +
#                     tourism + pes + south_north + landtenure + cattle,
#                   data = bol_dat_use,
#                   family = binomial(link = "logit"))
# 
# ### Run prediction
# pp_bolivia_3 <- raster::predict(bol_layers_all, fit_bol_m3,
#                                 na.rm = TRUE,
#                                 type = "response",
#                                 # filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m1_2018.tif",
#                                 overwrite = TRUE)
# 
# ### Mask cells that were not forested in 2008
# pp_bolivia_3_mask_08 <- mask(pp_bolivia_3, boli08,
#                              inverse= TRUE,
#                              maskvalue = 2)
# 
# writeRaster(pp_bolivia_3_mask_08,
#             filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m3_masked_2018.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
# ###########
# ### Model 4
# ###########
# ### Run regression
# fit_bol_m4 <- glm(trans_rc ~ elev + aspect + slope + soil +
#                     cities + mines + crop_suit + precip +
#                     prot_status + dist_ag + perc_diff +
#                     rivers + dist_fire + fire_dens + pes +
#                     south_north + landtenure + tourism + cattle,
#                   data = bol_dat_use,
#                   family = binomial(link = "logit"))
# 
# ### Run prediction
# pp_bolivia_4 <- raster::predict(bol_layers_all, fit_bol_m4,
#                                 na.rm = TRUE,
#                                 type = "response",
#                                 # filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m1_2018.tif",
#                                 overwrite = TRUE)
# 
# ### Mask cells that were not forested in 2008
# pp_bolivia_4_mask_08 <- mask(pp_bolivia_4, boli08,
#                              inverse= TRUE,
#                              maskvalue = 2)
# 
# writeRaster(pp_bolivia_4_mask_08,
#             filename = "/nfs/agfrontiers-data/luc_model/boliv_pp_m4_masked_2018.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

###############################
### Monte Carlo 2018
##############################

### Open PA shps
carra <- st_read("/nfs/agfrontiers-data/Case Study Info/Bolivia/wdpa_apr21_carrasco_102033.shp")
amb_np <- st_read("/nfs/agfrontiers-data/Case Study Info/Bolivia/wdpa_apr21_amboro_np_102033.shp")
amb_im <- st_read("/nfs/agfrontiers-data/Case Study Info/Bolivia/wdpa_apr21_amboro_im_102033.shp")

#######
## M1
#######
### Open PP1 map
pp_boli_mask_2018 <- raster("/nfs/agfrontiers-data/luc_model/boliv_pp_m1_masked_2018.tif")

### Vector for area converted from forest to ag
area_vec <- rep("", times = 1000)

### Vector for area converted from forest to ag within Carrasco
carra_area_vec <- rep("", times = 1000)

### Vector for area converted from forest to ag within Amboro NP
ambnp_area_vec <- rep("", times = 1000)

### Vector for area converted from forest to ag within Amboro IMNA
ambim_area_vec <- rep("", times = 1000)

for (i in 1:length(area_vec)) {

  ### Make blank raster, fill with random values
  j_blank <- random.raster(pp_boli_mask_2018,
                           min = 0,
                           max = 1,
                           distribution = "random")

  ### Set CRS
  j_crs <- crs(pp_boli_mask_2018)
  crs(j_blank) <- j_crs

  ### Set extent
  extent(j_blank) <- extent(pp_boli_mask_2018)

  ### Combine with model 1 PP map
  j_blank <- crop(j_blank, pp_boli_mask_2018)
  j_stack <- stack(j_blank, pp_boli_mask_2018)

  ### Reclassification function (1 = forest converts, 2 = forest stays forest)
  rc <- function(x1, x2) {
    ifelse(x1 < x2, 1, 0)
  }

  ### Make output raster
  j_classified <- overlay(j_stack, fun = rc)

  ### Save raster
  writeRaster(j_classified,
              filename = paste0("/nfs/agfrontiers-data/luc_model/boliv_1_project_2018/boli_projected_2018_",
                                i, ".tif"),
              format = "GTiff",
              overwrite = TRUE,
              options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

  ### Calculate area that converts
  area_change <- freq(j_classified, value = 1)
  forest_loss <- area_change * res(j_classified)[1] * res(j_classified)[2]

  ### Add to vector
  area_vec[i] <- forest_loss

  #######################
  ### Carrasco
  #######################
  ### Clip raster to shapefile
  rast_carra <- crop(j_classified, carra)
  rast_carra <- mask(rast_carra, carra)

  ### Calculate # pixels that convert
  forest_loss_carra <- freq(rast_carra, value = 1)

  ### Calculate area that converts
  forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]

  ### Add to vector
  carra_area_vec[i] <- forest_loss_carra

  #######################
  ### Amboro NP
  #######################
  ### Clip raster to shapefile
  rast_ambnp <- crop(j_classified, amb_np)
  rast_ambnp <- mask(rast_ambnp, amb_np)

  ### Calculate # pixels that convert
  forest_loss_ambnp <- freq(rast_ambnp, value = 1)

  ### Calculate area that converts
  forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]

  ### Add to vector
  ambnp_area_vec[i] <- forest_loss_ambnp

  #######################
  ### Amboro NP
  #######################
  ### Clip raster to shapefile
  rast_ambim <- crop(j_classified, amb_im)
  rast_ambim <- mask(rast_ambim, amb_im)

  ### Calculate # pixels that convert
  forest_loss_ambim <- freq(rast_ambim, value = 1)

  ### Calculate area that converts
  forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]

  ### Add to vector
  ambim_area_vec[i] <- forest_loss_ambim
}

### Save vector of area lost
write.csv(area_vec,
          file = "/nfs/agfrontiers-data/luc_model/boliv_1_project_2018/boli_projected_loss_2018.csv",
          row.names=FALSE)
write.csv(carra_area_vec,
          file = "/nfs/agfrontiers-data/luc_model/boliv_1_project_2018/boli_projected_loss_carrasco_2018.csv",
          row.names=FALSE)
write.csv(ambnp_area_vec,
          file = "/nfs/agfrontiers-data/luc_model/boliv_1_project_2018/boli_projected_loss_amboro_np_2018.csv",
          row.names=FALSE)
write.csv(ambim_area_vec,
          file = "/nfs/agfrontiers-data/luc_model/boliv_1_project_2018/boli_projected_loss_amboro_imna_2018.csv",
          row.names=FALSE)

# #######
# ## M2
# #######
# ### Open PP2 map
# pp_2_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_m2_mask.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro IMNA
# ambim_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
# 
#   ### Make blank raster, fill with random values
#   j_blank <- random.raster(pp_2_boli_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
# 
#   ### Set CRS
#   j_crs <- crs(pp_2_boli_mask)
#   crs(j_blank) <- j_crs
# 
#   ### Set extent
#   extent(j_blank) <- extent(pp_2_boli_mask)
# 
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_2_boli_mask)
#   j_stack <- stack(j_blank, pp_2_boli_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_",
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
#   ### Carrasco
#   #######################
#   ### Clip raster to shapefile
#   rast_carra <- crop(j_classified, carra)
#   rast_carra <- mask(rast_carra, carra)
# 
#   ### Calculate # pixels that convert
#   forest_loss_carra <- freq(rast_carra, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]
# 
#   ### Add to vector
#   carra_area_vec[i] <- forest_loss_carra
# 
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambnp <- crop(j_classified, amb_np)
#   rast_ambnp <- mask(rast_ambnp, amb_np)
# 
#   ### Calculate # pixels that convert
#   forest_loss_ambnp <- freq(rast_ambnp, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]
# 
#   ### Add to vector
#   ambnp_area_vec[i] <- forest_loss_ambnp
# 
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambim <- crop(j_classified, amb_im)
#   rast_ambim <- mask(rast_ambim, amb_im)
# 
#   ### Calculate # pixels that convert
#   forest_loss_ambim <- freq(rast_ambim, value = 1)
# 
#   ### Calculate area that converts
#   forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]
# 
#   ### Add to vector
#   ambim_area_vec[i] <- forest_loss_ambim
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss.csv",
#           row.names=FALSE)
# write.csv(carra_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss_carrasco.csv",
#           row.names=FALSE)
# write.csv(ambnp_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss_amboro_np.csv",
#           row.names=FALSE)
# write.csv(ambim_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_2_project/boli_2_projected_loss_amboro_imna.csv",
#           row.names=FALSE)

#######
## M3
#######
# ### Open PP3 map
# pp_3_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_m3_mask.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro IMNA
# ambim_area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank raster, fill with random values
#   j_blank <- random.raster(pp_3_boli_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp_3_boli_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp_3_boli_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_3_boli_mask)
#   j_stack <- stack(j_blank, pp_3_boli_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_",
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
#   ### Carrasco
#   #######################
#   ### Clip raster to shapefile
#   rast_carra <- crop(j_classified, carra)
#   rast_carra <- mask(rast_carra, carra)
#   
#   ### Calculate # pixels that convert
#   forest_loss_carra <- freq(rast_carra, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_carra <- forest_loss_carra * res(rast_carra)[1] * res(rast_carra)[2]
#   
#   ### Add to vector
#   carra_area_vec[i] <- forest_loss_carra
#   
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambnp <- crop(j_classified, amb_np)
#   rast_ambnp <- mask(rast_ambnp, amb_np)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambnp <- freq(rast_ambnp, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambnp <- forest_loss_ambnp * res(rast_ambnp)[1] * res(rast_ambnp)[2]
#   
#   ### Add to vector
#   ambnp_area_vec[i] <- forest_loss_ambnp
#   
#   #######################
#   ### Amboro NP
#   #######################
#   ### Clip raster to shapefile
#   rast_ambim <- crop(j_classified, amb_im)
#   rast_ambim <- mask(rast_ambim, amb_im)
#   
#   ### Calculate # pixels that convert
#   forest_loss_ambim <- freq(rast_ambim, value = 1)
#   
#   ### Calculate area that converts
#   forest_loss_ambim <- forest_loss_ambim * res(rast_ambim)[1] * res(rast_ambim)[2]
#   
#   ### Add to vector
#   ambim_area_vec[i] <- forest_loss_ambim
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss.csv",
#           row.names=FALSE)
# write.csv(carra_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss_carrasco.csv",
#           row.names=FALSE)
# write.csv(ambnp_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss_amboro_np.csv",
#           row.names=FALSE)
# write.csv(ambim_area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/bolivia_3_project/boli_3_projected_loss_amboro_imna.csv",
#           row.names=FALSE)

#######
## M4
#######
# ### Open PP4 map
# pp_4_boli_mask <- raster("/nfs/agfrontiers-data/luc_model/boli_pp_m4_mask.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Carrasco
# carra_area_vec <- rep("", times = 1000)
# 
# ### Vector for area converted from forest to ag within Amboro NP
# ambnp_area_v