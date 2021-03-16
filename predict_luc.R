### Code to make LUC models and apply to 2018 LUC maps

### Open packages
library(tidyverse)
library(sf)
library(raster)
library(spatialEco)

##################
### Brazil model 1 
##################

### Data prep
#############

# ### Open sample points
# bra_dat_all <- read.csv("/nfs/agfrontiers-data/luc_model/brazil_data_for_model.csv")
# 
# ### Rename columns
# colnames(bra_dat_all) <- c("uid", "lon", "lat", "lc_2008", "lc_2018",
#                         "id", "aspect", "mines",
#                         "roads", "rivers",
#                         "cities", "elev",
#                         "popdens", "poverty",
#                         "ppt", "prot_status",
#                         "slope", "soil",
#                         "crop_suit", "field_res",
#                         "paddd", "non_all_land",
#                         "dist_ill_mines", "dist_ag",
#                         "dist_fires", "fire_dens",
#                         "dist_pr_rr",
#                         "dist_pr_dams",
#                         "agref_sett", "cattle",
#                         "perc_diff")

# ### Make column for transition/none
# bra_dat_all <- bra_dat_all %>%
#   mutate(transition = ifelse(lc_2018 == "2",
#                              "no_change",
#                              ifelse(lc_2018 == "1",
#                                     "f_to_ag",
#                                     "other_change")),
#          trans_rc = ifelse(transition == "no_change",
#                            "0", ifelse(transition == "f_to_ag",
#                                        "1", "2")))

# ### Make columns factor
# fact_cols <- c("uid", "prot_status",
#                "paddd", "agref_sett",
#                "lc_2018", "trans_rc")
# bra_dat_all <- bra_dat_all %>%
#   mutate_each_(funs(factor(.)), fact_cols)
# 
# ### Restrict to transitions
# bra_dat_use <- bra_dat_all %>%
#   filter(trans_rc == "1" | trans_rc == "0")

### Model 1
###########

# ### Run regression
# fit_bra <- glm(trans_rc ~ aspect + slope + elev +
#                  roads + rivers + mines + cities +
#                  crop_suit + popdens + #poverty +
#                  # (1|uid) + lon +
#                  soil + perc_diff + prot_status,
#                data = bra_dat_use,
#                family = binomial())

### 2018 data
#############

# ### Open rasterbrick
# bra_layers_all <- brick("/nfs/agfrontiers-data/luc_model/brazil_all_layers.tif")
# 
# ### Rename columns
# names(bra_layers_all) <- c("aspect", "mines",
#                            "roads", "rivers",
#                            "cities", "elev",
#                            "popdens", "poverty",
#                            "ppt", "prot_status",
#                            "slope", "soil",
#                            "crop_suit", "field_res",
#                            "paddd", "non_all_land",
#                            "dist_ill_mines", "dist_ag",
#                            "dist_fires", "fire_dens",
#                            "dist_pr_rr",
#                            "dist_pr_dams",
#                            "agref_sett", "cattle")

# ### Open 2018 percent of neighboring pixels raster
# neigh_2018 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/bra18_perc_diff.tif")
# 
# ### Add to other layers
# bra_18_all <- addLayer(bra_layers_all, neigh_2018)
# names(bra_18_all)[25] <- "perc_diff"
# # writeRaster(bra_18_all,
# #             filename="/nfs/agfrontiers-data/Remote Sensing/KS files/bra18_all_data_w_perc_diff.tif",
# #             options="INTERLEAVE=BAND", overwrite=TRUE)
# 
# ### Convert NAs to 0 in ag reform layer
# bra_18_all[[23]][is.na(bra_18_all[[23]][])] <- 0
# 
# ### Drop precip and soil layers because of issue with extent
# bra_18_all <- dropLayer(bra_18_all, c(9, 12))
# 
# ### Open fixed soil layer
# soil <- raster("/nfs/agfrontiers-data/luc_model/sm_jaman_102033.tif")
# 
# ### Add to stack
# bra_18_all <- addLayer(bra_18_all, soil)
# names(bra_18_all)[24] <- "soil"
# 
# ### Open distance to ag in 2018 layer
# d_ag_18 <- raster("/nfs/agfrontiers-data/luc_model/bra_ag_dist_re_2018.tif")
# 
# ### Drop 2008 distance to ag, replace with 2018 distance to ag
# bra_18_all <- dropLayer(bra_18_all, 16)
# bra_18_all <- addLayer(bra_18_all, d_ag_18)
# names(bra_18_all)[24] <- "dist_ag"

### Run prediction
##################

# ### Run prediction
# pp_brazil <- raster::predict(bra_18_all, fit_bra,
#                              na.rm = TRUE,
#                              type = "response",
#                              filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_m1_logit.tif",
#                              overwrite = TRUE)

### Open 2018 LUC map
# brazi18 <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bra_dry_2018_102033.tif")

# ### Mask cells that were not forested in 2018
# pp_brazil_mask <- mask(pp_brazil, brazi18,
#                        # filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_m1_mask.tif",
#                        inverse= TRUE,
#                        maskvalue = 2)
# 
# writeRaster(pp_brazil_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_masked.tif",
#                         format = "GTiff",
#                         overwrite = TRUE,
#                         options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

##################
### Brazil model 2
##################
# ### Run regression
# fit_bra_2 <- glm(trans_rc ~ paddd + non_all_land +
#                    # cattle +
#                     # field_res +
#                     dist_ill_mines +
#                     dist_ag + dist_fires +
#                     fire_dens + dist_pr_rr +
#                     # dist_pr_dams +
#                     agref_sett,
#                   data = bra_dat_use,
#                   family = binomial(link = "logit"))

### Run prediction
##################
# pp_brazil_2 <- raster::predict(bra_18_all, fit_bra_2,
#                              na.rm = TRUE,
#                              type = "response")
#                              # filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_m2_logit.tif",
#                              # overwrite = TRUE)
# 
# ### Mask cells that were not forested in 2018
# pp_brazil_2_mask <- mask(pp_brazil_2, brazi18,
#                        inverse= TRUE,
#                        maskvalue = 2)
# 
# writeRaster(pp_brazil_2_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_2_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

##################
### Brazil model 3 
##################
# ### Run regression
# fit_bra_3 <- glm(trans_rc ~ paddd + non_all_land +
#                    # field_res +
#                    dist_ill_mines +
#                    dist_ag + dist_fires +
#                    fire_dens + dist_pr_rr +
#                    # dist_pr_dams +
#                    agref_sett + #cattle +
#                    aspect + slope + elev +
#                    roads + rivers + mines + cities +
#                    crop_suit + popdens + #poverty +
#                    # (1|uid) + lon +
#                    soil + perc_diff + prot_status,
#                  data = bra_dat_use,
#                  family = binomial(link = "logit"))
# 
# ### Run prediction
# ##################
# pp_brazil_3 <- raster::predict(bra_18_all, fit_bra_3,
#                                na.rm = TRUE,
#                                type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_brazil_3_mask <- mask(pp_brazil_3, brazi18,
#                          inverse= TRUE,
#                          maskvalue = 2)
# 
# writeRaster(pp_brazil_3_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_3_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

##################
### Brazil model 4 
##################
# ### Run regression
# fit_bra_4 <- glm(trans_rc ~ slope + elev + 
#                    roads + cities + mines +
#                    crop_suit + soil + perc_diff + prot_status +
#                    paddd + non_all_land +
#                    dist_ag + dist_fires +
#                    fire_dens + dist_pr_rr +
#                    agref_sett,
#                  data = bra_dat_use,
#                  family = binomial(link = "logit"))

### Run prediction
##################
# pp_brazil_4 <- raster::predict(bra_18_all, fit_bra_4,
#                                na.rm = TRUE,
#                                type = "response")
# 
# ### Mask cells that were not forested in 2018
# pp_brazil_4_mask <- mask(pp_brazil_4, brazi18,
#                          inverse= TRUE,
#                          maskvalue = 2)
# 
# writeRaster(pp_brazil_4_mask,
#             filename = "/nfs/agfrontiers-data/luc_model/brazil_pp_4_masked.tif",
#             format = "GTiff",
#             overwrite = TRUE,
#             options = c("INTERLEAVE=BAND","COMPRESS=LZW"))


######################
### Save model outputs 
######################

# # ### Model 1
# fit_bra_summ <- summary(fit_bra)
# fit_bra_coeff <- as.data.frame(fit_bra_summ$coefficients)
# write.csv(fit_bra_coeff,
#           file = "/nfs/agfrontiers-data/luc_model/fit_bra_model1.csv")
# 
# ### Model 2
# fit_bra_summ2 <- summary(fit_bra_2)
# fit_bra_coeff2 <- as.data.frame(fit_bra_summ2$coefficients)
# write.csv(fit_bra_coeff2,
#           file = "/nfs/agfrontiers-data/luc_model/fit_bra_model2.csv")
# 
# ### Model 3
# fit_bra_summ3 <- summary(fit_bra_3)
# fit_bra_coeff3 <- as.data.frame(fit_bra_summ3$coefficients)
# write.csv(fit_bra_coeff3,
#           file = "/nfs/agfrontiers-data/luc_model/fit_bra_model3.csv")

# ### Model 4
# fit_bra_summ4 <- summary(fit_bra_4)
# fit_bra_coeff4 <- as.data.frame(fit_bra_summ4$coefficients)
# write.csv(fit_bra_coeff4,
#           file = "/nfs/agfrontiers-data/luc_model/fit_bra_model4.csv")

##################
### Compare models 
##################

# ### anova
# anova(fit_bra, fit_bra_2, fit_bra_3, fit_bra_4, test ="Chisq")
# anova(fit_bra_3, fit_bra_4, test ="Chisq")
# 
# ### psuedo r2
# library(pscl)
# pR2(fit_bra)
# pR2(fit_bra_2)
# pR2(fit_bra_3)
# pR2(fit_bra_4)

##############################################################################
##############################################################################
##############################################################################

##################
### MC simulations 
##################

######################
### Model 1 projection
######################

# ### Open PP1 map
# pp_brazil_mask <- raster("/nfs/agfrontiers-data/luc_model/brazil_pp_masked.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
# 
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp_brazil_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
# 
#   ### Set CRS
#   j_crs <- crs(pp_brazil_mask)
#   crs(j_blank) <- j_crs
# 
#   ### Set extent
#   extent(j_blank) <- extent(pp_brazil_mask)
# 
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp_brazil_mask)
#   j_stack <- stack(j_blank, pp_brazil_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/brazil_1_project/brazil_projected_",
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
#   ### Calculate area that converts
#   area_change <- as.data.frame(j_classified) %>%
#     group_by(layer) %>%
#     tally() %>%
#     mutate(area = n * res(j_classified)[1] * res(j_classified)[2],
#            iteration = i)
# 
#   ### Subset to area that converted
#   forest_loss <- area_change[2, 3]
# 
#   ### Add to vector
#   area_vec[i] <- forest_loss
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/brazil_1_project/brazil_projected_loss.csv",
#           row.names=FALSE)

######################
### Model 2 projection
######################
# ### Open PP2 map
# pp2_brazil_mask <- raster("/nfs/agfrontiers-data/luc_model/brazil_pp_2_masked.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
# 
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp2_brazil_mask,
#                            min = 0,
#                            max = 1,
#                            distribution = "random")
# 
#   ### Set CRS
#   j_crs <- crs(pp2_brazil_mask)
#   crs(j_blank) <- j_crs
# 
#   ### Set extent
#   extent(j_blank) <- extent(pp2_brazil_mask)
# 
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp2_brazil_mask)
#   j_stack <- stack(j_blank, pp2_brazil_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/brazil_2_project/brazil_2_projected_",
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
# 
#   ### Calculate area that converts
#   area_change <- as.data.frame(j_classified) %>%
#     group_by(layer) %>%
#     tally() %>%
#     mutate(area = n * res(j_classified)[1] * res(j_classified)[2],
#            iteration = i)
# 
#   ### Subset to area that converted
#   forest_loss <- area_change[2, 3]
# 
#   ### Add to vector
#   area_vec[i] <- forest_loss
# }
# 
# ### Save vector of area lost
# write.csv(area_vec,
#           file = "/nfs/agfrontiers-data/luc_model/brazil_2_project/brazil_2_projected_loss.csv",
#           row.names=FALSE)

######################
### Model 3 projection
######################
# ### Open PP3 map
# pp3_brazil_mask <- raster("/nfs/agfrontiers-data/luc_model/brazil_pp_3_masked.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp3_brazil_mask,
#                            min = 0, 
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp3_brazil_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp3_brazil_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp3_brazil_mask)
#   j_stack <- stack(j_blank, pp3_brazil_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/brazil_3_project/brazil_3_projected_", 
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
#   
#   ### Calculate area that converts
#   area_change <- as.data.frame(j_classified) %>%
#     group_by(layer) %>%
#     tally() %>%
#     mutate(area = n * res(j_classified)[1] * res(j_classified)[2],
#            iteration = i)
#   
#   ### Subset to area that converted
#   forest_loss <- area_change[2, 3]
#   
#   ### Add to vector
#   area_vec[i] <- forest_loss
# }
# 
# ### Save vector of area lost
# write.csv(area_vec, 
#           file = "/nfs/agfrontiers-data/luc_model/brazil_3_project/brazil_3_projected_loss.csv", 
#           row.names=FALSE)


######################
### Model 4 projection
######################
# ### Open PP4 map
# pp4_brazil_mask <- raster("/nfs/agfrontiers-data/luc_model/brazil_pp_4_masked.tif")
# 
# ### Vector for area converted from forest to ag
# area_vec <- rep("", times = 1000)
# 
# for (i in 1:length(area_vec)) {
#   
#   ### Make blank Jamanxim raster, fill with random values
#   j_blank <- random.raster(pp4_brazil_mask,
#                            min = 0, 
#                            max = 1,
#                            distribution = "random")
#   
#   ### Set CRS
#   j_crs <- crs(pp4_brazil_mask)
#   crs(j_blank) <- j_crs
#   
#   ### Set extent
#   extent(j_blank) <- extent(pp4_brazil_mask)
#   
#   ### Combine with model 1 PP map
#   j_blank <- crop(j_blank, pp4_brazil_mask)
#   j_stack <- stack(j_blank, pp4_brazil_mask)
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
#               filename = paste0("/nfs/agfrontiers-data/luc_model/brazil_4_project/brazil_4_projected_", 
#                                 i, ".tif"),
#               format = "GTiff",
#               overwrite = TRUE,
#               options = c("INTERLEAVE=BAND","COMPRESS=LZW"))
#   
#   ### Calculate area that converts
#   area_change <- as.data.frame(j_classified) %>%
#     group_by(layer) %>%
#     tally() %>%
#     mutate(area = n * res(j_classified)[1] * res(j_classified)[2],
#            iteration = i)
#   
#   ### Subset to area that converted
#   forest_loss <- area_change[2, 3]
#   
#   ### Add to vector
#   area_vec[i] <- forest_loss
# }
# 
# ### Save vector of area lost
# write.csv(area_vec, 
#           file = "/nfs/agfrontiers-data/luc_model/brazil_4_project/brazil_4_projected_loss.csv", 
#           row.names=FALSE)

####################################################################################
####################################################################################

### Loop through simulated maps and: 
### 1. make raster that adds values of all map layers
### 2. calculate forest loss within Jamanxim NF and put that value in a vector
### 3. calculate FRAGSTATS and put in dataframe -- can't allocate vector...

# ### Add library for FRAGSTATS
# library(landscapemetrics)

### Open Jamanxim without buffer
jaman <- st_read("/nfs/agfrontiers-data/Case Study Info/Brazil/JamanBuffer/jamanxim_wdpa_march21_102033.shp")

### Try on subset of layers in Brazil Model 1
test_layers <- list.files(path = "/nfs/agfrontiers-data/luc_model/brazil_1_project", 
                          pattern = "*.tif", 
                          full.names = TRUE)

# test_layers <- test_layers[1:10]

### Vector for area converted from forest to ag within Jaman
jam_area_vec <- rep("", times = 1000)

### Blank raster for adding all rasters
raster_fill <- raster(test_layers[1])
raster_fill[raster_fill > 0] <- 0

### Write loop
for (i in 1:length(test_layers)) {
  
  ### Open raster
  ch_rast <- raster(test_layers[i])
  
  ###################################
  ### Add raster values
  ###################################
  ### Stack with raster_fill
  raster_fill <- stack(raster_fill, ch_rast)
  
  ### Add raster values
  raster_fill <- calc(raster_fill, sum)
  
  ###################################
  ### Calculate forest loss within JNP
  ###################################  
  ### Clip raster to shapefile
  rast_crop <- crop(ch_rast, jaman)
  rast_crop <- mask(rast_crop, jaman)
  
  ### Calculate # pixels that convert
  forest_loss <- freq(rast_crop, value = 1)
  
  ### Calculate area that converts
  forest_loss <- forest_loss * res(rast_crop)[1] * res(rast_crop)[2]
  
  ### Add to vector
  jam_area_vec[i] <- forest_loss
  
}

### Write out sum raster
writeRaster(raster_fill,
            filename = "/nfs/agfrontiers-data/luc_model/brazil_m1_stacked.tif",
            format = "GTiff",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"))

### Save vector of area lost
write.csv(jam_area_vec, 
          file = "/nfs/agfrontiers-data/luc_model/brazil_m1_jamanxim_forestloss.csv", 
          row.names=FALSE)