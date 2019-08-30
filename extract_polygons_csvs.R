### Load packages
library(tidyverse)
library(sf)
library(raster)

### Set wd
setwd("~/Documents/SESYNC/dinamica/remote_sensing")

### Open csv
data <- read.csv("export-agricultural-frontiers-in-the-amazon-basin.csv", stringsAsFactors = FALSE)

### Convert "-" to NA
data$wetlands[data$wetlands == "-"] <- NA

### Text to drop
drop_text <- data$soil[[1]]

### Drop from columns
data[, 13:19][data[, 13:19] == drop_text] <- NA

### Make data long
data_l <- gather(data, landclass, polygons, 
                 agriculture:wetlands, 
                 factor_key=TRUE)

# ### First five Sample IDs
# sample_ids <- data_l$Sample.id
# sample_ids <- sample_ids[1:5]
# 
# ### Filter by first five sample IDs to get smaller set
# data_subset <- data_l %>% filter(., Sample.id %in% sample_ids)

### Work with yucky long string
data_l$polygons <- gsub(".*:", "", data_l$polygons)


