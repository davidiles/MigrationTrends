setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'jagsUI',
  
  # For parallel processing
  'parallel','doParallel'
  
)

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())


#*******************************************************************************************************
# Model components:
#*******************************************************************************************************
#
# 1) Daily counts at stations
#
# 2) Locations of individual stations (for plotting and determining likely catchment areas)
#
# 3) "Sampling regions" or "catchment areas" for each station (can include mixtures)
#
# 4) Relative abundances in each BCR (which can be rolled up into larger jurisdictions/regions)
#
#*******************************************************************************************************


#-------------------------------------------------------------------------------------------------------
# Part 1: Load daily count data (pre-processed in file dat_combined.csv by script "2_format_data_analyze_separately.R")
#-------------------------------------------------------------------------------------------------------
dat_combined <- read.csv("./processed_data/dat_combined.csv")

#-------------------------------------------------------------------------------------------------------
# Part 2: Load locations of individual stations
#-------------------------------------------------------------------------------------------------------

cmmn_coordinates <- read.csv("../data/locations/CMMN_locations.csv")
colnames(cmmn_coordinates) <- c("results_code","station","name","statprov_code","lat","lon")
cmmn_coordinates$station <- as.character(cmmn_coordinates$station)

# Coordinates of US migration stations
us_coordinates <- rbind(data.frame(station = "MCCS", lat = 41.91, lon = -70.54),
                             data.frame(station = "AIMS", lat = 42.99, lon = -70.61),
                             data.frame(station = "KWRS", lat = 41.48, lon = -71.53),
                             data.frame(station = "BIBS", lat = 41.21, lon = -71.58),
                             data.frame(station = "BSBO", lat = 41.61, lon = -83.19),
                             data.frame(station = "FBBO", lat = 39.20, lon = -76.06),
                             data.frame(station = "PARC", lat = 40.16, lon = -79.27),
                             data.frame(station = "CFMS", lat = 64.86, lon = -147.74)
)

station_coordinates <- rbind(us_coordinates, cmmn_coordinates[,c("station","lat","lon")])

#--------------------
# Generate a plot of BCRs and associated range zones
#--------------------

lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS")
bcr2 <- st_transform(bcr1, crs = lcc)

map.BCR <- ggplot() +   theme_bw() +
  
  geom_sf(data = bcr2, fill = "gray85", col = "gray92") +
  geom_sf_text(data = bcr2, aes(label = BCR), size = 2, col = "black")

# Discrete range zones 
bcr2$range = "South"
bcr2$range[bcr2$PROVINCE_S %in% c("ALASKA","YUKON","BRITISH COLUMBIA","ALBERTA","SASKATCHEWAN","NORTHWEST TERRITORIES")] = "West"
bcr2$range[bcr2$PROVINCE_S %in% c("NUNAVUT","MANITOBA","ONTARIO")] = "Central"
bcr2$range[bcr2$range == "South" & bcr2$COUNTRY == "CANADA"] = "East"
bcr2$range = factor(bcr2$range, levels = c("West","Central","East","South"))

col.pal <- RColorBrewer::brewer.pal(length(unique(bcr2$range)),"Set2")
col.pal[4] <- "gray90"

map.BCR.range <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, aes(fill = range), col = "gray95")+
  scale_fill_manual(values=col.pal)

map.BCR.range

#-------------------------------------------------------------------------------------------------------
# Part 3: Assign stations to particular "catchment areas"
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# Part 4: Assign relative abundances to each BCR
#-------------------------------------------------------------------------------------------------------

# Load BAM relative density raster data
bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")
