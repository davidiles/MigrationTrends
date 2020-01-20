setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'jagsUI',
  
  # For parallel processing
  'parallel','doParallel',
  
  'sf','raster','sp','stars'
  
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

# Subset to Fall (for demonstration)
dat_combined <- subset(dat_combined, season == "Fall")

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
station_coordinates[station_coordinates == "NULL"] <- NA

station_coordinates <- subset(station_coordinates, station %in% dat_combined$station)
station_sf <-  st_as_sf(na.omit(station_coordinates), coords = c("lon", "lat"),crs = 4326, agr = "constant")


# Generate a plot of BCRs and associated range zones
lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS")
bcr2 <- st_transform(bcr1, crs = lcc)

map.BCR <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, fill = "gray85", col = "gray92") +
  geom_sf_text(data = bcr2, aes(label = BCR), size = 2, col = "black")

# Discrete range zones 
bcr2$zone = "South"
bcr2$zone[bcr2$PROVINCE_S %in% c("ALASKA","YUKON","BRITISH COLUMBIA","ALBERTA","SASKATCHEWAN","NORTHWEST TERRITORIES")] = "West"
bcr2$zone[bcr2$PROVINCE_S %in% c("NUNAVUT","MANITOBA","ONTARIO")] = "Central"
bcr2$zone[bcr2$zone == "South" & bcr2$COUNTRY == "CANADA"] = "East"
bcr2$zone = factor(bcr2$zone, levels = c("West","Central","East","South"))

col_pal <- RColorBrewer::brewer.pal(length(unique(bcr2$zone)),"Set2")
col_pal[4] <- "gray90"

map_BCR_zone <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, aes(fill = zone), col = "gray95")+
  scale_fill_manual(values=col_pal)
print(map_BCR_zone)

# Merge polygons based on "zone"

bcr_zone <- bcr2 %>% group_by(zone) %>% st_buffer(0) %>% summarize()

map_BCR_zone <- ggplot() +   theme_bw() +
  geom_sf(data = bcr_zone, aes(fill = zone), col = "gray95")+
  scale_fill_manual(values=col_pal)
print(map_BCR_zone)

#-------------------------------------------------------------------------------------------------------
# Part 3: Assign stations to particular "catchment areas"
#-------------------------------------------------------------------------------------------------------

# Step 1: subset to only a few relevant stations (temporary)
siteloc_map <- map_BCR_zone +
  #geom_sf(data = station_sf, size = 2, shape = 1) +
  geom_sf_text(data = station_sf, aes(label = station), size = 2, col = "black") +
  ggtitle("Sites with fall counts")
print(siteloc_map)

west_fall <- data.frame(station = c("CFMS",  # Creamer's field migration station
                                    "TLBBS", # Teslin Lake
                                    "MNO",   # Mackenzie nature observatory
                                    "LMBO"),  # Last Mountain Bird Observatory
                        zone = "West")
central_fall <- data.frame(station = c("MCCS", # Manomet
                                       "KWRS", # Kingston
                                       "BIBS"), # Block island
                           zone = "Central")
east_fall <- data.frame(station = c("MGBO", # McGill
                                    "PEPBO"), # Prince Edward Point
                        zone = "East")

station_zones <- rbind(west_fall,central_fall,east_fall)
station_sf <- full_join(station_sf, station_zones)

analysis_site_map <- map_BCR_zone +
  #geom_sf(data = station_sf, size = 2, shape = 1) +
  geom_sf_text(data = na.omit(station_sf), 
               aes(label = station), size = 2, col = "black") +
  ggtitle("Sites with fall counts")

print(analysis_site_map)

#-------------------------------------------------------------------------------------------------------
# Part 4: Assign relative abundances to each BCR
#-------------------------------------------------------------------------------------------------------

# Load BAM relative density raster data
bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")

# crop the raster
bam1_crop <- mask(bam1, bcr_zone) #takes a long time


# calculate mean density in each BCR
bcr_summary <-
  bcr_zone %>% mutate(
    mean_dens = raster_extract(bam1_crop, bcr_zone, fun = mean, na.rm = TRUE),
    n.pixel = raster_extract(bam1_crop, bcr_zone, fun = function(x) sum(is.na(x)), na.rm = TRUE)
  )

test = extract(bam1_crop, as(bcr_zone, 'Spatial'), fun = function(x) mean(x, na.rm = TRUE))
