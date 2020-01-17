# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'jagsUI',
  
  # For parallel processing
  'parallel','doParallel',
  
  # For analysis
  'sp','raster','sf','rgdal'
  
)

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS")
bcr2 <- st_transform(bcr1, crs = lcc)

map.BCR <- ggplot() +   theme_bw() +
  
  geom_sf(data = bcr2, fill = "gray85", col = "gray92") +
  geom_sf_text(data = bcr2, aes(label = BCR), size = 2, col = "black")

#------------------------------
# Define "western range"
#------------------------------

bcr2$range = 0
bcr2$range[bcr2$PROVINCE_S %in% c("ALASKA","YUKON","BRITISH COLUMBIA","ALBERTA","SASKATCHEWAN","NORTHWEST TERRITORIES")] = "West"

map.BCR.west <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, aes(fill = factor(origin.west)), col = "gray92")

map.BCR.west

#------------------------------
# Define "central range"
#------------------------------

bcr2$origin.central = 0
bcr2$origin.central[bcr2$PROVINCE_S %in% c("NUNAVUT","MANITOBA","ONTARIO")] = 1

map.BCR.central <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, aes(fill = factor(origin.central)), col = "gray92")

map.BCR.central

#------------------------------
# Define "eastern range"
#------------------------------

bcr2$origin.eastern = 0
bcr2$origin.eastern[bcr2$origin.west == 0 & bcr2$origin.central == 0 & bcr2$COUNTRY == "CANADA"] = 1

map.BCR.eastern <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, aes(fill = factor(origin.eastern)), col = "gray92")

map.BCR.eastern


