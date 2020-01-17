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
# Define discrete range zones
#------------------------------

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


