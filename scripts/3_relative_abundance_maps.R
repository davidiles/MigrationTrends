setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'sp','raster','sf','rgdal'
  
)

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 1: READ GEOTIFF FILE FROM BAM
#******************************************************************************************************************************************
#******************************************************************************************************************************************
GDALinfo("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")

bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")
plot(bam1, col = bpy.colors(n = 100), main = "BLPW predicted density")

# Specify a quantile-based map
qs = unique(round(quantile(bam1, c(seq(0,1,length.out = 50))),2))
plot(bam1, breaks = qs, col = viridis(length(qs)), main = "BLPW predicted density")

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 2: READ STABLE ISOTOPE BOUNDARIES FILE
#******************************************************************************************************************************************
#******************************************************************************************************************************************

bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY == "CANADA")

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 3: REPROJECT BOUNDARIES IN SAME CRS AS RELATIVE DENSITY MAP
#******************************************************************************************************************************************
#******************************************************************************************************************************************
bcr2 <- st_transform(bcr1, crs = crs(bam1))
bcr2.sp <- as(bcr2, "Spatial") # convert to spatial points

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 4: PLOT DENSITY MAP AND BCR BOUNDARIES
#******************************************************************************************************************************************
#******************************************************************************************************************************************

# Plot the map

bam1.mod <- bam1
#bam1.mod[bam1.mod > quantile(bam1.mod,0.99)] <- quantile(bam1.mod,0.99)

par(mar=c(5,5,5,5))

# Equal spacing for density ramp
brks = round(seq(0,cellStats(bam1.mod, max),length.out = 9),2)
brks[length(brks)] = brks[length(brks)] + 0.01

# Breaks placed at relevant quantiles
brks <- unique(round(quantile(bam1.mod, c(seq(0,1,length.out = 9))),3))
brks[length(brks)] = brks[length(brks)] + 0.1


col.pal <- c("gray85",RColorBrewer::brewer.pal(9,"Reds"))
col.pal <- colorRampPalette(c("gray85",RColorBrewer::brewer.pal(9,"Reds")[7]))(9)
#col.pal <- rev(viridis(length(brks)))
#col.pal <- rev(colorRampPalette(col.pal)(length(brks)))

plot(bcr2.sp, lwd=1, border = "white", col = "gray85") 

plot(bam1.mod, breaks = brks, col = col.pal, 
     add = TRUE,
     main = "BLPW predicted density",
     npretty = 2,
     ylim = extent(bam1.mod)[3:4],
     xlim = extent(bam1.mod)[1:2])

# Overlay region boundaries
plot(bcr2.sp, add=TRUE, lwd=1, border = "white") 

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 5: EXTRACT MEAN DENSITY IN EACH BCR
#******************************************************************************************************************************************
#******************************************************************************************************************************************

start_time <- Sys.time()

region.means <- extract(bam1, bcr2.sp, fun = mean, na.rm = TRUE)
region.sums <- extract(bam1, bcr2.sp, fun = sum, na.rm = TRUE)

end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time) # approx 16 minutes on ECCC laptop

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 6: PLOT SUMMED DENSITY IN EACH BCR
#******************************************************************************************************************************************
#******************************************************************************************************************************************
bcr2.sp$BLPW.SUM <- as.vector(region.sums)
bcr3.sp <- st_as_sf(bcr2.sp)

plot(bcr3.sp["BLPW.SUM"], key.pos = 4)
