setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'jagsUI',
  
  # For parallel processing
  'parallel','doParallel',
  
  'sf','raster','sp',
  'maptools','rgeos'
  
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
# Part 2: Define discrete geographic regions for stations
#-------------------------------------------------------------------------------------------------------

# Generate a plot of BCRs and associated range regions
lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS")

# Temporary subset for script testing
bcr1 <- subset(bcr1, COUNTRY == "CANADA")

bcr2 <- st_transform(bcr1, crs = lcc)

# map.BCR <- ggplot() +   theme_bw() +
#   geom_sf(data = bcr2, fill = "gray85", col = "gray92") +
#   geom_sf_text(data = bcr2, aes(label = BCR), size = 2, col = "black")

# Discrete range regions 
bcr2$region = "South"
bcr2$region[bcr2$PROVINCE_S %in% c("ALASKA","YUKON","BRITISH COLUMBIA","ALBERTA","SASKATCHEWAN","NORTHWEST TERRITORIES")] = "West"
bcr2$region[bcr2$PROVINCE_S %in% c("NUNAVUT","MANITOBA","ONTARIO")] = "Central"
bcr2$region[bcr2$region == "South" & bcr2$COUNTRY == "CANADA"] = "East"
bcr2$region = factor(bcr2$region, levels = c("West","Central","East","South"))

col_pal <- RColorBrewer::brewer.pal(length(unique(bcr2$region)),"Set2")
col_pal[4] <- "gray90"

map_BCR_region <- ggplot() +   theme_bw() +
  geom_sf(data = bcr2, aes(fill = region), col = "gray95")+
  scale_fill_manual(values=col_pal)
print(map_BCR_region)

#--------------------------------------------------------------------
# Merge polygons based on "region"
bcr_region <- bcr2 %>% group_by(region) %>% st_buffer(0) %>% summarize()
bcr2.sp <- as(bcr2, "Spatial")

# simplify polygons
bcr2.sp <- gSimplify(bcr2.sp, tol = 0.00001)
bcr2.sp <- gBuffer(bcr2.sp, byid=TRUE, width=0) #further simplify
region_union <- unionSpatialPolygons(bcr2.sp, bcr2$region)
#--------------------------------------------------------------------

# map_BCR_region <- ggplot() +   theme_bw() +
#   geom_sf(data = bcr_region, aes(fill = region), col = "gray95")+
#   scale_fill_manual(values=col_pal)
# print(map_BCR_region)


#-------------------------------------------------------------------------------------------------------
# Part 3: Assign relative abundances to each BCR
#-------------------------------------------------------------------------------------------------------

# Load BAM relative density raster data
bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")

# Merge polygons based on BCR
#bcr_merged <- bcr2 %>% group_by(BCR) %>% st_buffer(0) %>% summarize()

# Convert to sp object instead of sf
#bcr3 <- as(bcr_merged, "Spatial")
#bcr3 <- spTransform(bcr2.sp, crs(bam1))

# crop the raster
bam1_crop <- crop(bam1, region_union) 
bam1_mask <- mask(bam1_crop, region_union) #takes a long time

plot(bam1_mask)
plot(region_union, add = TRUE)

summary(bam1_mask)

# Extract pixels in each polygon; store in list
pixels = raster::extract(bam1_mask, region_union)
region_union$means <- unlist(lapply(pixels, mean, na.rm = TRUE))
region_union$npixel <- unlist(lapply(pixels, function(x) sum(!is.na(x)))) 
region_union$z <- 1:length(pixels)
region_union$abund <- region_union$means * region_union$npixel

# Convert to sf object for plotting
region_union_sf <- st_as_sf(region_union)

abund.region <- ggplot() +   theme_bw() +
  geom_sf(data = region_union_sf, aes(fill = abund))+
  geom_sf_text(data = region_union_sf, aes(label = z, col = factor(z)), size = 4)+
  scale_fill_gradientn(colors = c("white","black"), limits = c(0,max(region_union_sf$abund)))

plot(abund.region)

#-------------------------------------------------------------------------------------------------------
# Part 4: Assign stations to particular "catchment areas"
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




west_fall <- data.frame(station = c("CFMS"  # Creamer's field migration station
                                    #"TLBBS"), # Teslin Lake
                                    #"MNO",   # Mackenzie nature observatory
                                    #"LMBO"),  # Last Mountain Bird Observatory
),
                        region = "West")

central_fall <- data.frame(station = c("MCCS"), # Manomet
                                       #"KWRS", # Kingston
                                       #"BIBS"), # Block island
                           
                           region = "Central")

east_fall <- data.frame(station = c("MGBO"), # McGill
                                    #"PEPBO"), # Prince Edward Point
                        region = "East")

station_regions <- rbind(west_fall,central_fall,east_fall)
station_sf <- full_join(station_sf, station_regions)
station_sf <- na.omit(station_sf)


analysis_site_map <- map_BCR_region +
  #geom_sf(data = station_sf, size = 2, shape = 1) +
  geom_sf_text(data = na.omit(station_sf),
               aes(label = station), size = 2, col = "black") +
  ggtitle("Sites with fall counts")

print(analysis_site_map)

#-------------------------------------------------------------------------------------------------------
# Part 5: Package data for analysis
#-------------------------------------------------------------------------------------------------------
region_abund_df <- data.frame(region_number = 1:3,
                       region = c("Central","East","West"),
                       abund = region_union_sf$abund)

location_df <- na.omit(station_sf)
location_df$station_number <- 1:nrow(location_df) # assign numbers to each station
location_df <- merge(location_df, region_abund_df, all = TRUE)

# Merge region and station numbers into count dataframe
count_df <- subset(dat_combined, station %in% location_df$station & season == "Fall")
count_df <- merge(count_df, as.data.frame(location_df)[,c("station","region","station_number")], all = TRUE)
count_df <- merge(count_df, region_abund_df[,c("region_number","region")])

# Number of areas per station
narea = aggregate(area ~ station_number, data = count_df, FUN = function(x) length(unique(x)))
narea = narea[order(narea$station_number),]

# Relative date variables
count_df$doy_adjusted = count_df$doy - min(count_df$doy) + 1
count_df$year_adjusted = count_df$YearCollected - min(count_df$YearCollected) + 1
count_df = subset(count_df, !is.na(ObservationCount))

#-----------------------------
# Construct hypothetical catchment data for each station (25 birds)
# Number of individuals arriving at each station, in each year, from each region
N.isotope <- array(NA, dim = c(region = max(region_abund_df$region_number),station = max(location_df$station_number),year = max(count_df$year_adjusted)))
for (s in location_df$station_number){
    N.isotope[,s,5] <- rep(0,max(region_abund_df$region_number)) # zero counts
    N.isotope[location_df$region_number[which(location_df$station_number == s)],s,5] <- 25
}

# Matrix to contain number of birds sampled in each year at each station (fills in numbers for years with no data)
N.station.sampled <- array(25,dim = c(station = max(location_df$station_number),year = max(count_df$year_adjusted)))
#-----------------------------

jags.data = list(
  
  # Regional regions and initial abundances
  nregion = length(unique(count_df$region)),
  
  logN0 = rep(log(1),3), #log(region_abund_df$abund),
  
  # Number of stations in dataset
  nstation = length(unique(count_df$station_number)),
  
  # Number of "areas" per station (1 for all stations except LPBO)
  narea = narea$area,
  
  # Number of years in dataset
  nyear = max(count_df$year_adjusted),
  
  # Number of days in dataset
  nday = max(count_df$doy_adjusted),
  
  # Isotope/catchment data
  N.isotope = N.isotope,
  N.station.sampled = N.station.sampled,
  
  # Number of observations in total dataset
  nobs = nrow(count_df),
  
  # Daily counts at each station/area
  daily.count = count_df$ObservationCount,
  
  # Observation-level data
  day = count_df$doy_adjusted,
  station = count_df$station_number,
  year = count_df$year_adjusted,

  log.offset = log(count_df$net.hrs),
  pi = pi
)


sink("cmmn_national.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics in each region
      #---------------------------------------------

      # No uncertainty in initial abundance; logN0[r] is fed in as data
        
      for (r in 1:nregion){
        trend[r] ~ dnorm(0,1) # Regional trends
      }
      
      # Temporal variance in trend
      proc.sd ~ dunif(0,2)
      proc.tau <- pow(proc.sd,-2)
      
      # True (unobserved) population dynamics in each region
      for (y in 1:nyear){
        for (r in 1:nregion){
        
          # Exponential population model
          logN[r,y] <- logN0[r] + trend[r] * (y-1) + noise[r,y]
          noise[r,y] ~ dnorm(0,proc.tau)
          N[r,y] <- exp(logN[r,y])
        }
      }
      
      #---------------------------------------------
      # Observed numbers of birds at station [s] originating from region [r] in year [y]
      #---------------------------------------------
      
      # Proportion of birds from region [r] that are 
      # captured by station [s] (constant through time)
      for (r in 1:nregion){
        for (s in 1:nstation){
          rho[r,s] ~ dunif(0,1) 
        }
      }
      
      for (s in 1:nstation){
        for (y in 1:nyear){
          for (r in 1:nregion){
          
            # Seasonal total arriving at station [s] from region [r]
            mu[r,s,y] <- N[r,y] * rho[r,s]    
          
          }
        
          # Total seasonal count at station [s]
          true.total[s,y] <- sum(mu[1:nregion,s,y]) 
          
          # Multinomial to describe regional composition in this year, at this station
          N.isotope[1:nregion,s,y] ~ dmulti(mu[,s,y], N.station.sampled[s,y])
          
        }
      }
      
      #---------------------------------------------
      # Estimate totals at each station each year
      #---------------------------------------------
 
      for (s in 1:nstation){
      
      mean.migrate.HYPERMEAN[s] ~ dunif(1,nday)
      mean.migrate.HYPERSD[s] ~ dunif(0,nday/2)
      mean.migrate.HYPERTAU[s] <- pow(mean.migrate.HYPERSD[s], -2)
      
      # Parameter describing the width of the migration window (assume this window is constant)
      sd.migrate[s] ~ dunif(0,nday/2)
        
      daily.noise.sd[s] ~ dunif(0,2)
      daily.noise.tau[s] <- pow(daily.noise.sd[s],-2)
      
        for (y in 1:nyear){
        
        mean.migrate[s,y] ~ dnorm(mean.migrate.HYPERMEAN[s], mean.migrate.HYPERTAU[s])
        
          for (d in 1:nday){
            
            norm.density[d,s,y] <- 1/(sqrt(2*pi)*sd.migrate[s])*exp(-((d-mean.migrate[s,y])^2/(2*sd.migrate[s]^2)))
            
            # Expected count on each day
            expected.count[d,s,y] <- norm.density[d,s,y] * true.total[s,y]
            
            # Daily observation error
            daily.noise[d,s,y] ~ dnorm(0,daily.noise.tau[s])
            
            log.lambda[d,s,y] <- log(expected.count[d,s,y]) + daily.noise[d,s,y]

          }
        }
      }
      
      for (i in 1:nobs){
        log.lam[i] <- log.lambda[day[i],station[i],year[i]] + log.offset[i]
        lam[i] <- exp(log.lam[i])
        daily.count[i] ~ dpois(lam[i])
      }
      
      
      #---------------------------------------------
      # Derived quantities
      #---------------------------------------------
      
      # Total birds across all regions (used to determine national trend)
      for (y in 1:nyear){
        N.total[y] <- sum(N[,y])
      }
      
    }
    ",fill = TRUE)
sink()


parameters.to.save = c("trend",
                       "proc.sd",
                       
                       "rho",
                       
                       "true.total",
                       #"expected.count",
                       #"daily.noise",
                       
                       #"mu",
                       #"N.total",
                       
                       "N"
                       
)


#----------------------------------------------------------------
#Total counts of birds each season at each station
N_sy <- aggregate(ObservationCount ~ station_number + year_adjusted, data = count_df, FUN = sum, na.rm = TRUE)
effort_sy <- aggregate(net.hrs ~ station_number + year_adjusted, data = count_df, FUN = sum, na.rm = TRUE)
N_sy$N.standardized <- N_sy$ObservationCount / effort_sy$net.hrs

ggplot(data = N_sy) + geom_point(aes(y = N.standardized, x = year_adjusted))+facet_grid(station_number~., scales = "free")

rho.relative <- aggregate(N.standardized ~ station_number, data = N_sy, FUN = mean, na.rm = TRUE)
rho.init <- matrix(rep(rho.relative$N.standardized, 3), byrow = TRUE ,nrow =3)/3
#----------------------------------------------------------------


inits <- function() list(trend = rep(0,jags.data$nregion),
                         proc.sd = runif(1,0.1,0.5),
                         rho = rho.init)

out <- jags(data = jags.data,
            model.file = "cmmn_national.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 2,
            n.thin = 5,
            n.iter = 10000,
            n.burnin = 5000)
