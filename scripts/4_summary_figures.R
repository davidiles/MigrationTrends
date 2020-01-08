setwd("~/Projects/MigrationTrends/scripts")

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

#******************************************************************************************************************************************
# PART 1: LOAD PROCESSED MIGRATION COUNT DATA
#******************************************************************************************************************************************
dat_combined <- read.csv("./processed_data/dat_combined.csv")

#******************************************************************************************************************************************
# PART 2: EXTRACT AND SUMMARIZE MODEL RESULTS AT EACH STATION
#******************************************************************************************************************************************

station_season_combinations = unique(dat_combined[,c("station","season")])
station_season_combinations = subset(station_season_combinations, season == "Fall")

results_summary = data.frame()
for (i in 1:nrow(station_season_combinations)){
  
  dat = subset(dat_combined, station == station_season_combinations$station[i] & season == station_season_combinations$season[i])
  
  dat$doy_adjusted = dat$doy - min(dat$doy) + 1
  dat$year_adjusted = dat$YearCollected - min(dat$YearCollected) + 1
  dat = subset(dat, !is.na(ObservationCount))
  
  jags.data = list(daily.count = dat$ObservationCount,
                   nobs = nrow(dat),
                   
                   area = dat$area,
                   narea = max(dat$area),
                   
                   day = dat$doy_adjusted,
                   nday = max(dat$doy_adjusted),
                   
                   year = dat$year_adjusted,
                   nyear = max(dat$year_adjusted),
                   
                   upper_limit = max(aggregate(ObservationCount~year_adjusted + area, data = dat, FUN = sum)$ObservationCount)*5,
                   
                   log.offset = log(dat$net.hrs),
                   pi = pi)
  
  file = paste0("./jags_output/",station_season_combinations$station[i],"_",station_season_combinations$season[i],"_randompeak.RData")
  
  if (file.exists(file)){
    load(file = file)
    
    derived.trend = (log(out$sims.list$N.total[,ncol(out$sims.list$N.total)]) - log(out$sims.list$N.total[,1])) / (ncol(out$sims.list$N.total)-1)
    
    # Store time series of annual index estimates at this station
    N_annual_station = data.frame(station = dat$station[1],
                                  season = dat$season[1],
                                  
                                  year = (1:jags.data$nyear) + min(dat$YearCollected) - 1,
                                  
                                  index.500 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.500)),
                                  index.050 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.050)),
                                  index.950 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.950)),
                                  
                                  index.500.rescaled = apply(out$sims.list$N.total/out$sims.list$N.total[,1],2,function(x) quantile(x,0.500)),
                                  index.050.rescaled = apply(out$sims.list$N.total/out$sims.list$N.total[,1],2,function(x) quantile(x,0.050)),
                                  index.950.rescaled = apply(out$sims.list$N.total/out$sims.list$N.total[,1],2,function(x) quantile(x,0.950)),
                                  
                                  
                                  derived.trend.500 = quantile(derived.trend,0.5),
                                  derived.trend.050 = quantile(derived.trend,0.050),
                                  derived.trend.950 = quantile(derived.trend,0.950),
                                  
                                  mean.trend.500 = quantile(out$sims.list$log.trend,0.5),
                                  mean.trend.050 = quantile(out$sims.list$log.trend,0.050),
                                  mean.trend.950 = quantile(out$sims.list$log.trend,0.950),
                                  
                                  n.chains = out$mcmc.info$n.chains,
                                  n.samples = out$mcmc.info$n.samples,
                                  max.Rhat = max(unlist(out$Rhat),na.rm=TRUE),
                                  Rhat.1.1 = mean(unlist(out$Rhat)>1.1,na.rm=TRUE))
    
    results_summary = rbind(results_summary, N_annual_station)
  }
  
}

head(results_summary)

aggregate(max.Rhat ~ station + season, data = results_summary, FUN = mean)
aggregate(Rhat.1.1 ~ station + season, data = results_summary, FUN = mean)

results.fall = ggplot( data = subset(results_summary, season == "Fall" ) ) +
  #geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  #geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.050), ymax = log(index.950)), width = 0, col = "red")+
  geom_point(aes(x = year, y = log(index.500)), col = "red")+
  ylab("Population index")+
  xlab("Year")+
  facet_grid(station~season, scales = "free")+
  
  scale_color_manual(values = c("red","blue"))+
  #scale_y_continuous(breaks = log( c(1/5 , 1,  5)),
  #                   labels = c("5-fold decrease", "no change", "5-fold increase"),
  #                  minor_breaks = NULL)+
  
  #coord_cartesian(ylim = log( c(1/16 , 1,  16)))+
  theme_bw()

results.fall


scaled.results = ggplot( data = subset(results_summary, year > 2006) ) +
  geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.050.rescaled), ymax = log(index.950.rescaled), col = season), width = 0)+
  geom_point(aes(x = year, y = log(index.500.rescaled), col = season))+
  ylab("Population change relative to 2006")+
  xlab("Year")+
  facet_grid(station~season)+
  
  scale_color_manual(values = c("blue","red"))+
  scale_y_continuous(breaks = log( c(1/5 , 1,  5)),
                     labels = c("5-fold decrease", "no change", "5-fold increase"),
                     minor_breaks = NULL)+
  
  coord_cartesian(ylim = log( c(1/16 , 1,  16)))+
  theme_bw()

scaled.results

# Log-scale trends
# Note that derived trend is an endpoint analysis between 2006 and final year of data
trend.df = unique(results_summary[,c("station","season",
                                     "mean.trend.050","mean.trend.500","mean.trend.950",
                                     "derived.trend.050","derived.trend.500","derived.trend.950")])

trend.results = ggplot( data = trend.df) +
  
  geom_errorbarh(aes(y = station, xmin = derived.trend.050, xmax = derived.trend.950, col = season), height = 0)+
  geom_point(aes(y = station, x = derived.trend.500, col = season))+
  ylab("Station")+
  xlab("Year")+
  facet_grid(.~season)+
  coord_cartesian(xlim=c(-0.5,0.5))+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = c("blue","red"))+
  theme_bw()

print(trend.results)
trend.results


#******************************************************************************************************************************************
# PART 3: PLOT STATION TRENDS ON A MAP
#******************************************************************************************************************************************

# Load bam data
bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")

# Coordinates of migration stations
station_coordinates <- rbind(data.frame(station = "LPBO", lat = 42.58, lon = -80.39),
                             data.frame(station = "PEPBO", lat = 43.94, lon = -76.86),
                             data.frame(station = "TCBO", lat = 48.30, lon = -88.93),
                             data.frame(station = "MCCS", lat = 41.91, lon = -70.54),
                             data.frame(station = "AIMS", lat = 42.99, lon = -70.61),
                             data.frame(station = "KWRS", lat = 41.48, lon = -71.53),
                             data.frame(station = "BIBS", lat = 41.21, lon = -71.58),
                             data.frame(station = "BSBO", lat = 41.61, lon = -83.19),
                             data.frame(station = "FBBO", lat = 39.20, lon = -76.06),
                             data.frame(station = "PARC", lat = 40.16, lon = -79.27)
)
trend.df <- merge(trend.df, station_coordinates, all = TRUE)
trend.df <-  st_as_sf(trend.df, coords = c("lon", "lat"),crs = 4326, agr = "constant")

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS")
bcr2 <- st_transform(bcr1, crs = crs(bam1))

trend.df$trend.sign <- sign(trend.df$derived.trend.500)
trend.df$trend.sign[which(trend.df$derived.trend.050 < 0 & trend.df$derived.trend.950 > 0)] <- 0
trend.df$trend.sign[which(trend.df$derived.trend.050 > 0 & trend.df$derived.trend.950 < 0)] <- 0

ggplot() +
  geom_sf(data = bcr2, fill = "gray85", col = "gray65")+
  geom_sf(data = trend.df, aes(fill = derived.trend.500, shape = factor(trend.sign)), size = 5)+
  scale_fill_gradientn(colors = c("darkred","white","darkblue"), limits = c(-0.2,0.2), name = "trend")+
  scale_shape_manual(values = c(25,21,24), name = "trend direction", guide = FALSE)+
  coord_sf(xlim = c(-100, -50), ylim = c(35, 70), crs = crs(trend.df)) 
