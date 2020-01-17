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
#station_season_combinations = subset(station_season_combinations, season == "Fall")

results_summary = data.frame()
for (i in 1:nrow(station_season_combinations)){
  
  #i = which(station_season_combinations$station == "RUTH" & station_season_combinations$season == "Fall")
  
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
    
    # Examine the parameters that are not converging
    #Which parameters have not converged?
    #not_conv = unlist(out$Rhat)
    #not_conv[which(not_conv > 1.1)]
    #par(mfrow=c(5,2))
    #traceplot(out)
    #par(mfrow=c(1,1))
    
    #----------------------------------------------
    # Overlay observed counts on expected counts (save in separate files)
    #----------------------------------------------

    daily.est.500 = melt(apply(out$sims.list$expected.count, c(2,3,4), FUN = function(x) quantile(x, 0.5)),
                         value.name = "expected.500", varnames = c("area","doy_adjusted","year"))
    daily.est.025 = melt(apply(out$sims.list$expected.count, c(2,3,4), FUN = function(x) quantile(x, 0.025)),
                         value.name = "expected.025", varnames = c("area","doy_adjusted","year"))
    daily.est.975 = melt(apply(out$sims.list$expected.count, c(2,3,4), FUN = function(x) quantile(x, 0.975)),
                         value.name = "expected.975", varnames = c("area","doy_adjusted","year"))

    daily.est = merge(merge(daily.est.500,daily.est.025), daily.est.975)
    daily.est$doy = daily.est$doy_adjusted + min(dat$doy) - 1
    daily.est$YearCollected = daily.est$year + min(dat$YearCollected) - 1

    daily.plot = ggplot(data = daily.est)+
      geom_line(aes(x = doy, y = expected.500), col = "black")+
      geom_ribbon(aes(x = doy, ymin = expected.025, ymax = expected.975), alpha = 0.3)+

      xlab("Day of Year")+
      ylab("Daily Estimated Total")+
      ggtitle(dat$station[1])+
      facet_grid(area~YearCollected, scales = "free")+
      theme_bw()+
      theme(axis.text.x = element_text(size = 5))
    
   daily.plot = daily.plot +
        geom_point(data = dat, aes(x = doy, y = ObservationCount / net.hrs), col = "black", alpha = 0.3)
    
    pdf(file = paste0("../figures/results_station_plots/",dat$season[1],"_",dat$station[1],"_dailyexpected.pdf"), width = 12, height = length(unique(dat$area))*3)
    print(daily.plot)
    dev.off()
    
    # Compare expected and actual seasonal totals

    # dat$year_doy = paste(dat$YearCollected,dat$doy, sep="_")
    # daily.est$year_doy = paste(daily.est$YearCollected,daily.est$doy,sep="_")
    # daily.est.subset = subset(daily.est, year_doy %in% dat$year_doy)
    # 
    # seasonal.total.obs = aggregate(ObservationCount ~ YearCollected, data = dat, FUN = sum)
    # seasonal.total.exp = aggregate(expected.500 ~ YearCollected, data = daily.est.subset, FUN = sum)
    # seasonal.total = merge(seasonal.total.obs,seasonal.total.exp)
    # 
    # obs_vs_exp = ggplot(data = seasonal.total)+
    #   stat_smooth(aes(x = log(expected.500), y = log(ObservationCount)), method = "lm", alpha = 0.2, fill = "dodgerblue")+
    #   geom_point(aes(x = log(expected.500), y = log(ObservationCount)), col = "blue")+
    #   xlab("log(Expected Seasonal Total)")+
    #   ylab("log(Observed Seasonal Total)")+
    #   ggtitle(dat$station[1])+
    #   theme_bw()
    # print(obs_vs_exp)
    # 
    # summary(lm(log(ObservationCount) ~ log(expected.500), data = seasonal.total))

  }
  
}

head(results_summary)

# Remove results for stations that did not converge
max.Rhat.stations <- aggregate(max.Rhat ~ station + season, data = results_summary, FUN = mean)
results_summary <- subset(results_summary, station %in% subset(max.Rhat.stations, max.Rhat <= 1.2)$station)
max.Rhat.stations <- aggregate(max.Rhat ~ station + season, data = results_summary, FUN = mean)


#----------------------------------------------------------------------------
# Plot annual indices for each season
#----------------------------------------------------------------------------

results.fall = ggplot( data = subset(results_summary, season == "Fall" ) ) +
  #geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  #geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.050), ymax = log(index.950)), width = 0, col = "blue")+
  geom_point(aes(x = year, y = log(index.500)), col = "red")+
  ylab("Population index")+
  xlab("Year")+
  facet_grid(station~season, scales = "free")+
  
  scale_color_manual(values = c("red","blue"))+

  theme_bw()

pdf(file = "../figures/results_annual_indices/fall_indices.pdf", width = 8, height = 20)
print(results.fall)
dev.off()

results.spring = ggplot( data = subset(results_summary, season == "Spring" ) ) +
  #geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  #geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.050), ymax = log(index.950)), width = 0, col = "red")+
  geom_point(aes(x = year, y = log(index.500)), col = "red")+
  ylab("Population index")+
  xlab("Year")+
  facet_grid(station~season, scales = "free")+
  
  scale_color_manual(values = c("red","blue"))+
  
  theme_bw()

pdf(file = "../figures/results_annual_indices/spring_indices.pdf", width = 8, height = 20)
print(results.spring)
dev.off()

# scaled.results = ggplot( data = subset(results_summary, year >= 2006) ) +
#   geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
#   geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
#   
#   geom_errorbar(aes(x = year, ymin = log(index.050.rescaled), ymax = log(index.950.rescaled), col = season), width = 0)+
#   geom_point(aes(x = year, y = log(index.500.rescaled), col = season))+
#   ylab("Population change relative to 2006")+
#   xlab("Year")+
#   facet_grid(station~season)+
#   
#   scale_color_manual(values = c("blue","red"))+
#   scale_y_continuous(breaks = log( c(1/5 , 1,  5)),
#                      labels = c("5-fold decrease", "no change", "5-fold increase"),
#                      minor_breaks = NULL)+
#   
#   coord_cartesian(ylim = log( c(1/16 , 1,  16)))+
#   theme_bw()
# 
# scaled.results

# Log-scale trends
# Note that derived trend is an endpoint analysis between 2006 and final year of data
trend.df = unique(results_summary[,c("station","season",
                                     "mean.trend.050","mean.trend.500","mean.trend.950",
                                     "derived.trend.050","derived.trend.500","derived.trend.950")])

trend.results = ggplot( data = trend.df) +
  
  geom_errorbarh(aes(y = station, xmin = mean.trend.050, xmax = mean.trend.950, col = season), height = 0)+
  geom_point(aes(y = station, x = mean.trend.500, col = season))+
  ylab("Station")+
  xlab("Year")+
  facet_grid(.~season)+
  coord_cartesian(xlim=c(-0.35,0.35))+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = c("blue","red"))+
  theme_bw()



pdf(file = "../figures/results_trends/linear_trends.pdf", width = 8, height = 5)
print(trend.results)
dev.off()

#******************************************************************************************************************************************
# PART 3: PLOT STATION TRENDS ON A MAP
#******************************************************************************************************************************************

# Load bam data
bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")

cmmn_coordinates <- read.csv("../data/locations/CMMN_locations.csv")
colnames(cmmn_coordinates) <- c("results_code","station","name","statprov_code","lat","lon")
cmmn_coordinates$station <- as.character(cmmn_coordinates$station)

# Coordinates of US migration stations
station_coordinates <- rbind(data.frame(station = "MCCS", lat = 41.91, lon = -70.54),
                             data.frame(station = "AIMS", lat = 42.99, lon = -70.61),
                             data.frame(station = "KWRS", lat = 41.48, lon = -71.53),
                             data.frame(station = "BIBS", lat = 41.21, lon = -71.58),
                             data.frame(station = "BSBO", lat = 41.61, lon = -83.19),
                             data.frame(station = "FBBO", lat = 39.20, lon = -76.06),
                             data.frame(station = "PARC", lat = 40.16, lon = -79.27),
                             data.frame(station = "CFMS", lat = 64.86, lon = -147.74)
)

station_coordinates <- rbind(station_coordinates, cmmn_coordinates[,c("station","lat","lon")])


trend.df <- merge(trend.df, station_coordinates, all = TRUE)

# Range of years with data at each site
n.years.per.site <- aggregate(YearCollected ~ station + season, FUN = function(x) diff(range(x)), data = dat_combined)
colnames(n.years.per.site)[3] <- "n.years"
trend.df <- merge(trend.df,n.years.per.site, all = TRUE)


trend.df <- na.omit(trend.df)

trend.df <-  st_as_sf(trend.df, coords = c("lon", "lat"),crs = 4326, agr = "constant")

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp")
bcr1 <- subset(bcr1, COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS")
bcr2 <- st_transform(bcr1, crs = crs(bam1))

# Factor for the direction of trend
trend.df$trend.sign <- sign(trend.df$mean.trend.500)

# If 90% CRI overlap zero, set trend direction to 0
trend.df$trend.sign[which(trend.df$mean.trend.050 < 0 & trend.df$mean.trend.950 > 0)] <- 0
trend.df$trend.sign[which(trend.df$mean.trend.050 > 0 & trend.df$mean.trend.950 < 0)] <- 0
trend.df$trend.sign <- factor(trend.df$trend.sign, levels=c(-1,0,1))

#-------------------------
# Fall trends
#-------------------------
trend.map.fall <- ggplot() +   theme_bw() +
  
  geom_sf(data = bcr2, fill = "gray85", col = "gray92") +

  geom_sf(data = subset(trend.df, season == "Fall"), aes(fill = mean.trend.500, shape = trend.sign), alpha = 0.75, size = 3)+
 
  # geom_sf_text(data = subset(trend.df, season == "Fall"), aes(label = station), size = 2, col = "black") +
  
  scale_fill_gradientn(colors = c("darkred","white","darkblue"), limits = c(-0.18,0.18), name = "Station\ntrend")+
  
  scale_shape_manual(values = c(25,21,24), drop = FALSE,name = "trend direction", guide = FALSE)#+
  
  #scale_size_continuous(range = c(1,5), limits = c(1,max(trend.df$n.years)), guide = FALSE)
  
pdf(file = "../figures/results_trends/trend_map_fall.pdf", width = 8, height = 8)
print(trend.map.fall)
dev.off()

#-------------------------
# Spring trends
#-------------------------
trend.map.spring <- ggplot() +   theme_bw() +
  
  geom_sf(data = bcr2, fill = "gray85", col = "gray92") +
  
  geom_sf(data = subset(trend.df, season == "Spring"), aes(fill = mean.trend.500, shape = trend.sign), alpha = 0.75, size = 4)+
  
  #geom_sf_text(data = subset(trend.df, season == "Spring"), aes(label = station), size = 2, col = "black") +
  
  scale_fill_gradientn(colors = c("darkred","white","darkblue"), limits = c(-0.18,0.18), name = "Station\ntrend")+
  
  scale_shape_manual(values = c(25,21,24), drop = FALSE,name = "trend direction", guide = FALSE)#+

#scale_size_continuous(range = c(1,5), limits = c(1,max(trend.df$n.years)), guide = FALSE)


pdf(file = "../figures/results_trends/trend_map_spring.pdf", width = 8, height = 8)
print(trend.map.spring)
dev.off()