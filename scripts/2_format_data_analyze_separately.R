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

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 1: FORMAT MIGRATION COUNT DATA
#******************************************************************************************************************************************
#******************************************************************************************************************************************

#------------------------------------------------------------------------------------------------------------------------------------------
# Read/format data from BSC (Canadian data)
#------------------------------------------------------------------------------------------------------------------------------------------

#Data provided by Danielle Ethier (BSC).  Has been pre-processed/cleaned.  Does not include offsets or effort (i.e., hours nets were open).

dat_can <- rbind(read.csv("../data/migration_counts/CAN/LPBO.BLPW.2018.csv") %>% add_column(., station = "LPBO"),
                 read.csv("../data/migration_counts/CAN/PEPBO.BLPW.2018.csv") %>% add_column(., station = "PEPBO"),
                 read.csv("../data/migration_counts/CAN/TCBO.BLPW.2018.csv") %>% add_column(., station = "TCBO"))

# Create a column to distinguish specific sites within a particular station
dat_can$area <- 1
dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO2")] <- 2
dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO3")] <- 3

# Number of days with data each year at each station
day_range_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = range)
ndays_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = length)
years_per_station <- aggregate(YearCollected ~ SurveyAreaIdentifier + season, data = dat_can, FUN = range)

# Ensure that all stations/days/years/seasons are included as data
area_season_combinations_can <- unique(dat_can[,c("SurveyAreaIdentifier","season")])
dat_combined_can = data.frame()
for (i in 1:nrow(area_season_combinations_can)){
  dat <- subset(dat_can, SurveyAreaIdentifier == area_season_combinations_can$SurveyAreaIdentifier[i] & season == area_season_combinations_can$season[i])
  if (nrow(dat) == 0) next
  min_doy <- min(dat$doy)
  max_doy <- max(dat$doy)
  min_year <- min(dat$YearCollected)
  max_year <- max(dat$YearCollected)
  
  # Create a "full" dataframe to store counts on all days (including NAs)
  dat_full <- expand.grid(YearCollected = seq(min_year,max_year),
                          doy = seq(min_doy,max_doy))
  
  # Fill with counts
  dat_full <- merge(dat_full, dat, all.x = TRUE)
  
  # Ensure relevant data is filled in
  dat_full$SurveyAreaIdentifier = dat$SurveyAreaIdentifier[1]
  dat_full$station = dat$station[1]
  dat_full$area = dat$area[1]
  dat_full$season = dat$season[1]
  dat_full$min_doy = min_doy
  dat_full$max_doy = max_doy
  dat_full$min_year = min_year
  dat_full$max_year = max_year
  
  
  # Plot daily counts in each season
  daily_count_plot <- ggplot(data = dat_full) +
    geom_line(aes(x = doy, y = ObservationCount), col = "blue")+
    facet_wrap(.~YearCollected)+
    ggtitle(dat_full$SurveyAreaIdentifier[1])+
    theme_bw()
  
  #pdf(file = paste0("../figures/station_plots/",dat_full$season[1],"_",dat_full$SurveyAreaIdentifier,".pdf"), width = 20,height=10)
  #print(daily.count.plot)
  #dev.off()
  
  dat_combined_can = rbind(dat_combined_can, dat_full)
}
dat_combined_can$country = "CAN"

#------------------------------------------------------------------------------------------------------------------------------------------
# Read/format data from Ricky Dunn (USA data)
#------------------------------------------------------------------------------------------------------------------------------------------
dat_usa <- rbind(readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - AIMS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - BIBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - BSBO fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - BSBO spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - FBBS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - MCCS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - MCCS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - FBBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - KWRS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"))

# datasets with different column names
tmp = readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - PARC fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")

colnames(tmp) <- colnames(dat_usa)

dat_usa <- rbind(dat_usa, tmp)
rm(tmp)    

# Fill net N/100 net-hr column
dat_usa$`N/100 net-hr` = dat_usa$N/dat_usa$`Net-hrs`*100
dat_usa = na.omit(dat_usa)
summary(dat_usa)

# Change column names
colnames(dat_usa) = c("station","YearCollected","doy","net.hrs","ObservationCount","N.per.100.net.hrs","season")
dat_usa$area = 1

# Load analysis windows for each station and restrict data to those dates
us_windows = read.csv("../data/migration_counts/USA/US_station_windows.csv")
dat_usa2 = data.frame()
for (i in 1:nrow(us_windows)){
  x = us_windows[i,]
  station = us_windows$station[i]
  area = us_windows$area[i]
  season = us_windows$season[i]
  
  # data inside that range
  dat_usa2 = rbind(dat_usa2, subset(dat_usa, station == x$station & area == x$area & season == x$season & doy >= x$start_date & doy <= x$end_date))
}
dat_usa = dat_usa2
dat_usa = merge(dat_usa, us_windows, all.x = TRUE)


area_season_combinations_usa <- unique(dat_usa[,c("station","season")])
dat_combined_usa = data.frame()
for (i in 1:nrow(area_season_combinations_usa)){
  dat <- subset(dat_usa, station == area_season_combinations_usa$station[i] & season == area_season_combinations_usa$season[i])
  if (nrow(dat) == 0) next
  min_doy <- min(dat$doy)
  max_doy <- max(dat$doy)
  min_year <- min(dat$YearCollected)
  max_year <- max(dat$YearCollected)
  
  # Create a "full" dataframe to store counts on all days (including NAs)
  dat_full <- expand.grid(YearCollected = seq(min_year,max_year),
                          doy = seq(min_doy,max_doy))
  
  # Fill with counts
  dat_full <- merge(dat_full, dat, all.x = TRUE)
  
  # Ensure relevant data is filled in
  dat_full$station = dat$station[1]
  dat_full$station = dat$station[1]
  dat_full$area = dat$area[1]
  dat_full$season = dat$season[1]
  dat_full$min_doy = min_doy
  dat_full$max_doy = max_doy
  dat_full$min_year = min_year
  dat_full$max_year = max_year
  
  # Plot daily counts in each season
  daily_count_plot <- ggplot(data = dat_full) +
    geom_line(aes(x = doy, y = ObservationCount), col = "blue")+
    facet_wrap(.~YearCollected)+
    ggtitle(dat_full$station[1])+
    theme_bw()
  
  dat_combined_usa = rbind(dat_combined_usa, dat_full)
}
dat_combined_usa$country = "USA"

#------------------------------------------------------------------------------------------------------------------------------------------
# Combine CAN and USA data into a single dataframe; generate plots
#------------------------------------------------------------------------------------------------------------------------------------------
dat_combined = dplyr::bind_rows(dat_combined_can, dat_combined_usa)

#----------
# For any sites without recorded net hours, fill in with median
dat_combined$net.hrs[which(is.na(dat_combined$net.hrs))] = median(dat_combined$net.hrs, na.rm = TRUE)
#----------

rm(list=setdiff(ls(), c("dat_combined")))

# Limit to data collected after 2005
dat_combined = subset(dat_combined, YearCollected >= 2006)

#----------------------------------
# Plots of daily counts
#----------------------------------

# # CAN Spring
# CAN_Spring_plot <- ggplot(data = subset(dat_combined, country == "CAN" & season == "Spring")) +
#   geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
#   facet_grid(station~YearCollected, scales = "free_y")+
#   xlab("Day of Year")+
#   ylab("Count")+
#   theme_bw()
# pdf(file = paste0(file = "../figures/CAN_Spring_plot.pdf"), width = 30,height=6)
# print(CAN_Spring_plot )
# dev.off()
# 
# # CAN Fall
# CAN_Fall_plot <- ggplot(data = subset(dat_combined, country == "CAN" & season == "Fall")) +
#   geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
#   facet_grid(station~YearCollected, scales = "free_y")+
#   xlab("Day of Year")+
#   ylab("Count")+
#   theme_bw()
# pdf(file = paste0(file = "../figures/CAN_Fall_plot.pdf"), width = 30,height=6)
# print(CAN_Fall_plot )
# dev.off()
# 
# # USA Spring
# USA_Spring_plot <- ggplot(data = subset(dat_combined, country == "USA" & season == "Spring")) +
#   geom_point(aes(x = doy, y = N.per.100.net.hrs), col = "blue")+
#   facet_grid(station~YearCollected, scales = "free_y")+
#   xlab("Day of Year")+
#   ylab("Count(N/100 net hours)")+
#   theme_bw()
# pdf(file = paste0(file = "../figures/USA_Spring_plot.pdf"), width = 30,height=6)
# print(USA_Spring_plot )
# dev.off()
# 
# # USA Fall
# USA_Fall_plot <- ggplot(data = subset(dat_combined, country == "USA" & season == "Fall")) +
#   geom_point(aes(x = doy, y = N.per.100.net.hrs), col = "blue")+
#   facet_grid(station~YearCollected, scales = "free_y")+
#   xlab("Day of Year")+
#   ylab("Count(N/100 net hours)")+
#   theme_bw()
# pdf(file = paste0(file = "../figures/USA_Fall_plot.pdf"), width = 30,height=6)
# print(USA_Fall_plot )
# dev.off()

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 2: SET UP ANALYSIS
#******************************************************************************************************************************************
#******************************************************************************************************************************************

#---------------------------------------------------------------------------------------------
# Bayesian analysis (separate for each station)
#---------------------------------------------------------------------------------------------

# Model with random effects for annual phenology

sink("cmmn_separate_randompeak.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics
      #---------------------------------------------

      # Trend and annual noise
      log.trend ~ dnorm(0,1)
      proc.sd ~ dunif(0,2)
      proc.tau <- pow(proc.sd,-2)
      
      # Separate intercepts for each sub-area within a station (for stations with multiple sub-areas)
      for (A in 1:narea){
        intercept[A] ~ dunif(0,upper_limit)
        log.intercept[A] <- log(intercept[A])
        
        for (y in 1:nyear){
        
          mu[A,y] <- log.intercept[A] + log.trend*(y-1)
          logN[A,y] <- mu[A,y] + noise[A,y]
          noise[A,y] ~ dnorm(0,proc.tau)
          N[A,y] <- exp(logN[A,y])
        
        } # close year loop
      }
      
      for (y in 1:nyear){
        N.total[y] <- sum(N[,y])
        
      }
      
      #---------------------------------------------
      # Model for daily counts
      #---------------------------------------------

      # Parameters describing the mean date of migration, and variation in that peak date among years
      mean.migrate.HYPERMEAN ~ dunif(1,nday)
      mean.migrate.HYPERSD ~ dunif(0,nday)
      mean.migrate.HYPERTAU <- pow(mean.migrate.HYPERSD, -2)
      
      # Parameter describing the width of the migration window (assume this window is constant)
      sd.migrate ~ dunif(0,nday)
         
      # Magnitude of observation error can differ among sub-areas at each station
      for (A in 1:narea){
        daily.noise.sd[A] ~ dunif(0,2)
        daily.noise.tau[A] <- pow(daily.noise.sd[A],-2)
      }
      
        for (y in 1:nyear){
        
        # Migration period is assumed to be the same at each sub-area within a year
        mean.migrate[y] ~ dnorm(mean.migrate.HYPERMEAN, mean.migrate.HYPERTAU)
        
          for (A in 1:narea){
          
            for (d in 1:nday){
              
              # Each day, at each sub-area, within each year, estimate the proportion of total annual detections.
              norm.density[A,d,y] <- 1/(sqrt(2*pi)*sd.migrate)*exp(-((d-mean.migrate[y])^2/(2*sd.migrate^2)))
              
              # Expected count on each day (probability density * total abundance)
              expected.count[A,d,y] <- norm.density[A,d,y] * N[A,y]
              
              # Daily observation error
              daily.noise[A,d,y] ~ dnorm(0,daily.noise.tau[A])
              
              log.lambda[A,d,y] <- log(expected.count[A,d,y]) + daily.noise[A,d,y]
              
            } # close day loop
            
          } # close area loop
      
      } # close year loop 
      
      for (i in 1:nobs){
        log.lam[i] <- log.lambda[area[i],day[i],year[i]] + log.offset[i]
        lam[i] <- exp(log.lam[i])
        daily.count[i] ~ dpois(lam[i])
      }
      
      
      #---------------------------------------------
      # Goodness-of-fit
      #---------------------------------------------
      # for (i in 1:nobs){
      # 
      #   #-----------------------------------------------
      #   # Assess fit at level 1 (deviations from lambda)
      #   #-----------------------------------------------
      #   sim.count.1[i] ~ dpois(lam[i])
      # 
      #   sqerror.obs.1[i] <- pow(daily.count[i] - lam[i], 2)
      #   sqerror.sim.1[i] <- pow(sim.count.1[i] - lam[i], 2)
      # 
      #   X2.obs.1[i]      <- sqerror.obs.1[i]/lam[i]
      #   X2.sim.1[i]      <- sqerror.sim.1[i]/lam[i]
      # 
      #   # #-----------------------------------------------
      #   # # Assess fit at level 2 (deviations from expected)
      #   # #-----------------------------------------------
      #   # sim.noise[i] ~ dnorm(0,daily.noise.tau)
      #   # sim.count.2[i] ~ dpois(exp(log(expected[i]) + sim.noise[i]))
      #   # 
      #   # sqerror.obs.2[i] <- pow(daily.count[i] - expected[i], 2)
      #   # sqerror.sim.2[i] <- pow(sim.count.2[i] - expected[i], 2)
      #   # 
      #   # X2.obs.2[i]      <- sqerror.obs.2[i]/expected[i]
      #   # X2.sim.2[i]      <- sqerror.sim.2[i]/expected[i]
      # 
      # }
      # 
      # chi2.obs.1 <- sum(X2.obs.1[])
      # chi2.sim.1 <- sum(X2.sim.1[])
      # 
      # # chi2.obs.2 <- sum(X2.obs.2[])
      # # chi2.sim.2 <- sum(X2.sim.2[])
      
    }
    ",fill = TRUE)
sink()


write.csv(dat_combined, file = "./processed_data/dat_combined.csv", row.names = FALSE)

# ******************************************************************************************************************************************
# ******************************************************************************************************************************************
# PART 3: RUN ANALYSIS
# ******************************************************************************************************************************************
# ******************************************************************************************************************************************

numCores <- detectCores() # Detect number of cores on machine
numCores <- numCores - 2  # Reserve 2 cores for other tasks
registerDoParallel(numCores) # Number of cores to use for parallel processing

station_season_combinations = unique(dat_combined[,c("station","season")])
station_season_combinations = subset(station_season_combinations, season == "Fall")
allresults = foreach(i = (1:nrow(station_season_combinations)), .combine = list, .packages = c("jagsUI")) %dopar% { #nrow(station_season_combinations)
#for (i in which(station_season_combinations$station == "LPBO" & station_season_combinations$season == "Fall")){
  
  start_time <- Sys.time()
  
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
                   
                   upper_limit = max(aggregate(ObservationCount~year_adjusted + area, data = dat, FUN = sum)$ObservationCount)*10,
                   
                   log.offset = log(dat$net.hrs),
                   pi = pi
  )
  
  inits <- function() list(log.trend = rnorm(1,0,0.05),
                           proc.sd = runif(1,0.1,0.5))
  
  parameters.to.save = c("log.trend",
                         "proc.sd",
                         
                         "log.intercept",
                         "daily.noise.sd",
                         "mean.migrate.HYPERMEAN",
                         "mean.migrate.HYPERSD",
                         
                         "mean.migrate",
                         "sd.migrate",
                         
                         "expected.count",
                         "mu",
                         "N.total",
                         "N"#,
                         
                         # Goodness of fit testing
                         #"chi2.obs.1",
                         #"chi2.sim.1"
         
  )
  
  out <- jags(data = jags.data,
              model.file = "cmmn_separate_randompeak.jags",
              parameters.to.save = parameters.to.save,
              inits = inits,
              n.chains = 2,
              n.thin = 50,
              n.iter = 150000,
              n.burnin = 100000)
    
  # Use this code to automatically select suitable number of iterations
  # modelFit <- autorun.jags(model="cmmn_separate_randompeak.jags", 
  #                          monitor=parameters.to.save, 
  #                          data=jags.data, n.chains=3,
  #                          method="parallel", 
  #                          startburnin = 25000,
  #                          psrf.target=1.02)
  
  max(unlist(out$Rhat),na.rm = TRUE)
  mean(unlist(out$Rhat) > 1.10,na.rm = TRUE)
  
  end_time <- Sys.time()
  run_time <- end_time - start_time
  
  save(out, file = paste0("./jags_output/",station_season_combinations$station[i],"_",station_season_combinations$season[i],"_randompeak.RData"))
  
  out
  
  # #******************************************************************************************************************************************
  # #******************************************************************************************************************************************
  # # PART 4: SUMMARIZE RESULTS
  # #******************************************************************************************************************************************
  # #******************************************************************************************************************************************
  # 
  # # Annual indices
  # N.area = melt(apply(out$sims.list$N,c(2,3),function(x) quantile(x,0.500)), varnames = c("area","year"), value.name = "index.500")
  # N.area$index.025 = melt(apply(out$sims.list$N,c(2,3),function(x) quantile(x,0.025)), varnames = c("area","year"), value.name = "index.025")$index.025
  # N.area$index.975 = melt(apply(out$sims.list$N,c(2,3),function(x) quantile(x,0.975)), varnames = c("area","year"), value.name = "index.975")$index.975
  # 
  # N.overall = data.frame(area = "Overall",
  #                        year = 1:jags.data$nyear,
  #                        index.500 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.500)),
  #                        index.025 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.025)),
  #                        index.975 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.975)))
  # 
  # N.area$area = factor(N.area$area)
  # 
  # N.overall = rbind(N.overall,N.area)
  # N.overall$year = N.overall$year + min(dat$YearCollected) - 1
  # 
  # # area.plot = ggplot( data = N.overall ) +
  # #   geom_errorbar(aes(x = year, ymin = log(index.025), ymax = log(index.975), col = factor(area)), width = 0)+
  # #   geom_point(aes(x = year, y = log(index.500), col = factor(area)))+
  # #   facet_grid(area~., scales = "free")+
  # #   ylab("Index")+
  # #   xlab("Year")+
  # #   scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(N.overall$area)), "Dark2"), name = "Station/Area")+
  # #   theme_bw()
  # #
  # # print(area.plot)
  # #
  # overall.plot = ggplot( data = subset(N.overall, area == "Overall") ) +
  #   geom_errorbar(aes(x = year, ymin = log(index.025), ymax = log(index.975)), width = 0)+
  #   geom_point(aes(x = year, y = log(index.500)))+
  #   ylab("Index")+
  #   xlab("Year")+
  #   ggtitle(paste0(dat$station[1]," (Overall)"))+
  #   theme_bw()
  # 
  # print(overall.plot)
  # 
  # #----------------------------------------------
  # # Overlay observed counts on expected counts
  # #----------------------------------------------
  # 
  # daily.est.500 = melt(apply(out$sims.list$expected.count, c(2,3,4), FUN = function(x) quantile(x, 0.5)),
  #                      value.name = "expected.500", varnames = c("area","doy_adjusted","year"))
  # daily.est.025 = melt(apply(out$sims.list$expected.count, c(2,3,4), FUN = function(x) quantile(x, 0.025)),
  #                      value.name = "expected.025", varnames = c("area","doy_adjusted","year"))
  # daily.est.975 = melt(apply(out$sims.list$expected.count, c(2,3,4), FUN = function(x) quantile(x, 0.975)),
  #                      value.name = "expected.975", varnames = c("area","doy_adjusted","year"))
  # 
  # daily.est = merge(merge(daily.est.500,daily.est.025), daily.est.975)
  # daily.est$doy = daily.est$doy_adjusted + min(dat$doy) - 1
  # daily.est$YearCollected = daily.est$year + min(dat$YearCollected) - 1
  # 
  # daily.plot = ggplot(data = subset(daily.est, YearCollected >= 2015))+
  #   #geom_line(aes(x = doy, y = expected.500), col = "black")+
  #   #geom_ribbon(aes(x = doy, ymin = expected.025, ymax = expected.975), alpha = 0.3)+
  #   
  #   geom_point(data = subset(dat, YearCollected >= 2015), aes(x = doy, y = ObservationCount), col = "blue", alpha = 0.3)+
  #   xlab("Day of Year")+
  #   ylab("Daily Estimated Total")+
  #   ggtitle(dat$station[1])+
  #   facet_grid(.~YearCollected, scales = "free")+
  #   theme_bw()+
  #   theme(axis.text.x = element_text(size = 5))
  # print(daily.plot)
  # 
  # # Compare expected and actual seasonal totals
  # 
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
  # 
}

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 5: PLOT RESULTS
#******************************************************************************************************************************************
#******************************************************************************************************************************************
rm(list=setdiff(ls(), c("dat_combined")))

station_season_combinations = unique(dat_combined[,c("station","season")])

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
                                  index.025 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.025)),
                                  index.975 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.975)),
                                  
                                  index.500.rescaled = apply(out$sims.list$N.total/out$sims.list$N.total[,1],2,function(x) quantile(x,0.500)),
                                  index.025.rescaled = apply(out$sims.list$N.total/out$sims.list$N.total[,1],2,function(x) quantile(x,0.025)),
                                  index.975.rescaled = apply(out$sims.list$N.total/out$sims.list$N.total[,1],2,function(x) quantile(x,0.975)),
                                  
                                  
                                  derived.trend.500 = quantile(derived.trend,0.5),
                                  derived.trend.025 = quantile(derived.trend,0.025),
                                  derived.trend.975 = quantile(derived.trend,0.975),
                                  
                                  mean.trend.500 = quantile(out$sims.list$log.trend,0.5),
                                  mean.trend.025 = quantile(out$sims.list$log.trend,0.025),
                                  mean.trend.975 = quantile(out$sims.list$log.trend,0.975),
                                  
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


#----------------------
# Summarized results
#----------------------

results.fall = ggplot( data = subset(results_summary, season == "Fall" ) ) +
  #geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  #geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.025), ymax = log(index.975)), width = 0, col = "red")+
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


results.spring = ggplot( data = subset(results_summary, season == "Spring" ) ) +
  #geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  #geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.025), ymax = log(index.975)), width = 0, col = "blue")+
  geom_point(aes(x = year, y = log(index.500)), col = "blue")+
  ylab("Population index")+
  xlab("Year")+
  facet_grid(station~season, scales = "free")+
  
  #scale_y_continuous(breaks = log( c(1/5 , 1,  5)),
  #                   labels = c("5-fold decrease", "no change", "5-fold increase"),
  #                  minor_breaks = NULL)+
  
  #coord_cartesian(ylim = log( c(1/16 , 1,  16)))+
  theme_bw()

results.spring


scaled.results = ggplot( data = subset(results_summary, year > 2006) ) +
  geom_hline(yintercept = 0, linetype = 1, col = "gray85", size = 2)+
  geom_hline(yintercept = log(c(1/5,5)), linetype = 1, col = "gray85", size = 1)+
  
  geom_errorbar(aes(x = year, ymin = log(index.025.rescaled), ymax = log(index.975.rescaled), col = season), width = 0)+
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

# Log trends
trend.df = unique(results_summary[,c("station","season",
                                     "mean.trend.025","mean.trend.500","mean.trend.975",
                                     "derived.trend.025","derived.trend.500","derived.trend.975")])

trend.results = ggplot( data = trend.df) +
  
  geom_errorbarh(aes(y = station, xmin = mean.trend.025, xmax = mean.trend.975, col = season), height = 0)+
  geom_point(aes(y = station, x = mean.trend.500, col = season))+
  ylab("Station")+
  xlab("Year")+
  facet_grid(.~season)+
  coord_cartesian(xlim=c(-0.5,0.5))+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = c("blue","red"))+
  theme_bw()

print(trend.results)
trend.results