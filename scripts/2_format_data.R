setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  'tidyverse','reshape2','viridis','jagsUI',
  
  'knitr', 'kableExtra')

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

# For any sites without recorded net hours, fill in with 0.001 (doesnt change results)
dat_combined$net.hrs[which(is.na(dat_combined$net.hrs))] = 0.001

rm(list=setdiff(ls(), c("dat_combined")))

# Limit to data collected after 1995
dat_combined = subset(dat_combined, YearCollected >= 1995 & YearCollected <= 2015)


#----------------------------------
# Plots of daily counts
#----------------------------------

# CAN Spring
CAN_Spring_plot <- ggplot(data = subset(dat_combined, country == "CAN" & season == "Spring")) +
  geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count")+
  theme_bw()
pdf(file = paste0(file = "../figures/CAN_Spring_plot.pdf"), width = 30,height=6)
print(CAN_Spring_plot )
dev.off()

# CAN Fall
CAN_Fall_plot <- ggplot(data = subset(dat_combined, country == "CAN" & season == "Fall")) +
  geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count")+
  theme_bw()
pdf(file = paste0(file = "../figures/CAN_Fall_plot.pdf"), width = 30,height=6)
print(CAN_Fall_plot )
dev.off()

# USA Spring
USA_Spring_plot <- ggplot(data = subset(dat_combined, country == "USA" & season == "Spring")) +
  geom_point(aes(x = doy, y = N.per.100.net.hrs), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count(N/100 net hours)")+
  theme_bw()
pdf(file = paste0(file = "../figures/USA_Spring_plot.pdf"), width = 30,height=6)
print(USA_Spring_plot )
dev.off()

# USA Fall
USA_Fall_plot <- ggplot(data = subset(dat_combined, country == "USA" & season == "Fall")) +
  geom_point(aes(x = doy, y = N.per.100.net.hrs), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count(N/100 net hours)")+
  theme_bw()
pdf(file = paste0(file = "../figures/USA_Fall_plot.pdf"), width = 30,height=6)
print(USA_Fall_plot )
dev.off()

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 2: SET UP ANALYSIS
#******************************************************************************************************************************************
#******************************************************************************************************************************************

#---------------------------------------------------------------------------------------------
# Bayesian analysis (separate for each station)
#---------------------------------------------------------------------------------------------

sink("cmmn_separate.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics
      #---------------------------------------------

      # Shared trend and annual noise
      log.trend ~ dnorm(0,1)
      proc.sd ~ dunif(0,2)
      proc.tau <- pow(proc.sd,-2)
      
      # Separate intercepts for each sub-area within a station
      for (A in 1:narea){
        intercept[A] ~ dunif(0,upper_limit*5)
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

      for (A in 1:narea){
        daily.noise.sd[A] ~ dunif(0,2)
        daily.noise.tau[A] <- pow(daily.noise.sd[A],-2)
        
        mean.migrate[A] ~ dunif(1,nday)
        sd.migrate[A] ~ dunif(0,nday)
          
        for (y in 1:nyear){
            for (d in 1:nday){
              
              norm.density[A,d,y] <- 1/(sqrt(2*pi)*sd.migrate[A])*exp(-((d-mean.migrate[A])^2/(2*sd.migrate[A]^2)))
              
              # Expected count on each day
              expected.count[A,d,y] <- norm.density[A,d,y] * N[A,y]
              
              # Daily observation error
              daily.noise[A,d,y] ~ dnorm(0,daily.noise.tau[A])
              
              log.lambda[A,d,y] <- log(expected.count[A,d,y]) + daily.noise[A,d,y]
              
            } # close day loop
            
          } # close year loop
      
      } # close area loop
      
      for (i in 1:nobs){
        log.lam[i] <- log.lambda[area[i],day[i],year[i]] + log.offset[i]
        lam[i] <- exp(log.lam[i])
        daily.count[i] ~ dpois(lam[i])
      }
      
      
      # #---------------------------------------------
      # # Goodness-of-fit
      # #---------------------------------------------
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
      #   #-----------------------------------------------
      #   # Assess fit at level 2 (deviations from expected)
      #   #-----------------------------------------------
      #   sim.noise[i] ~ dnorm(0,daily.noise.tau)
      #   sim.count.2[i] ~ dpois(exp(log(expected[i]) + sim.noise[i]))
      #   
      #   sqerror.obs.2[i] <- pow(daily.count[i] - expected[i], 2)
      #   sqerror.sim.2[i] <- pow(sim.count.2[i] - expected[i], 2)
      #   
      #   X2.obs.2[i]      <- sqerror.obs.2[i]/expected[i]
      #   X2.sim.2[i]      <- sqerror.sim.2[i]/expected[i]
      #   
      # }
      # 
      # chi2.obs.1 <- sum(X2.obs.1[])
      # chi2.sim.1 <- sum(X2.sim.1[])
      # 
      # chi2.obs.2 <- sum(X2.obs.2[])
      # chi2.sim.2 <- sum(X2.sim.2[])
      
    }
    ",fill = TRUE)
sink()

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 3: RUN ANALYSIS
#******************************************************************************************************************************************
#******************************************************************************************************************************************

station_season_combinations = unique(dat_combined[,c("station","season")])

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
                   pi = pi
  )
  
  inits <- function() list(intercept = runif(jags.data$narea,0,100))
  out <- jags(data = jags.data,
              model.file = "cmmn_separate.jags",
              parameters.to.save = c(
                "log.trend",
                "proc.sd",
                
                "log.intercept",
                "daily.noise.sd",
                "mean.migrate",
                "sd.migrate",
                
                "N.total",
                "N"
                
                
              ),
              inits = inits,
              n.chains = 2,
              n.thin = 10,
              n.iter = 2000,
              n.burnin = 1000)
  
  max(unlist(out$Rhat),na.rm = TRUE)
  mean(unlist(out$Rhat) > 1.10,na.rm = TRUE)
  
  save(out, file = paste0("./jags_output/",station_season_combinations$station[i],"_",station_season_combinations$season[i],".RData"))
  #******************************************************************************************************************************************
  #******************************************************************************************************************************************
  # PART 4: SUMMARIZE RESULTS
  #******************************************************************************************************************************************
  #******************************************************************************************************************************************
  
  # Annual indices
  N.area = melt(apply(out$sims.list$N,c(2,3),function(x) quantile(x,0.500)), varnames = c("area","year"), value.name = "index.500")
  N.area$index.025 = melt(apply(out$sims.list$N,c(2,3),function(x) quantile(x,0.025)), varnames = c("area","year"), value.name = "index.025")$index.025
  N.area$index.975 = melt(apply(out$sims.list$N,c(2,3),function(x) quantile(x,0.975)), varnames = c("area","year"), value.name = "index.975")$index.975
  
  N.overall = data.frame(area = "Overall",
                         year = 1:jags.data$nyear,
                         index.500 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.500)),
                         index.025 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.025)),
                         index.975 = apply(out$sims.list$N.total,2,function(x) quantile(x,0.975)))
  
  N.area$area = factor(N.area$area)
  
  N.overall = rbind(N.overall,N.area)
  N.overall$year = N.overall$year + dat$min_year[1] - 1
  
  # area.plot = ggplot( data = N.overall ) +
  #   geom_errorbar(aes(x = year, ymin = log(index.025), ymax = log(index.975), col = factor(area)), width = 0)+
  #   geom_point(aes(x = year, y = log(index.500), col = factor(area)))+
  #   facet_grid(area~., scales = "free")+
  #   ylab("Index")+
  #   xlab("Year")+
  #   scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(N.overall$area)), "Dark2"), name = "Station/Area")+
  #   theme_bw()
  # 
  # print(area.plot)
  # 
  overall.plot = ggplot( data = subset(N.overall, area == "Overall" & year <= 2014) ) +
    geom_errorbar(aes(x = year, ymin = log(index.025), ymax = log(index.975)), width = 0)+
    geom_point(aes(x = year, y = log(index.500)))+
    ylab("Index")+
    xlab("Year")+
    ggtitle(paste0(dat$station[1]," (Overall)"))+
    theme_bw()
  
  print(overall.plot)
  
  # N.sums = data.frame(year = unique(N.overall$year))
  # for (y in unique(dat$year_adjusted)){
  #   seasonal.sum.est = apply(out$sims.list$lam[,which(dat$year_adjusted == y)],1,sum)
  #   N.sums$sum.500[y] = quantile(seasonal.sum.est, 0.5)
  #   N.sums$sum.025[y] = quantile(seasonal.sum.est, 0.025)
  #   N.sums$sum.975[y] = quantile(seasonal.sum.est, 0.975)
  #   
  # }
  # 
  # sum.plot = ggplot( data = subset(N.sums, year <= 2014) ) +
  #   geom_errorbar(aes(x = year, ymin = log(sum.025), ymax = log(sum.975)), width = 0)+
  #   geom_point(aes(x = year, y = log(sum.500)))+
  #   ylab("Index")+
  #   xlab("Year")+
  #   ggtitle(paste0(dat$station[1]," (Overall)"))+
  #   theme_bw()
  # 
  # print(sum.plot)
  # 
  # daily.plot = ggplot(data = dat)+
  #   geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
  #   xlab("Day of Year")+
  #   ylab("Daily Estimated Total")+
  #   facet_grid(area~YearCollected)+
  #   theme_bw()
  # print(daily.plot)
  # 
  # annual.plot = ggplot(data = aggregate(ObservationCount~YearCollected, data = dat, FUN = sum))+
  #   geom_point(aes(x = YearCollected, y = ObservationCount), col = "blue")+
  #   xlab("Year")+
  #   ylab("Seasonal Total")+
  #   theme_bw()
  # print(annual.plot)
}



