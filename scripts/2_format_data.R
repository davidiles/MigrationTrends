setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  'tidyverse','reshape2','viridis','jagsUI',
  
  'knitr', 'kableExtra')

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

#---------------------------------------------------------------------
# Read/format data from BSC
#---------------------------------------------------------------------

#Data provided by Danielle Ethier (BSC).  Has been pre-processed/cleaned.  Does not include offsets or effort (i.e., hours nets were open).

dat_can <- rbind(read.csv("../data/migration_counts/CAN/LPBO.BLPW.2018.csv") %>% add_column(., station = "LPBO"),
                 read.csv("../data/migration_counts/CAN/PEPBO.BLPW.2018.csv") %>% add_column(., station = "PEPBO"),
                 read.csv("../data/migration_counts/CAN/TCBO.BLPW.2018.csv") %>% add_column(., station = "TCBO"))

# Create a column to distinguish specific sites within a particular station
dat_can$area <- 1
dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO2")] <- 2
dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO3")] <- 3

#---------------------------------------------------------------------
# Limit to data collected after 1995

dat_can = subset(dat_can, YearCollected >= 2008)

#---------------------------------------------------------------------

# Number of days with data each year at each station
day_range_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = range)
ndays_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = length)

area_season_combinations <- unique(dat_can[,c("SurveyAreaIdentifier","season")])
dat_combined = data.frame()
for (i in 1:nrow(area_season_combinations)){
  dat <- subset(dat_can, SurveyAreaIdentifier == area_season_combinations$SurveyAreaIdentifier[i] & season == area_season_combinations$season[i])
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
  
  dat_combined = rbind(dat_combined, dat_full)
}

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
              
              lambda[A,d,y] <- exp(log(expected.count[A,d,y]) + daily.noise[A,d,y])
              
            } # close day loop
            
          } # close year loop
      
      } # close area loop
      
      for (i in 1:nobs){
        expected[i] <- expected.count[area[i],day[i],year[i]]
        lam[i] <- lambda[area[i],day[i],year[i]]
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

#---------------------------------------------------------------------------------------------
# Prepare data for analysis
#---------------------------------------------------------------------------------------------

dat = subset(dat_combined, station == "LPBO" & season == "Fall")

dat$doy_adjusted = dat$doy - dat$min_doy + 1
dat$year_adjusted = dat$YearCollected - dat$min_year + 1

jags.data = list(daily.count = dat$ObservationCount,
                 nobs = nrow(dat),
                 
                 area = dat$area,
                 narea = max(dat$area),
                 
                 day = dat$doy_adjusted,
                 nday = max(dat$doy_adjusted),
                 
                 year = dat$year_adjusted,
                 nyear = max(dat$year_adjusted),
                 
                 upper_limit = max(dat$ObservationCount,na.rm=TRUE)*max(dat$doy_adjusted),
                 pi = pi
)

inits <- function() list(intercept = runif(jags.data$narea,0,1000))
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
                                   "N",
                                   "expected",
                                   "lam"
                                   
                                   
            ),
            inits = inits,
            n.chains = 2,
            n.thin = 5,
            n.iter = 2000,
            n.burnin = 1000)

max(unlist(out$Rhat),na.rm = TRUE)
mean(unlist(out$Rhat) > 1.10,na.rm = TRUE)

