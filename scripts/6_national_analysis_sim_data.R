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

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Read in seasonal counts at stations
#--------------------------------------------------------------------
#--------------------------------------------------------------------
nregion = 3

# Read in data and subset to Fall (for demonstration)
dat_combined <- read.csv("./processed_data/dat_combined.csv")
dat_combined <- subset(dat_combined, season == "Fall")

stations_to_keep <- c("LPBO","MGBO","MCCS","PEPBO")
dat_combined <- subset(dat_combined, station %in% stations_to_keep & area == 1)


min_date <- min(dat_combined$doy, na.rm = TRUE)
min_year <- min(dat_combined$YearCollected, na.rm = TRUE)

dat_combined$doy_adj <- dat_combined$doy - min_date + 1
dat_combined$year_adj <- dat_combined$YearCollected - min_year + 1
dat_combined$station_number <- as.numeric(factor(dat_combined$station, levels = stations_to_keep))

nday <- max(dat_combined$doy_adj)
nstation <- max(dat_combined$station_number)
nyear <- max(dat_combined$year_adj)
      
#-------------------------
unique(dat_combined[,c("station","station_number")]) %>% arrange(. , station_number)
#-------------------------

# Create a nday x nstation x nyear array to store daily counts and effort offsets
daily.count <- array(NA, dim = c(nday,nstation,nyear))
daily.offsets <- array(NA, dim = c(nday,nstation,nyear))


for (i in 1:nrow(dat_combined)){
  daily.count[dat_combined$doy_adj[i],dat_combined$station_number[i],dat_combined$year_adj[i]] = dat_combined$ObservationCount[i]
  daily.offsets[dat_combined$doy_adj[i],dat_combined$station_number[i],dat_combined$year_adj[i]] = dat_combined$net.hrs[i]
  
  }
daily.offsets[is.na(daily.offsets)] <- median(daily.offsets, na.rm=TRUE)

mean(daily.count / daily.offsets, na.rm = TRUE)

# Plausible values of rho (per season)
rho.daily = daily.count / daily.offsets


#********************
# Isotope sampling
#********************

N.isotope = array(0, dim = c(nregion,nstation,nyear))

#Station 1 gets birds from region 1
N.isotope[1,1,] = 25
N.isotope[1,2,] = 25
N.isotope[2,3,] = 25
N.isotope[3,4,] = 25

# Remove all isotope data except for 5th year
N.isotope[,1:nstation,1:4] = NA
N.isotope[,1:nstation,6:nyear] = NA

N.station.sampled <- array(25, dim = c(nstation,nyear))


sink("cmmn_part2.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics in each region
      #---------------------------------------------
 
      for (r in 1:nregion){
        trend[r] ~ dnorm(0,0.1) # Regional trends
      }
      
      # Temporal variance in trend
      #proc.sd ~ dunif(0,2)
      #proc.tau <- pow(proc.sd,-2)
      
      # True (unobserved) population dynamics in each region
      for (y in 1:nyear){
        for (r in 1:nregion){
        
          # Exponential population model
          logN[r,y] <- logN0[r] + trend[r] * (y-1) #+ noise[r,y]
          #noise[r,y] ~ dnorm(0,proc.tau)
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
           rho.variable[r,s] ~ dgamma(0.001,0.001) 
           rho[r,s] <- rho.variable[r,s] * rho.fix[r,s]
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
      
        mean.migrate[s] ~ dunif(1,nday)
        sd.migrate[s] ~ dunif(0,nday/2)
      
        daily.noise.sd[s] ~ dunif(0,2)
        daily.noise.tau[s] <- pow(daily.noise.sd[s],-2)
      
        for (y in 1:nyear){
          for (d in 1:nday){
            
            norm.density[d,s,y] <- 1/(sqrt(2*pi)*sd.migrate[s])*
                exp(-((d-mean.migrate[s])^2/(2*sd.migrate[s]^2)))
            
            # Expected count on each day
            expected.count[d,s,y] <- norm.density[d,s,y] * true.total[s,y] * exp(log.offset[d,s,y])
            
            # Daily observation error
            daily.noise[d,s,y] ~ dnorm(0,daily.noise.tau[s])
            
            lambda[d,s,y] <- exp(log(expected.count[d,s,y]) + daily.noise[d,s,y])
            daily.count[d,s,y] ~ dpois(lambda[d,s,y])
            
          }
        }
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
                       
                       "mean.migrate",
                       "sd.migrate",
                       "daily.noise.sd",
                       
                       "N"
                       
)

#************************************
# Package data for analysis
#************************************
rho.fix <- matrix(0, nrow = nregion, ncol = nstation)

rho.fix <- (apply(N.isotope,c(1,2),sum, na.rm = TRUE) > 0) * 1 #In cases where a transition was never observed, fix it to zero

jags.data = list(
  
  # Regional regions and initial abundances
  nregion = nregion,
  logN0 = rep(0,nregion),
  
  # Number of stations in dataset
  nstation = nstation,
  
  # Number of years in dataset
  nyear = nyear,
  
  # Number of days in each season
  nday = nday,
  
  # Isotope/catchment data
  N.isotope = N.isotope,
  N.station.sampled = N.station.sampled,
  
  daily.count = daily.count,
  
  log.offset = log(daily.offsets),
  pi = pi,
  
  rho.fix = rho.fix #Can set certain transitions to zero if necessary
)

inits <- function() list(trend = rep(0,nregion),
                         rho.variable = matrix(mean(apply(rho.daily,c(2,3),sum, na.rm=TRUE)),nrow = nregion, ncol = nstation))

out <- jags(data = jags.data,
            model.file = "cmmn_part2.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 2,
            n.thin = 5,
            n.iter = 50000,
            n.burnin = 20000)

print(out)


# par(mfrow=c(3,1))
# # Station total estimates
# mu.matrix.est = apply(out$sims.list$true.total, c(2,3), mean)
# mu.matrix
# mu_df_est = melt(mu.matrix.est, varnames = c("station_number","year_adjusted"), value.name = "mu.est")
# mu_df_true = melt(mu.matrix, varnames = c("station_number","year_adjusted"), value.name = "mu.true")
# mu_df = merge(mu_df_est,mu_df_true)
# 
# #plot(mu.est~mu.true, data = mu_df, main = iteration) # Should be a strong positive relationship
# #abline(a = 0, b = 1)
# 
# # Regional total estimates
# N.est = apply(out$sims.list$N, c(2,3), mean)
# N_df_est = melt(N.est, varnames = c("region","year_adjusted"), value.name = "N.est")
# N_df_true = melt(N.matrix, varnames = c("region","year_adjusted"), value.name = "N.true")
# N_df = merge(N_df_est,N_df_true)
# plot(N.est~N.true, data = N_df, main = iteration) # Should be a strong positive relationship
# abline(a = 0, b = 1)
# 
# #Estimates of rho
# rho.est = apply(out$sims.list$rho, c(2,3), mean)
# rho_df_est = melt(rho.est, varnames = c("from_region","to_station"), value.name = "rho.est")
# rho_df_true = melt(rho, varnames = c("from_region","to_station"), value.name = "rho.true")
# rho_df = merge(rho_df_est,rho_df_true)
# #plot(rho.est~rho.true, data = rho_df, main = iteration) # Should be a strong positive relationship
# #abline(a = 0, b = 1)
# 
# par(mfrow=c(1,1))
# 


#Relevant quantities to check

#1) estimate of intitial abundance (regional and overall)
N.est.regional <- apply(out$sims.list$N, c(2,3), mean)
N.est.total <- apply(apply(out$sims.list$N, c(1,3), sum),2,mean)

#2) estimates of rho
rho.est <- melt(apply(out$sims.list$rho, c(2,3), mean), varnames = c("from_region","to_station"))
rho.truth <- melt(rho,varnames = c("from_region","to_station"))

#3) estimates of trend
trend.est <- apply(out$sims.list$trend, c(2), mean)






iteration.summary <- data.frame() 

#Regional abundances (initial)
iteration.summary = rbind(iteration.summary, data.frame(param = paste0("N.intial.region_",1:nregion),
                                                        param.type = "N.region",
                                                        estimate = N.est.regional[,1],
                                                        truth = N.matrix[,1]))

#Regional abundances (final)
iteration.summary = rbind(iteration.summary,data.frame(param = paste0("N.final.region_",1:nregion),
                                                       param.type = "N.region",
                                                       estimate = N.est.regional[,nyear],
                                                       truth = N.matrix[,nyear]))

#Estimates of rho 
iteration.summary = rbind(iteration.summary,data.frame(param = paste0("from.region.",rho.est$from_region,"_to.station.",rho.est$to_station),
                                                       param.type = "rho",
                                                       estimate = rho.est$value,
                                                       truth = rho.truth$value))

#Estimates of trend 
iteration.summary = rbind(iteration.summary,data.frame(param = paste0("trend.",1:nregion),
                                                       param.type = "trend",
                                                       estimate = trend.est,
                                                       truth = region.trend))

iteration.summary$seed = seed

iteration.summary$maxRhat = max(unlist(out$Rhat),na.rm = TRUE)
iteration.summary$meanRhat = mean(unlist(out$Rhat) > 1.2,na.rm = TRUE)

iteration.summary$bias <- iteration.summary$estimate - iteration.summary$truth


