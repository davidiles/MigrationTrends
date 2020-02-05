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

# Read in data and subset to Fall (for demonstration)
dat_combined <- read.csv("./processed_data/dat_combined.csv")
dat_combined <- subset(dat_combined, season == "Fall")

stations_to_keep <- c("CFMS","MGBO","MCCS")
dat_combined <- subset(dat_combined, station %in% stations_to_keep)

min_date <- min(dat_combined$doy, na.rm = TRUE)
min_year <- min(dat_combined$YearCollected, na.rm = TRUE)

dat_combined$doy_adj <- dat_combined$doy - min_date + 1
dat_combined$year_adj <- dat_combined$YearCollected - min_year + 1

dat_combined$station_number <- as.numeric(factor(dat_combined$station))

# Create a nday x nstation x nyear array to store daily counts



nregion = 3
nyear = 10
nday = 30

mean.migrate <- 15 #peak migration occurs on day 15
sd.migrate <- 5    #migration window has a normal shape, with sd of 6 days (95% of migration occurs in 20 day window)
daily.noise.sd <- 0.1

#---------------------------------------
# Simulate dynamics
#---------------------------------------

region.abund = rep(1,nregion)
region.trend = seq(-0.1,0.1, length.out = nregion)

# Assume there are at least 2 stations primarily monitoring each region
rho = matrix(0,nrow=nregion, ncol = nstation) #FROM region TO station

rho[1,1] = runif(1,0.005,0.015)
rho[1,2] = runif(1,0.005,0.015)
rho[2,3] = runif(1,0.005,0.015)
rho[2,4] = runif(1,0.005,0.015)
rho[3,5] = runif(1,0.005,0.015)
rho[3,6] = runif(1,0.005,0.015)

#rho[rho == 0] = runif(length(rho[rho == 0]),0,0.001)

# Regional dynamics
proc.sd = 0.2

N.matrix = matrix(NA, nrow = nregion, ncol = nyear) # Relative abundance in each year in each region
mu = array(NA, dim = c(nregion, nstation, nyear)) # Counts arriving at each station in each year from each region

for (r in 1:nregion){
  
  for (y in 1:nyear){
    N.matrix[r,y] = exp(log(region.abund[r]) + region.trend[r]*(y-1) + rnorm(1,0,proc.sd))
    
    for (s in 1:nstation){
      mu[r,s,y] = N.matrix[r,y]*rho[r,s]
      
    }
    
  }
}

mu.matrix = apply(mu,c(2,3), FUN = sum)

#********************
# Isotope sampling
#********************

N.isotope = mu * NA

for (y in 1:nyear){
  for (s in 1:nstation){
    N.isotope[,s,y] = rmultinom(1,25, mu[,s,y])
  }
}
N.station.sampled = apply(N.isotope,c(2,3), sum) #How many birds were sampled at each station in each year

# Remove all isotope data except for first year
N.isotope[,1:nstation,1:4] = NA
N.isotope[,1:nstation,6:nyear] = NA

#********************
# Seasonal totals at each station (effort offset and poisson error)
#********************

log.offset <- log(10000)

#********************
#********************
# Distribute seasonal total across season
#********************
#********************

#daily proportion of total seasonal count that passes each station
migration.density <- dnorm(1:nday, mean = mean.migrate, sd = sd.migrate)

daily.expected.count <- array(NA, dim=c(nday,nstation,nyear))
daily.count <- array(NA, dim=c(nday,nstation,nyear))

for (d in 1:nday){
  for (s in 1:nstation){
    for (y in 1:nyear){
      
      daily.expected.count[d,s,y] <- mu.matrix[s,y] * exp(log.offset) * migration.density[d]
      daily.noise <- rnorm(1,0,daily.noise.sd)
      daily.count[d,s,y] <- rpois(1,exp( log(daily.expected.count[d,s,y]) + daily.noise))
    }
  }
}


sink("cmmn_part2.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics in each region
      #---------------------------------------------
 
      for (r in 1:nregion){
        trend[r] ~ dnorm(0,0.01) # Regional trends
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
           rho.variable[r,s] ~ dunif(0,1) 
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
        sd.migrate[s] ~ dunif(0,20)
      
        daily.noise.sd[s] ~ dunif(0,2)
        daily.noise.tau[s] <- pow(daily.noise.sd[s],-2)
      
        for (y in 1:nyear){
          for (d in 1:nday){
            
            norm.density[d,s,y] <- 1/(sqrt(2*pi)*sd.migrate[s])*
                exp(-((d-mean.migrate[s])^2/(2*sd.migrate[s]^2)))
            
            # Expected count on each day
            expected.count[d,s,y] <- norm.density[d,s,y] * true.total[s,y] * exp(log.offset)
            
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
rho.fix <- (rho > 0)*1

rho.fix <- (apply(N.isotope,c(1,2),sum, na.rm = TRUE) > 0) * 1 #In cases where a transition was never observed, fix it to zero

jags.data = list(
  
  # Regional regions and initial abundances
  nregion = nregion,
  logN0 = log(region.abund),
  
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
  
  log.offset = log.offset,
  pi = pi,
  
  rho.fix = rho.fix #Can set certain transitions to zero if necessary
)

inits <- function() list(trend = region.trend,
                         rho.variable = rho,
                         proc.sd = proc.sd,
                         mean.migrate = rep(mean.migrate,nstation),
                         sd.migrate = rep(sd.migrate,nstation),
                         daily.noise.sd = rep(daily.noise.sd,nstation))

out <- jags(data = jags.data,
            model.file = "cmmn_part2.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 2,
            n.thin = 5,
            n.iter = 25000,
            n.burnin = 10000)

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


