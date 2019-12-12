setwd("~/Projects/MigrationTrends/scripts/model_assumptions/")

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



#----------------------------------------------------------------
# Simulation parameters
#----------------------------------------------------------------
nregion <- 2  # Number of distinct regions
nstation <- 3 # Number of monitoring stations

nyear <- 20    # Length of simulation
proc.sd <- 0.1 # Interannual sd of log-scale abundance

isotope.years <- 5 # Year in which isotope data was collected
N.sampled <- 15    # Number of birds sampled at each station

nday <- 30 # Number of days that each station is open

peak.migration.mean <- 10 # Mean date of peak migration at each station
peak.migration.sd   <- 3  # SD of peak migration date (among-years)
sd.migration.mean <- 10    # Width of migration period

daily.obs.error <- 0.3    # Magnitude of overdispersion in counts

rho.sd <- 0.5  # Magnitude of inter-annual variation in migration strength

col.pal <- RColorBrewer::brewer.pal(nregion,"Paired")[1:nregion] # Define nice colour palette

sink("cmmn_national_sim.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics in each region
      #---------------------------------------------
      
      # Trends and initial abundance in each region [r]
      mean.trend ~ dnorm(0,1) # Grand mean of temporal trend
      sd.trend.spatial ~ dunif(0,1) # Spatial variance in trends
      tau.trend.spatial <- pow(sd.trend.spatial, -2)
      
      for (r in 1:nregion){
      
        logN0.uncert[r] <- logN0[r]
        trend[r] ~ dnorm(0,1)

      }
      
      # Temporal variance in trend
      proc.sd ~ dunif(0,2)
      proc.tau <- pow(proc.sd,-2)
      
      # True (unobserved) population dynamics in each region
      for (y in 1:nyear){
        for (r in 1:nregion){
        
          # Exponential population model
          logN[r,y] <- logN0.uncert[r] + trend[r] * (y-1) + noise[r,y]
          noise[r,y] ~ dnorm(0,proc.tau)
          N[r,y] <- exp(logN[r,y])
        }
      }
      
      #---------------------------------------------
      # Observed numbers of birds at station [s] originating from region [r] in year [y]
      #---------------------------------------------
      

      for (r in 1:nregion){
        for (s in 1:nstation){
        
          rho[r,s] ~ dunif(0,1) # Proportion of birds from region [r] that are captured by station [s] (constant through time)
          
        }
      }
      
      for (s in 1:nstation){
        for (y in 1:nyear){
          for (r in 1:nregion){
          
            mu[r,s,y] <- N[r,y] * rho[r,s] # Seasonal total arriving at station [s] from region [r]  
          
          }
      
          true.total[s,y] <- sum(mu[1:nregion,s,y]) # Total seasonal count at station [s]
          N.isotope[1:nregion,s,y] ~ dmulti(mu[,s,y], N.station.sampled[s,y]) # Multinomial to describe regional composition in this year, at this station
          
        }
      }
      
      #---------------------------------------------
      # Estimate totals at each station each year
      #---------------------------------------------
      
        
      for (s in 1:nstation){
      
        daily.noise.sd[s] ~ dunif(0,2)
        daily.noise.tau[s] <- pow(daily.noise.sd[s],-2)
        
        mean.migrate.HYPERMEAN[s] ~ dunif(1,nday)
        mean.migrate.HYPERSD[s] ~ dunif(0,nday)
        mean.migrate.HYPERTAU[s] <- pow(mean.migrate.HYPERSD[s], -2)
      
        sd.migrate[s] ~ dunif(0,nday)
          
        for (y in 1:nyear){
        
          mean.migrate[s,y] ~ dnorm(mean.migrate.HYPERMEAN[s], mean.migrate.HYPERTAU[s])

          for (d in 1:nday){
            
            norm.density[d,s,y] <- 1/(sqrt(2*pi)*sd.migrate[s])*
                exp(-((d-mean.migrate[s,y])^2/(2*sd.migrate[s]^2)))
            
            # Expected count on each day
            expected.count[d,s,y] <- norm.density[d,s,y] * true.total[s,y]
            
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

#------------------------------------------------------------
# Iterate many times, evaluate results
#------------------------------------------------------------

numCores <- detectCores() # Detect number of cores on machine
numCores <- 4  
registerDoParallel(numCores) # Number of cores to use for parallel processing

foreach(iteration = (1:10000), .combine = rbind, .packages = c("jagsUI","tidyverse","reshape2")) %dopar% {
  
  results.df = data.frame()
  
  #---------------------
  seed = round(runif(1,1,10000))

  file = "./variable_rho.csv"
  if (file.exists(file)){
    results.df = read.csv(file)
    seed = (1:10000)[which(!(1:10000 %in% results.df$seed))]
    seed = sample(seed,1)
  }
  set.seed(seed)
  
  #---------------------
  
  #---------------------------
  # Simulate true population dynamics
  #---------------------------
  # Define trends and relative abundance in each region
  trend <- data.frame(region = 1:nregion,trend = seq(-0.1,0.1,length.out = nregion)) # Trend in each region
  logN0 <- data.frame(region = 1:nregion, logN0 = log(rep(100000,times = nregion))) # Initial abundance
  
  N.df <- expand.grid(year = 1:nyear,region = 1:nregion) %>%
    full_join(x = ., y = trend) %>%
    full_join(x = ., y = logN0)
  
  for (i in 1:nrow(N.df)) N.df$logN[i] <- N.df$logN0[i] + (N.df$y[i]-1)*N.df$trend[i] + rnorm(1,0,proc.sd)
  
  N.df$N <- round(exp(N.df$logN)) # convert to integers for later use in binomial
  N.total <- aggregate(N~year, data = N.df, FUN = sum)
  #---------------------------
  # Simulate migration process
  #---------------------------
  
  # nregion x nstation matrix.  each cell describes the proportion of birds from region [r] that arrive at station [s]
  min.rho = 0.0005
  max.rho = 0.007
  
  migrant.mat <- matrix(0, nrow = nregion, ncol = nstation)
  while (sum(colSums(migrant.mat)==0)>0){ #while loop ensures each station recieves migrants from *somewhere*
    for (r in 1:nregion){
      
      rho <- rep(0,nstation) # Empty vector to store probabilities
      #stations.used <- sample(1:nstation,round(runif(1,1,3))) # Which stations have non-zero probabilities?
      stations.used <- c(r,sample(1:nstation,sample(1:3,1,replace=TRUE))) # Which stations have non-zero probabilities?
      stations.used <- stations.used[which(stations.used > 0 & stations.used <= nstation)]
      
      rho[stations.used] <- runif(length(stations.used),min.rho,max.rho) # Insert random probabilities into those stations
      migrant.mat[r,] <- rho
      
    }
  }
  
  
  # Reorganize matrix
  migrant.mat <- t(migrant.mat)
  colnames(migrant.mat) <- seq(1:nregion)
  rownames(migrant.mat) <- seq(1:nstation)
  
  # Convert to dataframe
  migrant.df <- melt(migrant.mat)
  colnames(migrant.df) <- c("station","region","rho")
  
  N.region.station <- expand.grid(year = 1:nyear, station = 1:nstation, region = 1:nregion) %>%
    full_join(x = . , y = N.df) %>%
    full_join(x = . , y = migrant.df)
  
  # Number of birds captured at station is a binomial process
  N.region.station$rho.variable <- plogis(qlogis(N.region.station$rho) + rnorm(nrow(N.region.station),0,rho.sd))
  N.region.station$N.captured <- rbinom(nrow(N.region.station), N.region.station$N,N.region.station$rho.variable)
  N.region.station$logN <- log(N.region.station$N)
  
  # Total number of birds caught each year at each station
  N.station <- aggregate(N.captured ~ year + station, data = N.region.station, FUN = sum)
  colnames(N.station)[3] <- "N.captured.total"
  N.region.station <- merge(N.region.station, N.station)
  
  # Proportion caught at each station arriving from each region
  N.region.station$prop <- N.region.station$N.captured/N.region.station$N.captured.total
  
  # Organize dataframe
  N.region.station <- N.region.station[order(N.region.station$station,N.region.station$y,N.region.station$region),]
  
  #---------------------------
  # Simulate daily migration counts at each station
  #---------------------------
  migration.phenology.df <- expand.grid(station = 1:nstation,year = 1:nyear)
  
  # Fill in date of peak migration and sd of migration in each year at each station
  migration.phenology.df$peak.migration = rnorm(nrow(migration.phenology.df),peak.migration.mean, peak.migration.sd)
  migration.phenology.df$sd.migration = rep(sd.migration.mean, times = nrow(migration.phenology.df))
  
  N.daily.df <- expand.grid(year = 1:nyear, station = 1:nstation, day = 1:nday) %>%
    full_join(x = . , y = N.station) %>%
    full_join(x = . , y = migration.phenology.df)
  
  N.daily.df$dnorm <- dnorm(N.daily.df$day, mean = N.daily.df$peak.migration, sd = N.daily.df$sd.migration)
  
  N.daily.df$lambda.mu <- N.daily.df$N.captured.total * N.daily.df$dnorm
  N.daily.df$lambda <- exp(log(N.daily.df$N.captured.total * N.daily.df$dnorm) + rnorm(nrow(N.daily.df),0,daily.obs.error))
  N.daily.df$daily.count <- rpois(nrow(N.daily.df), N.daily.df$lambda)
  
  #---------------------------
  # Simulate isotope data
  #---------------------------
  N.isotope.array <- array(NA, dim = c(nregion,nstation,nyear))
  for (s in 1:nstation){
    for (y.iso in isotope.years){
      for (r in 1:nregion){
        
        N.isotope.array[,s,y.iso] = subset(N.region.station, year == y.iso & station == s)$N.captured %>%
          rmultinom(1,size = N.sampled, prob = .)
        
      }
    }
  }
  
  isotope.df <- melt(N.isotope.array, varnames = c("region","station","year"), value.name = "count")
  
  #---------------------------
  # Package data for JAGS analysis
  #---------------------------
  
  # Package data for analysis in JAGS
  daily.count.array <- array(NA, dim = c(nday,nstation,nyear))
  for (s in 1:nstation){
    for (y in 1:nyear){
      for (d in 1:nday){
        
        daily.count.array[d,s,y] = subset(N.daily.df, year == y & station == s & day == d)$daily.count
        
      }
    }
  }
  
  jags.data <- list(nyear = nyear,
                    
                    nregion = nregion,
                    
                    nstation = nstation,
                    
                    # Estimate of initial abundance (from BAM/BBS/e-bird)
                    logN0 = subset(N.df, year == 1)$logN0,
                    
                    # Isotope samples
                    N.station.sampled = array(N.sampled, dim = c(nstation,nyear)),
                    N.isotope = N.isotope.array,
                    
                    # Daily totals each day of the season for each year
                    nday = nday,
                    daily.count = daily.count.array,
                    pi = pi
  )
  
  dimnames(jags.data$N.isotope)[[1]] <- paste0("region",1:nregion)
  dimnames(jags.data$N.isotope)[[2]] <- paste0("station",1:nstation)
  dimnames(jags.data$N.isotope)[[3]] <- paste0("year",1:nyear)
  
  dimnames(jags.data$daily.count)[[1]] <- paste0("day",1:nday)
  dimnames(jags.data$daily.count)[[2]] <- paste0("station",1:nstation)
  dimnames(jags.data$daily.count)[[3]] <- paste0("year",1:nyear)
  
  #---------------------------
  # Analyze with JAGS
  #---------------------------
  
  inits <- function() list(trend = trend$trend,
                           rho = t(migrant.mat))
  
  out <- jags(data = jags.data,
              model.file = "cmmn_national_sim.jags",
              parameters.to.save = c("logN0.uncert",
                                     "trend",
                                     "proc.sd",
                                     
                                     #"daily.noise.sd",
                                     
                                     #"mean.migrate.HYPERMEAN",
                                     #"mean.migrate.HYPERSD",
                                     #"sd.migrate",
                                     #"mean.migrate",
                                     
                                     "rho",
                                     
                                     "N.total",
                                     "N"
              ),
              inits = inits,
              n.chains = 2,
              n.thin = 5,
              n.iter = 50000,
              n.burnin = 30000)
  
  max(unlist(out$Rhat),na.rm=TRUE)
  
  #------------------------------------------------
  # Evaluate bias/precision in trend estimates
  #------------------------------------------------
  
  ########################
  # Overall national trend
  ########################
  # Actual overall 20-year rate of change
  trend.overall.true = log(N.total$N[nyear]/N.total$N[1])/ (nyear-1)
  
  # Estimated overall 20-year rate of change
  trend.overall.est = log(out$sims.list$N.total[,nyear]/out$sims.list$N.total[,1])/ (nyear-1)
  
  # Mean bias in overall trend (derived)
  bias.trend.overall = mean(trend.overall.est - trend.overall.true)
  
  # Precision (width of 95% credible interval)
  precision.trend.overall = diff(quantile(trend.overall.est,c(0.025,0.975)))
  
  ########################
  # Regional trends - 2 ways of calculating: using 1) the linear slope parameter or 2) the derived trend
  ########################
  bias.trend.region.derived = bias.trend.region.linear = c()
  precision.trend.region.derived = precision.trend.region.linear = c()
  
  for (i in 1:nregion){
    trend.region.est.derived = log(out$sims.list$N[,i,nyear]/out$sims.list$N[,i,1])/ (nyear-1)
    trend.region.true.derived = (log(subset(N.df,year == nyear & region == i)$N / subset(N.df,year == 1 & region == i)$N)/(nyear-1))
    
    trend.region.est.linear = out$sims.list$trend[,i]
    trend.region.true.linear = trend$trend[i]
    
    bias.trend.region.derived[i] = mean(trend.region.est.derived - trend.region.true.derived)
    bias.trend.region.linear[i] = mean(trend.region.est.linear - trend.region.true.linear)
    
    precision.trend.region.derived[i] = diff(quantile(trend.region.est.derived,c(0.025,0.975)))
    precision.trend.region.linear[i] = diff(quantile(trend.region.est.linear,c(0.025,0.975)))
    
  }
  
  df = data.frame(iteration = iteration,
                  seed = seed,
                  
                  region = c("Overall",paste0("Region_",1:nregion)),
                  
                  # Derived trends
                  bias.trend.derived = c(bias.trend.overall, bias.trend.region.derived),
                  precision.trend.derived = c(precision.trend.overall,precision.trend.region.derived),
                  
                  # Linear trends (no fitted linear trend at national scale)
                  bias.trend.linear = c(NA, bias.trend.region.linear),
                  precision.trend.linear = c(NA,precision.trend.region.linear),
                  
                  n.chains = out$mcmc.info$n.chains,
                  n.samples = out$mcmc.info$n.samples,
                  max.Rhat = max(unlist(out$Rhat),na.rm=TRUE),
                  Rhat.1.1 = mean(unlist(out$Rhat)>1.1,na.rm=TRUE))
  
  df

  if (file.exists(file)){
    results.df = read.csv(file)
   }
  
  results.df = rbind(results.df,df)

  write.csv(results.df, file = file, row.names = FALSE)
}