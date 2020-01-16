# library(ggplot2)
# library(jagsUI)
# library(reshape2)
# library(bbsBayes)
# library(dplyr)
# library(ggrepel)


# Required packages
my_packs <- c(
  
  # For data management and plotting
  'ebirdst' , 'raster' , 'velox' , 'sf', 'smoothr', 'rnaturalearth',
  'dplyr', 'tidyr', 'stringr', 'ggplot2',
  
  'fields',
  
  'bbsBayes','jagsUI','ggprepel'
)


rm(list=ls())

#Function that determines when the station is open
open.fn = function(dat, sp){
  
  site.summary = data.frame()
  
  for (L in unique(dat$Locality)){
    
    L.dat = subset(dat, Locality == L)
    if (sp %in% L.dat$CommonName == FALSE) next #Skip if target species is never observed in this Locality
    
    locality.summary = expand.grid(JulianDay = seq(min(L.dat$JulianDay),max(L.dat$JulianDay)),
                                   YearCollected = seq(min(L.dat$YearCollected),max(L.dat$YearCollected)))
    
    # Calculate number of species with non-zero counts on each day of each year
    number.species.counted = aggregate(CommonName~JulianDay + YearCollected,
                                       data = subset(L.dat, ObservationCount >= 1 & !is.na(ObservationCount)),
                                       FUN = function(x) length(unique(x)))
    colnames(number.species.counted)[3] = "Number.of.Species.Counted"
    locality.summary = merge(locality.summary, number.species.counted, all = TRUE)
    
    # Extract counts of focal species
    focal.species.counted = subset(L.dat, ObservationCount >= 1 & !is.na(ObservationCount) & CommonName == sp, select = c("JulianDay","YearCollected","ObservationCount"))
    colnames(focal.species.counted)[3] = "Focal.Species.Counted"
    locality.summary = merge(locality.summary, focal.species.counted, all = TRUE)
    
    # Assume site is open if >= 5 species were counted on that day
    locality.summary$Open = !is.na(locality.summary$Number.of.Species.Counted) & locality.summary$Number.of.Species.Counted >= 5
    locality.summary$Focal.Species.Counted[which(locality.summary$Open == FALSE)] = NA
    
    locality.summary$Focal.Species.Counted[which(locality.summary$Open == TRUE & is.na(locality.summary$Focal.Species.Counted))] = 0
    
    locality.summary$Locality = L
    locality.summary$SurveyAreaIdentifier = L.dat$SurveyAreaIdentifier[1]
    
    locality.summary$DecimalLatitude = L.dat$DecimalLatitude[1]
    locality.summary$DecimalLongitude = L.dat$DecimalLongitude[1]
    locality.summary$ObservationDate = L.dat$ObservationDate[1]
    
    
    locality.summary$Focal.Species = sp
    
    site.summary = rbind(site.summary, locality.summary)
    
  }
  
  return(site.summary)
}

focal.spec = "Blackpoll Warbler"
#focal.spec = "Wilson's Warbler"

#-----------------------------------------------------------
# Read and plot BBS data
#-----------------------------------------------------------
strat_data <- stratify(by = "bbs_cws")

jags_data <- prepare_jags_data(strat_data, species_to_run = focal.spec, model = "firstdiff")

#Run model and plot results
# mod <- run_model(jags_data = jags_data,n_burnin = 5000,n_iter=15000)
# 
# strat_indices <- generate_strata_indices(mod)
# strat_trend <- generate_strata_trends(indices = strat_indices)
# generate_map(strat_trend, stratify_by = "bbs_cws")
# 
# cont_indices <- generate_cont_indices(mod)
# plot_cont_indices(indices = cont_indices)

# BBS records of BLPW
route.info = unique(strat_data$route_strat[,c("rt.uni","Latitude","Longitude","strat_name")])

blpw.bbs = data.frame(count = jags_data$count,route = jags_data$route,year = jags_data$r_year)
blpw.bbs = merge(blpw.bbs, route.info, by.x = "route", by.y = "rt.uni", all.x = TRUE)

#Time series of counts on each route
ggplot(data = blpw.bbs) +
  geom_line(aes(x = year, y = log(count), linetype = route), alpha = 0.3)+
  geom_point(aes(x = year, y = log(count)), alpha = 0.3)+
  scale_linetype_manual(values=rep(1,length(unique(blpw.bbs$route))), guide = FALSE)+
  facet_grid(strat_name~., scales = "free")

#Extract mean count from each route, plot on map
mc = aggregate(count ~ route + Latitude + Longitude + strat_name, data = blpw.bbs, FUN = mean)
WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify

#Number of years where BLPW was detected on each route
ny = aggregate(year ~ route + Latitude + Longitude + strat_name, data = subset(blpw.bbs, count > 0), FUN = length)

p1 = ggplot(data = mc)+
  
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = "white", colour = "#7f7f7f", size=0.5) +
  
  geom_point(data = strat_data$route_strat, aes(x = Longitude, y = Latitude), col = "gray90", size = 1)+
  geom_point(aes(x = Longitude, y = Latitude), col = "dodgerblue", size = 1)+
  coord_cartesian(xlim=c(-175,-50), ylim=c(40,75))+
  xlab("Lon")+ylab("Lat")+
  ggtitle("BBS Routes that detected BLPW")+
  theme_bw()

#print(p1)

pdf("../figures/data_summaries/BLPW_bbs_mean.pdf", width = 8, height = 6)
print(p1)
dev.off()


# BBS Routes that detected BLPW in more than 5 years
p2 = ggplot(data = subset(ny, year > 5))+
  
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = "white", colour = "#7f7f7f", size=0.5) +
  
  geom_point(data = strat_data$route_strat, aes(x = Longitude, y = Latitude), col = "gray90", size = 1)+
  geom_point(aes(x = Longitude, y = Latitude), col = "dodgerblue", size = 1)+
  coord_cartesian(xlim=c(-175,-50), ylim=c(40,75))+
  xlab("Lon")+ylab("Lat")+
  ggtitle("BBS Routes that detected BLPW in more than 5 years")+
  theme_bw()

print(p2)

pdf("../figures/data_summaries/BLPW_bbs_ny5.pdf", width = 8, height = 6)
print(p2)
dev.off()

# #-----------------------------------------------------------
# # All code below this point is deprecated.  It is retained in this script
# # in case pieces of it are useful for subsequent analyses (e.g., to 
# # implement the same site selection criteria as in Crewe et al.)
# #-----------------------------------------------------------
# 

# #-----------------------------------------------------------
# # Read and format CMMN data
# #-----------------------------------------------------------
# 
# #Read in all the stations
# acbo = read.delim("../data/acbo_data.txt") #Albert creek
# bbo = read.delim("../data/bbo_data.txt") #Beaverhill bird observatory
# bpbo = read.delim("../data/bpbo_data.txt") #Bruce peninsula
# hbo = read.delim("../data/hbo_data.txt") #Haldimand
# lpbo = read.delim("../data/lpbo_data.txt") #Long point
# lslbo = read.delim("../data/lslbo_data.txt") #Lesser slave lake
# mgbo = read.delim("../data/mgbo_data.txt") #McGill Bird Observatory
# mibo = read.delim("../data/mibo_data.txt") #McKellar Island Bird Observatory
# mno = read.delim("../data/mno_data.txt") #Mackenzie Nature Observatory
# tcbo = read.delim("../data/tcbo_data.txt") #Thunder cape
# tlbbs = read.delim("../data/tlbbs_data.txt") #Teslin Lake Bird Banding Station
# ttpbrs = read.delim("../data/ttpbrs_data.txt") #Tommy Thompson Park Bird Research Station
# 
# min.day = 0
# max.day = 360
# min.year = 1980
# max.year = 2016
# 
# #Albert creek
# acbo.fs = open.fn(acbo, sp = focal.spec)
# acbo.fs$Site = "acbo"
# acbo.fs = subset(acbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# acbo.fs$locality.number = as.numeric(factor(as.character(acbo.fs$Locality)))
# 
# #Beaverhill Bird Observatory
# bbo.fs = open.fn(bbo, sp = focal.spec)
# bbo.fs$Site = "bbo"
# bbo.fs = subset(bbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# bbo.fs$locality.number = as.numeric(factor(as.character(bbo.fs$Locality)))
# 
# #Bruce peninsula
# bpbo.fs = open.fn(bpbo, sp = focal.spec)
# bpbo.fs$Site = "bpbo"
# bpbo.fs = subset(bpbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# bpbo.fs$locality.number = as.numeric(factor(as.character(bpbo.fs$Locality)))
# 
# #Haldimand
# hbo.fs = open.fn(hbo, sp = focal.spec)
# hbo.fs$Site = "hbo"
# hbo.fs = subset(hbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# hbo.fs$locality.number = as.numeric(factor(as.character(hbo.fs$Locality)))
# 
# #Long point
# lpbo.fs = open.fn(lpbo, sp = focal.spec)
# lpbo.fs$Site = "lpbo"
# lpbo.fs = subset(lpbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# lpbo.fs$locality.number = as.numeric(factor(as.character(lpbo.fs$Locality)))
# 
# #Lesser slave lake
# lslbo.fs = open.fn(lslbo, sp = focal.spec)
# lslbo.fs$Site = "lslbo"
# lslbo.fs = subset(lslbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# lslbo.fs$locality.number = as.numeric(factor(as.character(lslbo.fs$Locality)))
# 
# #McGill Bird Observatory
# mgbo.fs = open.fn(mgbo, sp = focal.spec)
# mgbo.fs$Site = "mgbo"
# mgbo.fs = subset(mgbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# mgbo.fs$locality.number = as.numeric(factor(as.character(mgbo.fs$Locality)))
# 
# #McKellar Island Bird Observatory
# mibo.fs = open.fn(mibo, sp = focal.spec)
# mibo.fs$Site = "mibo"
# mibo.fs = subset(mibo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# mibo.fs$locality.number = as.numeric(factor(as.character(mibo.fs$Locality)))
# 
# #Mackenzie Nature Observatory
# mno.fs = open.fn(mno, sp = focal.spec)
# mno.fs$Site = "mno"
# mno.fs = subset(mno.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# mno.fs$locality.number = as.numeric(factor(as.character(mno.fs$Locality)))
# 
# #Thunder cape
# tcbo.fs = open.fn(tcbo, sp = focal.spec)
# tcbo.fs$Site = "tcbo"
# tcbo.fs = subset(tcbo.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# tcbo.fs$locality.number = as.numeric(factor(as.character(tcbo.fs$Locality)))
# 
# #Teslin Lake Bird Banding Station
# tlbbs.fs = open.fn(tlbbs, sp = focal.spec)
# tlbbs.fs$Site = "tlbbs"
# tlbbs.fs = subset(tlbbs.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# tlbbs.fs$locality.number = as.numeric(factor(as.character(tlbbs.fs$Locality)))
# 
# #Tommy Thompson Park Bird Research Station
# ttpbrs.fs = open.fn(ttpbrs, sp = focal.spec)
# ttpbrs.fs$Site = "ttpbrs"
# ttpbrs.fs = subset(ttpbrs.fs, JulianDay >= min.day & JulianDay <= max.day & YearCollected >= min.year & YearCollected <= max.year & Open == TRUE)
# ttpbrs.fs$locality.number = as.numeric(factor(as.character(ttpbrs.fs$Locality)))
# 
# 
# 
# #Combine all data into a single dataframe
# dat = rbind(acbo.fs,
#             bbo.fs,
#             bpbo.fs,
#             hbo.fs,
#             lpbo.fs,
#             lslbo.fs,
#             mgbo.fs,
#             mibo.fs,
#             mno.fs,
#             tcbo.fs,
#             tlbbs.fs,
#             ttpbrs.fs)
# 
# dat$Site.Locality = paste(dat$Site,dat$locality.number,sep=".")
# 
# #Lat Long of each station
# site.locs = unique(dat[,c("Locality","Site.Locality","DecimalLatitude","DecimalLongitude")])
# 
# #mean daily spring and fall counts at each CMMN station
# cmmn.spring.mean = aggregate(Focal.Species.Counted~Site.Locality, data = subset(dat, JulianDay >= 100 & JulianDay <= 180), FUN = function(x) mean(x, na.rm=TRUE))
# cmmn.fall.mean = aggregate(Focal.Species.Counted~Site.Locality, data = subset(dat, JulianDay >= 180 & JulianDay <= 360), FUN = function(x) mean(x, na.rm=TRUE))
# 
# site.locs = merge(site.locs,cmmn.spring.mean, all.x = TRUE)
# colnames(site.locs)[ncol(site.locs)] = "mean.spring.count"
# site.locs = merge(site.locs,cmmn.fall.mean, all.x = TRUE)
# colnames(site.locs)[ncol(site.locs)] = "mean.fall.count"
# 
# fall.count.plot = ggplot(data = mc)+
#   
#   geom_map(data = WorldData, map = WorldData,
#            aes(x = long, y = lat, group = group, map_id=region),
#            fill = "white", colour = "#7f7f7f", size=0.5) +
#   
#   geom_point(data = strat_data$route_strat, aes(x = Longitude, y = Latitude), col = "gray90")+
#   geom_point(aes(x = Longitude, y = Latitude, size = count), col = "dodgerblue")+
#   
#   geom_point(data = site.locs, aes(x = DecimalLongitude, 
#                                    y = DecimalLatitude, 
#                                    col = mean.fall.count), 
#              shape = 18, size = 4)+
#   geom_label_repel(data = site.locs, 
#                    aes(x = DecimalLongitude, y = DecimalLatitude, label = Site.Locality,col = mean.fall.count), 
#                    force = 6)+
#   
#   scale_color_gradientn(colors = c("black","darkred","goldenrod"))+
#   coord_cartesian(xlim=c(-175,-50), ylim=c(40,75))+
#   ggtitle(paste(focal.spec, "max counts on BBS routes"))+
#   theme_bw()
# 
# pdf("../figures/fig2_cmmn_fallmeans.pdf", width = 16, height = 10)
# print(fall.count.plot)
# dev.off()
# 
# spring.count.plot = ggplot(data = mc)+
#   
#   geom_map(data = WorldData, map = WorldData,
#            aes(x = long, y = lat, group = group, map_id=region),
#            fill = "white", colour = "#7f7f7f", size=0.5) +
#   
#   geom_point(data = strat_data$route_strat, aes(x = Longitude, y = Latitude), col = "gray90")+
#   geom_point(aes(x = Longitude, y = Latitude, size = count), col = "dodgerblue")+
#   
#   geom_point(data = site.locs, aes(x = DecimalLongitude, 
#                                    y = DecimalLatitude, 
#                                    col = mean.spring.count), 
#              shape = 18, size = 4)+
#   geom_label_repel(data = site.locs, 
#                    aes(x = DecimalLongitude, y = DecimalLatitude, label = Site.Locality,col = mean.spring.count), 
#                    force = 6)+
#   
#   scale_color_gradientn(colors = c("black","darkred","goldenrod"))+
#   coord_cartesian(xlim=c(-175,-50), ylim=c(40,75))+
#   ggtitle(paste(focal.spec, "max counts on BBS routes"))+
#   theme_bw()
# 
# pdf("../figures/fig3_cmmn_springmeans.pdf", width = 16, height = 10)
# print(spring.count.plot)
# dev.off()
# 
# #-----------------------------------------------------
# # Restrict analysis to migration stations with suitable data
# #-----------------------------------------------------
# min.day = 180
# max.day = 300
# min.year = 1980
# max.year = 2016
# 
# dat = subset(dat, JulianDay >= min.day & JulianDay <= max.day)
# 
# #To be included in analysis, site must record:
# # more than 10 birds per season on average
# yearly.sums = aggregate(Focal.Species.Counted ~ YearCollected + Locality + Site, data = dat, FUN = sum)
# site.means = aggregate(Focal.Species.Counted ~ Locality + Site, data = yearly.sums, FUN = mean)
# sites.to.keep = subset(site.means, Focal.Species.Counted >= 10)
# dat = subset(dat, Locality %in% sites.to.keep$Locality)
# 
# # detected more than 5 times per season on average
# number.days.observed = aggregate((Focal.Species.Counted > 0) ~ YearCollected + Locality + Site, data = dat, FUN = function(x) sum(x,na.rm=TRUE))
# colnames(number.days.observed)[4] = "Times.Observed"
# site.means = aggregate(Times.Observed ~ Locality + Site, data = number.days.observed, FUN = mean)
# sites.to.keep = subset(site.means, Times.Observed >= 5)
# 
# dat = subset(dat, Locality %in% sites.to.keep$Locality)
# 
# daily.means = aggregate(Focal.Species.Counted ~ JulianDay + Site.Locality, data = dat, FUN = function(x) mean(x,na.rm=TRUE))
# 
# daily.n = aggregate(YearCollected ~ JulianDay + Site.Locality, data = dat, FUN = function(x) length(unique(x)))
# daily.means$Focal.Species.Counted[which(daily.n$YearCollected<=5)] = NA
# 
# dat$site.number = as.numeric(factor(as.character(dat$Site)))
# 
# yearly.daily.counts.plot = ggplot(data = dat)+
#   geom_line(aes(x = JulianDay, y = Focal.Species.Counted), col = "blue")+
#   facet_grid(Site.Locality~YearCollected, scales = "free")
# 
# pdf("../figures/fig4_yearly_daily_counts.pdf", width = length(unique(dat$YearCollected))*5, height = length(unique(dat$Locality))*3)
# print(yearly.daily.counts.plot)
# dev.off()
# 
# yearly.daily.counts.plot = ggplot(data = subset(dat, Site.Locality == "lpbo.4" & YearCollected %in% seq(2005,2016)))+
#   geom_line(aes(x = JulianDay, y = Focal.Species.Counted), col = "blue")+
#   facet_grid(.~YearCollected)+
#   ylab("Number Counted")+
#   xlab("Day")+
#   theme_bw()+
#   theme(axis.text = element_text(size = 5))
# print(yearly.daily.counts.plot)
# 
# yearly.daily.counts.plot = ggplot(data = subset(dat, Site.Locality == "bpbo.1" & YearCollected %in% seq(2005,2016)))+
#   geom_line(aes(x = JulianDay, y = Focal.Species.Counted), col = "blue")+
#   facet_grid(.~YearCollected)+
#   ylab("Number Counted")+
#   xlab("Day")+
#   theme_bw()+
#   theme(axis.text = element_text(size = 5))
# print(yearly.daily.counts.plot)

