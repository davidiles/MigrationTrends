

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'ebirdst' , 'raster' , 'velox' , 'sf', 'smoothr', 'rnaturalearth',
  'dplyr', 'tidyr', 'stringr', 'ggplot2',
  
  'fields'
)

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

# resolve namespace conflicts
select <- dplyr::select

rm(list=ls())

sp_path <- ebirdst_download(species = "Barn Swallow")

# load the abundance data
# this automatically labels layers with their dates (load_raster is an ebird command)
abd <- load_raster("abundance_umean", path = sp_path)

mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
ne_scale <- 50

#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
# Data download
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------

# land polygon
ne_land <- ne_countries(scale = ne_scale, returnclass = "sf") %>%
  filter(continent %in% c("North America", "South America")) %>%
  st_set_precision(1e6) %>%
  st_union() %>% 
  st_geometry()

# function to subset other features to those  within this land area
wh_subset <- function(x) {
  in_wh <- as.logical(st_intersects(x, ne_land, sparse = FALSE))
  st_transform(x[in_wh], crs = mollweide)
}

# country lines
ne_country_lines <- ne_download(scale = ne_scale, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()

# state lines
ne_state_lines <- ne_download(scale = ne_scale, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()

# rivers
ne_rivers <- ne_download(scale = ne_scale, category = "physical",
                         type = "rivers_lake_centerlines",
                         returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()

# lakes
ne_lakes <- ne_download(scale = ne_scale, category = "physical",
                        type = "lakes",
                        returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()

ne_land <- st_transform(ne_land, crs = mollweide)

#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
# Seasonal abundance
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------

# subset to the yellow-bellied sapsucker season definitions
belvir_dates <- filter(ebirdst_runs, species_code == "belvir") %>% 
  
  # just keep the seasonal definition columns
  select(setdiff(matches("(start)|(end)"), matches("year_round"))) %>% 
  
  # transpose
  gather("label", "date") %>% 
  
  # spread data so start and end dates are in separate columns
  separate(label, c("season", "start_end"), "_(?=s|e)") %>% 
  spread(start_end, date) %>% 
  select(season, start_dt, end_dt)

# did the season pass review
belvir_dates <- mutate(belvir_dates, pass = !(is.na(start_dt) | is.na(end_dt)))
belvir_dates


#-------------------------------------------------------
# Use these season definitions to assign each of the weekly abundance layers to a season
#-------------------------------------------------------

# dates for each abundance layer
weeks <- parse_raster_dates(abd)

# assign to seasons
weeks_season <- rep(NA_character_, length(weeks))
for (i in seq_len(nrow(belvir_dates))) {
  s <- belvir_dates[i, ]
  
  # skip seasona assignment if season failed
  if (!s$pass) {
    next()
  }
  
  # handle seasons cross jan 1 separately
  if (s$start_dt <= s$end_dt) {
    in_season <- weeks >= s$start_dt & weeks <= s$end_dt
  } else {
    in_season <- weeks >= s$start_dt | weeks <= s$end_dt
  }
  weeks_season[in_season] <- s$season
}
table(weeks_season)

#-------------------------------------------------------
# Average all the weeks within each season to produce seasonal relative abundance rasters.
#-------------------------------------------------------

# drop weeks not assigned to season
week_pass <- !is.na(weeks_season)
abd <- abd[[which(week_pass)]]
weeks <- weeks[week_pass]
weeks_season <- weeks_season[week_pass]

# average over weeks in season
mean_season <- function(s) {
  calc(abd[[which(weeks_season == s)]], mean, na.rm = TRUE)
}

seasons <- unique(weeks_season)

abd_season <- lapply(seasons, mean_season) %>% 
  stack() %>% 
  setNames(seasons)
abd_season

#-------------------------------------------------------
# ABUNDANCE MAPS
#-------------------------------------------------------

migration_threshold <- 0.4

mig_seasons <- c("prebreeding_migration", "postbreeding_migration")
if (all(mig_seasons %in% names(abd_season))) {
  
  # identify areas with abundance in only one season
  abd_nz <- abd_season[[mig_seasons]] > 0
  
  just_pre <- mask(abd_nz[["prebreeding_migration"]],
                   abd_nz[["postbreeding_migration"]], 
                   maskvalue = 1)
  
  just_post <- mask(abd_nz[["postbreeding_migration"]],
                    abd_nz[["prebreeding_migration"]], 
                    maskvalue = 1)
  
  # count the number of cells with abundance in only one season
  n_just <- cellStats(stack(just_pre, just_post), sum)
  n_all <- cellStats(abd_nz, sum)
  
  # is the proportion of one season cells above the 40% threshold
  split_migration <- max(n_just / n_all, na.rm = TRUE) >= migration_threshold
} else {
  split_migration <- FALSE
}

n_just / n_all

split_migration

#-----------------------------------------

threshold_yearround <- 0.01

# decide whether to show year-round layer
if (nlayers(abd_season) == 4) {
  
  # annual abundance
  abd_yr <- calc(abd, fun = mean, na.rm = TRUE)
  
  # mask out cells that aren't occupied year-round
  year_round <- calc(abd_season > 0, fun = sum, na.rm = TRUE) == 4
  abd_yr_mask <- mask(abd_yr, year_round, maskvalue = 0)
  
  # determine proportion of celss that are occupied year round
  n_yr <- cellStats(abd_yr_mask > 0, sum)
  n_an <- cellStats(abd_yr > 0, sum)
  
  # only show year round abundance if it's above 1% of range threshold
  show_yearround <- ((n_yr / n_an) >= threshold_yearround)
  
} else {
  show_yearround <- FALSE
}

show_yearround

#------------------------------------------------
bin_breaks <- calc_bins(abd_season)
#------------------------------------------------

#-----------------------------------------
#Now that everything is in place, we can actually make the seasonal relative abundance map! Note that the Status and Trends maps distinguish between regions of zero abundance and regions with no prediction, which are displayed in different shades of gray, and regions with non-zero abundance, shown in color. To account for this, we'll need to generate a raster layer delineating the region within which predictions were made.

# project the abundance data to mollweide
# use nearest neighbour resampling to preserve true zeros
abd_season_proj <- projectRaster(abd_season, crs = mollweide, method = "ngb")

# determine spatial extent for plotting
ext <- calc_full_extent(abd_season_proj)

# set the plotting order of the seasons
season_order <- c("postbreeding_migration", "prebreeding_migration", 
                  "nonbreeding", "breeding")

# prediction region, cells with predicted value in at least one week
pred_region <- calc(abd_season_proj, mean, na.rm = TRUE)

# mask to land area
ne_land_buffer <- st_buffer(ne_land, dist = max(res(pred_region)) / 2)
pred_region <- mask(pred_region, as_Spatial(ne_land_buffer))

# remove zeros from abundnace layers
abd_no_zero <- subs(abd_season_proj, data.frame(from = 0, to = NA), 
                    subsWithNA = FALSE)

# set up plot area
par(mar = c(0 , 0, 0, 0))
plot(ne_land, col = "#eeeeee", border = NA, 
     xlim = c(ext@xmin, ext@xmax),
     ylim = c(ext@ymin, ext@ymax))

# prediction region and explicit zeros
plot(pred_region, col = "#dddddd", maxpixels = raster::ncell(pred_region),
     legend = FALSE, add = TRUE)

# lakes
plot(ne_lakes, col = "#ffffff", border =  "#444444", lwd = 0.5, add = TRUE)

# land border
plot(ne_land, col = NA, border = "#444444", lwd = 0.5, add = TRUE)

# seasonal layer
plot_seasons <- intersect(season_order, names(abd_no_zero))
for (s in plot_seasons) {
  # handle splitting of migration seasons into different colors
  if (!split_migration && s %in% c("prebreeding_migration", 
                                   "postbreeding_migration")) {
    pal_season <- "migration"
    
  } else {
    pal_season <- s
  }
  pal <- abundance_palette(length(bin_breaks$bins) - 1, pal_season)
  plot(abd_no_zero[[s]], col = pal, breaks = bin_breaks$bins,
       maxpixels = ncell(abd_no_zero[[s]]),
       legend = FALSE, add = TRUE)
}

# year round
if (show_yearround) {
  year_round_proj <- projectRaster(year_round, crs = mollweide, method = "ngb")
  plot(year_round_proj, 
       col = abundance_palette(length(bin_breaks$bins) - 1, "year_round"), 
       breaks = bin_breaks$bins,
       maxpixels = ncell(year_round_proj),
       legend = FALSE, add = TRUE)
}

# linework
plot(ne_rivers, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 2, add = TRUE)

# legends
legend_seasons <- plot_seasons
if (split_migration) {
  legend_seasons[legend_seasons %in% c("prebreeding_migration", 
                                       "postbreeding_migration")] <- "migration"
  legend_seasons <- unique(legend_seasons)
}
if (show_yearround) {
  legend_seasons <- c(legend_seasons, "year_round")
}

# thin out labels
lbl_at <- bin_breaks$bins^bin_breaks$power
lbl_at <- c(min(lbl_at), median(lbl_at), max(lbl_at))
lbl <- lbl_at^(1 / bin_breaks$power)
lbl <- format(round(lbl, 2), nsmall = 2)

# plot legends
for (i in seq_along(legend_seasons)) {
  pal <- abundance_palette(length(bin_breaks$bins) - 1, legend_seasons[i])
  if (i == 1) {
    axis_args <- list(at = lbl_at, labels = lbl, line = -1,
                      cex.axis = 0.75, lwd = 0)
  } else {
    axis_args <- list(at = lbl_at, labels = rep("", 3),
                      cex.axis = 0.75, lwd = 0)
  }
  legend_title <- legend_seasons[i] %>% 
    str_replace_all("_", " ") %>% 
    str_to_title()
  fields::image.plot(zlim = range(bin_breaks$bins^bin_breaks$power), 
                     legend.only = TRUE, 
                     breaks = bin_breaks$bins^bin_breaks$power, col = pal,
                     smallplot = c(0.05, 0.35, 0.01 + 0.06 * i, 0.03 + 0.06 * i),
                     horizontal = TRUE,
                     axis.args = axis_args,
                     legend.args = list(text = legend_title, side = 3, 
                                        cex = 0.9, col = "black", line = 0.1))
}

title("Barn Swallow Relative Abundance (birds per km/hr)", 
      line = -1, cex.main = 1)
