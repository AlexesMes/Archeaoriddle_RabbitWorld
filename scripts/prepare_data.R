# Load Libraries and Data ----
library(rcarbon)
library(nimbleCarbon)
library(maptools)
library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(parallel)
library(stars)
library(raster)
library(terra)

rm(list = ls())


#-------------------------------------------------------------------------------
## Read site location data  ----
raw_dat <- read.csv(here("data", "Biblio_data.csv"))

#Extra site location data
raw_dat8 <- read.csv(here("data", "square_8.csv"))
raw_dat12 <- read.csv(here("data", "square_12.csv"), colClasses = c("integer", "character", "numeric", "numeric", "character", "character"))
raw_dat48 <- read.csv(here("data", "square_48.csv"))
raw_dat62 <- read.csv(here("data", "square_62.csv"))
raw_dat98 <- read.csv(here("data", "square_98.csv"))

#Merge data
raw_dat <- rbind(raw_dat, raw_dat8, raw_dat12, raw_dat48, raw_dat62, raw_dat98)

#-------------------------------------------------------------------------------
## Data filtering and cleaning ----

#All dates should be a separate row (observation) in the dataframe
len <- ifelse(grepl("\\|", raw_dat$dates), sapply(gregexpr("\\|", raw_dat$dates), length) + 1, 1)
raw_dat$len <- len
raw_dat <- data.frame(cbind(raw_dat, str_split_fixed(raw_dat$dates, "\\|", max(len))))

dateID_cols <- paste0("X", as.character(1:max(len)))
  
#Pivot data to tidy form
dat <- raw_dat %>% 
  arrange(economy) %>% #For phasemodel useful to have site_ids for Farmers running from 1:f_n_sites
  mutate(siteID = seq(1, nrow(raw_dat), 1)) %>% 
  dplyr::select(-dates, -len, -X) %>% 
  pivot_longer(dateID_cols, names_to = "dateID", values_to = "dates") %>% 
  filter(dates!="") %>% 
  mutate(dateID=row_number())

#Correct column types and separate c14 dates and errors
dat <- dat %>% 
  separate_wider_delim(dates, "\u00B1", names=c("c14age", "c14error")) %>% #where \u00B1 is the unicode for the '+-' symbol
  mutate(economy = as.factor(economy), 
         c14age = as.numeric(c14age), 
         c14error=as.numeric(str_remove(c14error, "BP")))

#-------------------------------------------------------------------------------
## Determining which calibration curve should be used----

# Assign calibration curve ----
dat$calCurve <- ifelse((dat$lat>=0), 'intcal20', 'shcal20') #Assign the calibration curve to use based on the site's position relative to the equator 
#table(dat$calCurve)

#-------------------------------------------------------------------------------
# Restructure Data for Bayesian Analyses ----

# Compute median calibrated dates ----
dat$median_dates = medCal(calibrate(dat$c14age, 
                                    dat$c14error, 
                                    calCurve=dat$calCurve))

# Collect site level information ----
earliest_dates <- aggregate(median_dates~siteID, data=dat, FUN=max) #Earliest medCal Date for Each Site
latest_dates <- aggregate(median_dates~siteID, data=dat, FUN=min) #Latest medCal Date for Each Site
n_dates <- aggregate(median_dates~siteID, data=dat, FUN=length) #Number of medCal Date for Each Site (equivalent to len column previously examined)

siteInfo <- data.frame(siteID = earliest_dates$siteID,
                       earliest = earliest_dates$median_dates,
                       latest = latest_dates$median_dates,
                       diff = earliest_dates$median_dates - latest_dates$median_dates,
                       n_dates = n_dates$median_dates) %>% unique()

siteInfo <- siteInfo %>% 
  left_join(unique(dplyr::select(dat, 
                                 lat, 
                                 lon, 
                                 siteID,
                                 sitename,
                                 economy,
                                 calCurve)))

#------
#Filter for farming sites
f_siteInfo <- siteInfo %>% filter(economy=="F")
#------

# Collect date level information ----
dateInfo <- unique(dplyr::select(dat,
                          dateID,
                          siteID,
                          economy,
                          cra=c14age,
                          cra_error=c14error,
                          median_dates=median_dates,
                          calCurve=calCurve)) %>% arrange(dateID) 
dateInfo$earliestAtSite  <- FALSE #initialize 

for (i in unique(siteInfo$siteID))
{
  tmp_index  <-  which(dateInfo$siteID==i) #Identify the row indexes of all sites that share siteID x
  ii  <- tmp_index[which.max(dateInfo$median_dates[tmp_index])] #Select the index of the row with the maximum median date (earliest date) within the sub-group of rows selected in the previous line (given by tmp.index)
  dateInfo$earliestAtSite[ii]  <- TRUE #Designate this observation as the earliest date
} 
#table(dateInfo$earliestAtSite) #Check

#------
#Filter for farming sites
f_dateInfo <- dateInfo %>% filter(economy=="F")


#-------------------------------------------------------------------------------
## Designating approximate origin of poppychewers ----

#Since we will be tracking the dispersal of the poppychewers (F), and we have no additional archaeological knowledge to rely on, we assume the approximate origin as the oldest site we have access to
f_possible_origin_dat <- dat %>%
  filter(economy=="F") %>% 
  arrange(desc(c14age)) %>%
  slice(1L)

f_origin_siteID <- f_possible_origin_dat$siteID

#-------------------------------------------------------------------------------
## Compute Great-Arc Distances in km for farming sites ---- #TODO: could there be a way to incorporate height values and/or environmental fitness into these calculations of distances?
#Poppychewers
f_sites <- f_siteInfo
coordinates(f_sites) <- c('lon','lat')
proj4string(f_sites)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
f_dist_mat  <- spDists(f_sites, longlat=TRUE) #inter-site distance matrix: each site's distance from every other site (i.e. with n sites, this matrix is n^2)
f_origin_point  <- c(f_possible_origin_dat$lon, f_possible_origin_dat$lat) #Origin: Fallsail site
f_dist_org  <-  spDistsN1(f_sites, f_origin_point, longlat=TRUE) #distance from Fallsail site

#All sites
sites <- siteInfo
coordinates(sites) <- c('lon','lat')
proj4string(sites)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#-------------------------------------------------------------------------------
##Sample window: Rabbithole

#Read in basemap ----
rabbithole_map <- rast(here('maps', 'resources.tiff'))

#Define sample window ----
boundary_win <- rabbithole_map %>% st_bbox() 
sample_win_sf <- as.polygons(rabbithole_map > -Inf) %>% st_as_sf()
sample_win_sp <- sample_win_sf %>% as("Spatial")

#Divide sample window into square polygons
sq_grid <- st_make_grid(sample_win_sp, cellsize = 1, square = TRUE) %>% st_sf() #cellsize=0.5 makes 100 sq areas

#Check if the grid contains land (alternatively only sea) ----
contains_land <- as.logical(st_intersects(sq_grid$geometry, sample_win_sf$geometry))
contains_land[is.na(contains_land)] <- FALSE

#Proportion of land contained in an area - will be used to normalise area-average parameter values  ----
area_sq <- max(st_area(sq_grid)) #Area of 1 square grid block (in m^2) completely covered by land
land_within_sq <- st_intersection(sq_grid$geometry, sample_win_sf$geometry) #Creates geometry of the portion of the land within each square
#Check: ggplot(data = land_within_sq) + geom_sf() 
sq_grid <- sq_grid %>% mutate(land_within_poly = land_within_sq)#save for later
prop_land <- st_area(land_within_sq)/area_sq #proportion of each square area that is covered by land


#Assign square IDs ----
sq_grid <- sq_grid %>%
  mutate(area_id = row_number(),
         area_center = st_centroid(geometry),
         contains_land = contains_land,
         prop_land = as.numeric(prop_land))

#--------
##Assign hex area id to each site ----
#Poppychewers
f_sites_sf <- as(f_sites, 'sf')
f_siteInfo$area_id <- as.integer(st_within(f_sites_sf$geometry, sq_grid$geometry))

#All sites
sites_sf <- as(sites, 'sf')
siteInfo$area_id <- as.integer(st_within(sites_sf$geometry, sq_grid$geometry))

##CHECK ---
# area_freq  <- plyr::count(siteInfo, 'area_id')
# f_area_freq  <- plyr::count(f_siteInfo, 'area_id')
# To check that this lines up visually with how many sites are in each hex area -- see map_fig_grid

#-------------------------------------------------------------------------------  
## Create list with constants and data ----

# Data ----
#Poppychewers date information
f_dat <- list(cra=f_dateInfo$cra, cra_error=f_dateInfo$cra_error)
#All date information
dat <- list(cra=dateInfo$cra, cra_error=dateInfo$cra_error, economy=dateInfo$economy) 

# Constants ----
data(intcal20)
data(shcal20)

#Poppychewers
f_constants <- list()
f_constants$n_sites <- nrow(f_siteInfo)
f_constants$n_dates  <- nrow(f_dateInfo)
f_constants$id_sites <- f_dateInfo$siteID
f_constants$dist_mat  <- f_dist_mat
f_constants$dist_org  <- f_dist_org
f_constants$n_areas  <- nrow(sq_grid) #All areas (even empty ones) are included #Only occupied areas: length(unique(f_siteInfo$area_id))
f_constants$id_area  <- f_siteInfo$area_id 
f_constants$origin_point <- f_origin_point
f_constants$origin_siteID <- f_origin_siteID
#Calibration curves
f_constants$calBP <- intcal20$CalBP #Same for intcal20 and shcal20 
f_constants$C14BP  <- cbind(intcal20$C14Age, shcal20$C14Age) #Northern and southern hemisphere calibration curves
f_constants$C14err  <- cbind(intcal20$C14Age.sigma, shcal20$C14Age.sigma)

#All sites
constants <- list()
constants$n_sites <- nrow(siteInfo)
constants$n_dates  <- nrow(dateInfo)
constants$id_sites <- dateInfo$siteID
constants$n_areas  <- nrow(sq_grid)
constants$id_area  <- siteInfo$area_id 
#Calibration curves
constants$calBP <- intcal20$CalBP #Same for intcal20 and shcal20 
constants$C14BP  <- cbind(intcal20$C14Age, shcal20$C14Age) #Northern and southern hemisphere calibration curves
constants$C14err  <- cbind(intcal20$C14Age.sigma, shcal20$C14Age.sigma)

#-------------------------------------------------------------------------------  
## Save everything on a R image file ----
# Save poppychewer data
save(f_sites, f_constants, f_dat, f_siteInfo, f_dateInfo, sample_win_sp, sq_grid, file=here('data','c14data_farmers.RData'))

# Save poppychewer and rabbitskinner data
save(sites, constants, dat, siteInfo, dateInfo, sample_win_sp, sq_grid, file=here('data','c14data.RData'))
