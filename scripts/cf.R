# Load Library and Data ----
library(here)
library(dplyr)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)
library(tidyr)
library(stringr)
library(raster)
library(terra)
library(stars)
library(ggplot2)
library(viridis)
library(rgeos)
library(ggforce)
library(deldir)

rm(list = ls())

`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','c14data_farmers.RData'))


#Read in elevation and resource data ----
fitness_sf <- read_stars("maps/resources.tiff") %>% 
  st_as_sf() %>% 
  rename("fitness_value" = "resources.tiff")

height_sf <- read_stars("maps/east_narnia4x.tif") %>% 
  st_as_sf() %>% 
  rename("height_value" = "east_narnia4x.tif")

#--------------
#Grid concentrations
f_area_freq  <- plyr::count(f_siteInfo, 'area_id') %>% rename(site_freq = freq)

sq_grid_conc <- sq_grid %>% 
  left_join(f_area_freq, by=join_by(area_id)) %>% 
  mutate(site_freq = case_when((contains_land==TRUE & is.na(site_freq)) ~ 0, is.na(contains_land) ~ NA, !is.na(site_freq) ~ site_freq)) #0 if there are no sites, NA if the square area is in the sea

#--------------  
#Aggregate fitness and height values for the land portion of each square area ----
#Fitness
fitness_sf$area_id <- as.integer(st_within(fitness_sf$geometry, sq_grid$land_within_poly)) #Normalised based on the % of land cover vs. ocean cover. Using sq_grid$land_within_poly rather than sq_grid$geometry (i.e. the whole square area) effectively normalises for the portion of land that is contained within the square. Do not take sea values into consideration as human's can't settle in these areas.

fitness_df <- fitness_sf %>% 
  group_by(area_id) %>% 
  summarize(sq_fitness_value = mean(fitness_value)) %>% 
  as.data.frame() %>% 
  dplyr::select(area_id, sq_fitness_value) %>% 
  na.omit()

#Height
height_sf$area_id <- as.integer(st_within(height_sf$geometry, sq_grid$land_within_poly)) #TODO: check why there are so many NAs generated here... 

height_df <- height_sf %>% 
  group_by(area_id) %>% 
  summarize(sq_height_value = mean(height_value)) %>% 
  as.data.frame() %>% 
  dplyr::select(area_id, sq_height_value) %>% 
  na.omit()

#-----
#Join area-level data
sq_grid_conc <- sq_grid_conc %>% 
  left_join(fitness_df, by=join_by(area_id)) %>% 
  left_join(height_df, by=join_by(area_id)) %>% 
  mutate(lat=st_coordinates(area_center)[,1], 
         long=st_coordinates(area_center)[,2])

sq_grid_conc_sp <- sq_grid_conc %>% as('Spatial')


#-------------
## Compute Great-Arc Distances in km between area centers ---
#We take the center of each area k to be a point representing that whole area
sq_area_centers <- sq_grid$area_center %>% as('Spatial')
sq_dist_mat <- spDists(sq_area_centers, longlat=TRUE) #Inter-area distance matrix: each area's distance from every other area.
sq_dist <- spDistsN1(sq_area_centers[1], sq_area_centers[2], longlat=TRUE) #distance between centroids of two adjacent squares

#The square area containing our putative origin
origin_sf <- as(f_sites, 'sf') %>% 
  filter(siteID == f_constants$origin_siteID)
origin_sf$area_id <- as.integer(st_within(origin_sf$geometry, sq_grid$geometry))  
sq_origin_area_center <- sq_area_centers[origin_sf$area_id] #The center of the square containing our putative origin

#-------------
#Delaunay triangulation 
del <- deldir(sq_area_centers@coords, id=sq_grid$area_id)
tiles <- tile.list(del)

transition_matrix <- del$delsgs %>% 
  mutate(height_diff = 
           height_df$sq_height_value[match(ind2, height_df$area_id)] - height_df$sq_height_value[match(ind1, height_df$area_id)],
         fitness_diff = 
           fitness_df$sq_fitness_value[match(ind2, fitness_df$area_id)] - fitness_df$sq_fitness_value[match(ind1, fitness_df$area_id)])

#-------------
## Plot delaunay triangulation
ggplot(data = sq_grid_conc) +
  geom_sf(data = as(sample_win_sp, 'sf')) +
  geom_sf(alpha=0.5, aes(fill=site_freq)) + #Alternatively: fill=cf #fill=sq_fitness_value #fill=sq_height_value #fill=site_freq
  geom_sf(data = as(f_sites, 'sf'), size=3, alpha=0.5, aes(colour="purple")) + #sites
  geom_sf_label(aes(label = area_id), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
  geom_delaunay_segment(aes(x=sq_area_centers@coords[,1], y=sq_area_centers@coords[,2]), 
                        alpha=0.5, 
                        colour='purple',
                        size=0.8) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis()


##Plot edges with fitness and height values 
#Fitness --
ggplot() +
  geom_sf(data = as(sample_win_sp, 'sf')) + #map
  geom_delaunay_segment(data = sq_grid_conc, aes(lat, long, colour=sq_fitness_value, group=-1L),
                        alpha=0.5, 
                        size=2.2) +
  #geom_sf(data = as(sample_win_sp, 'sf')) +
  labs(x = "Longitude", y = "Latitude") +
  scale_color_viridis()

#Height --
ggplot() +
  geom_sf(data = as(sample_win_sp, 'sf')) + #map
  geom_delaunay_segment(data = sq_grid_conc, aes(lat, long, colour=sq_height_value, group=-1L),
                        alpha=0.5, 
                        size=2.2) +
  #geom_sf(data = as(sample_win_sp, 'sf')) +
  labs(x = "Longitude", y = "Latitude") +
  scale_color_viridis()




#===============================================================================
##FRICTION CALCULATION
#Distance between the center of every square area and our origin area center
sq_grid_conc <- sq_grid_conc %>% 
  mutate(sq_dist_org = spDistsN1(sq_area_centers, sq_origin_area_center, longlat=TRUE)) #distance from center of square containing Fallsail site

#------------
## Compute difference in fitness value between square areas and our origin area ---
#Origin square area fitness value
sq_origin_fitness <- fitness_df %>% 
  filter(area_id==origin_sf$area_id)
sq_origin_fitness <-sq_origin_fitness$sq_fitness_value

#Difference in fitness from origin
sq_grid_conc <- sq_grid_conc %>% 
  mutate(sq_fit_org = sq_fitness_value - sq_origin_fitness) #difference in fitness from center of square containing Fallsail site

#------------
## Compute difference in elevation in meters between square areas and our origin area ---
#Origin square area height value
sq_origin_elev <- height_df %>% 
  filter(area_id==origin_sf$area_id)
sq_origin_elev <-sq_origin_elev$sq_height_value

#Difference in elevation from origin
sq_grid_conc <- sq_grid_conc %>% 
  mutate(sq_elev_org = sq_height_value - sq_origin_elev) #difference in elevation from center of square containing Fallsail site  


#-----------
#Friction of square area k from origin square area
#Initially: cf = sq_dist_org - w_f*sq_fit_org + w_e*sq_elev_org
#The signs are based on the following initial assumptions:
#(i) Moving to lower elevation is good -- leads to a decrease in friction
#(ii) Moving to higher fitness values is good -- leads to a decrease in friction
#The weights, w_f and w_e, currently ensure that sq_fit_org and sq_eleve_org are of the same order of magnitude
#TODO: explore these assumptions and give weights to each of these parameters based on the relative importance

w_e <- 1 #To refine. Arrived at by looking at (max(sq_grid_conc$sq_elev_org, na.rm=TRUE)-min(sq_grid_conc$sq_elev_org, na.rm=TRUE)) compared to (max(sq_grid_conc$sq_dist_org, na.rm=TRUE)-min(sq_grid_conc$sq_dist_org, na.rm=TRUE))
w_f <- 600 #To refine. Arrived at by looking at (max(sq_grid_conc$sq_fit_org, na.rm=TRUE)-min(sq_grid_conc$sq_fit_org, na.rm=TRUE)) compared to (max(sq_grid_conc$sq_dist_org, na.rm=TRUE)-min(sq_grid_conc$sq_dist_org, na.rm=TRUE))

sq_grid_conc <- sq_grid_conc %>% 
  mutate(cf = sq_dist_org - w_f*sq_fit_org + w_e*sq_elev_org) %>% 
  mutate(cf = case_when(cf<0 ~ 50, cf>=0 ~ cf)) #If friction is negative -- assign a small positive value. It takes a minimum amount of effort to move. #TODO: adjust amount. Currently chosen at random.


#-------------------------------------------------------------------------------
## Store Output ----
save(sq_grid_conc, file=here('output', 'friction_originaldata.RData')) #save as friction.RData if all data is used (i.e. after the excavation of additional squares)












