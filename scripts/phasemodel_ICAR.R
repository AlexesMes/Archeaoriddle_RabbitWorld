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

rm(list = ls())

`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','c14data_farmers.RData'))


#===============================================================================
# MCMC RunScript (Model 1 Basic without site interdependence or aggregation up to an area level) ----
# model1 <- nimbleCode({
#   for (i in 1:n_dates)
#   {
#     theta[i] ~ dunif(max=a, min=b); #where a and b are the start and end dates of the focus area (switched around because BP dates go in the positive direction for nimbleCarbon)
#     # Calibration
#     mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #C14 age on the relevant calibration curve
#     sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]); #error on the calibration curve
#     error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2); #the samples' C14 error + error on calibration curve
#     cra[i] ~ dnorm(mean=mu[i], sd=error[i]); #observed radiocarbon age of the sample
#   }
#   a ~ dunif(5000,12000);
#   b ~ dunif(5000,12000);
#   unif.const ~ dconstraint(b<a); #date a is earlier than date b
# })
# 
# 
# #Initialise parameters ---- 
# constants1 <- list(n_dates = f_constants$n_dates,
#                    calBP = intcal20$CalBP,
#                    C14BP = intcal20$C14Age,
#                    C14err = intcal20$C14Age.sigma)
# 
# d1 <- list(cra = f_dateInfo$cra,
#            cra_error = f_dateInfo$cra_error,
#            unif.const=1)
# 
# theta.init = medCal(calibrate(d1$cra, 
#                               d1$cra_error, 
#                               verbose = FALSE))
# 
# inits1 <- list(a=12000,
#                b=7000,
#                theta=theta.init)
# 
# 
# #Run MCMC ----
# mcmc.samples1<- nimbleMCMC(code = model1,
#                            constants = constants1,
#                            data = d1,
#                            niter = 200, 
#                            nchains = 3, 
#                            thin=10, 
#                            nburnin = 100,
#                            monitors=c('a','b','theta'), 
#                            inits=inits1, 
#                            samplesAsCodaMCMC=TRUE)
# 
# #Diagnostics ----
# rhat1  <- gelman.diag(mcmc.samples1, multivariate = FALSE)
# ess1  <- effectiveSize(mcmc.samples1)

#===============================================================================
##MODEL 2 -- Add in areas

# General Setup ----
# Data --
dat <- list(cra = f_dateInfo$cra,
            cra_error = f_dateInfo$cra_error,
            constraint_uniform = rep(1, f_constants$n_areas),
            cra_constraint = rep(1, f_constants$n_dates)) # Set-up constraint for ignoring inference outside calibration range

origin_site <- f_siteInfo %>% filter(siteID==f_constants$origin_siteID)

# Initial parameters --
buffer <- 100
theta_init <- f_dateInfo$median_dates

# Initialise regional parameters ----
# Initialise sq areas which contain sites
init_a  <- aggregate(earliest~area_id, FUN=max, data=f_siteInfo) #find earliest date in each region k
init_b  <- aggregate(latest~area_id, FUN=min, data=f_siteInfo) #find latest date in each area k

#Initialise sq areas which do not contain sites
init_empty_area <- function(init_df) {
  for(i in 1:f_constants$n_area){
    
    area_ids <- init_df$area_id #List of sq areas ids with sites
    
    if (i %!in% area_ids){
      empty_sq_id <- i #Id of empty hex
      neighbour_sq <- which.min(abs(i - area_ids)) #Determine closest sq neighbor which has sites. If there are more than one neighbour sq area with sites, it selects the first observation (i.e. sq with the smallest id, since ids are in ascending order -- which favours neighbours South West to the sq area in question (proving a good lower bound to age estimate). This makes sense, since our putative origin is in Sq 12 and we think the direction of spread is North East). #TODO: refine this selection of neighbours, once general direction of spread is determined through GPQR
      neighbour_sq_id <- area_ids[neighbour_sq] #Determine area id of closest sq neighbor
      neighbour_date <- init_df[neighbour_sq , 2] #Select the date associated with the neighbour sq
      init_df <- rbind(init_df, c(i, neighbour_date)) #Assign this date to the empty sq area
    }
  }
  return(init_df)
}
init_a <- init_a %>% init_empty_area() %>%  arrange(area_id)
init_b <- init_b %>% init_empty_area() %>%  arrange(area_id)

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

f_constants$id_date_area <- f_constants$id_area[f_constants$id_sites] #indexing determines which area k date i belongs to. 
  
#------------------------------------------------------------
# model2 <- nimbleCode({
#   for (i in 1:n_dates)
#   {
#     theta[i] ~ dunif(max=a[id_date_area[i]], min=b[id_date_area[i]]); # Alternatively k=id_area[id_sites[i]] -- but, dynamic indexing not allowed
#     # Calibration
#     mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #C14 age on the relevant calibration curve
#     sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]); #error on the calibration curve
#     error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2); #the samples' C14 error + error on calibration curve
#     cra[i] ~ dnorm(mean=mu[i], sd=error[i]); #observed radiocarbon age of the sample
#   }
# 
#   # Set Prior for Each Region
#   for (k in 1:n_areas){
#     a[k] ~ dunif(4000, 12000);
#     b[k] ~ dunif(4000, 12000);
#     constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
#   }
# })
# 
# 
# #Initialise parameters ---- 
# constants2 <- list(n_dates = f_constants$n_dates,
#                    n_areas = f_constants$n_areas,
#                    id_date_area = f_constants$id_date_area,
#                    calBP = intcal20$CalBP,
#                    C14BP = intcal20$C14Age,
#                    C14err = intcal20$C14Age.sigma)
# 
# 
# 
# inits2 <- list(a=init_a,
#                b=init_b,
#                theta=theta_init)
# 
# 
# #Run MCMC ----
# mcmc.samples2<- nimbleMCMC(code = model2,
#                            constants = constants2,
#                            data = dat,
#                            niter = 200, 
#                            nchains = 3, 
#                            thin=10, 
#                            nburnin = 100,
#                            monitors=c('a','b','theta'), 
#                            inits=inits2, 
#                            samplesAsCodaMCMC=TRUE)
# 
# #Diagnostics ----
# rhat2  <- gelman.diag(mcmc.samples2, multivariate = FALSE)
# ess2  <- effectiveSize(mcmc.samples2)

#===============================================================================
##MODEL 3 -- Model arrival time with additional variables

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
height_sf$area_id <- as.integer(st_within(height_sf$geometry, sq_grid$land_within_poly))

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

# #Plot square areas coloured by site frequency, fitness values or elevation
# ggplot(data = sq_grid_conc) +
#   geom_sf(data = as(sample_win_sp, 'sf')) +
#   geom_sf(alpha=0.5, aes(fill=site_freq)) + #Alternatively: fill=sq_fitness_value #fill=sq_height_value
#   geom_sf(data = as(f_sites, 'sf'), size=3, alpha=0.5, aes(colour="purple")) + #sites
#   geom_sf_label(aes(label = area_id), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
#   labs(x = "Longitude", y = "Latitude") +
#   scale_fill_viridis()

#-------------
## Compute Great-Arc Distances in km between area centers ---
#We take the center of each area k to be a point representing that whole area
sq_area_centers <- sq_grid$area_center %>% as('Spatial')
sq_dist_mat <- spDists(sq_area_centers, longlat=TRUE) #Inter-area distance matrix: each area's distance from every other area.

#The square area containing our putative origin
origin_sf <- as(f_sites, 'sf') %>% 
  filter(siteID == f_constants$origin_siteID)
origin_sf$area_id <- as.integer(st_within(origin_sf$geometry, sq_grid$geometry))  
sq_origin_area_center <- sq_area_centers[origin_sf$area_id] #The center of the square containing our putative origin

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
#Cumulative friction of square area k from origin square area
#Initially: cf = sq_dist_org - w_f*sq_fit_org + w_e*sq_elev_org
#The signs are based on the following initial assumptions:
#(i) Moving to lower elevation is good -- leads to a decrease in cumulative friction
#(ii) Moving to higher fitness values is good -- leads to a decrease in cumulative friction
#The weights, w_f and w_e, currently ensure that sq_fit_org and sq_eleve_org are of the same order of magnitude
#TODO: explore these assumptions and give weights to each of these parameters based on the relative importance

w_e <- 1 #To refine. Arrived at by looking at (max(sq_grid_conc$sq_elev_org, na.rm=TRUE)-min(sq_grid_conc$sq_elev_org, na.rm=TRUE)) compared to (max(sq_grid_conc$sq_dist_org, na.rm=TRUE)-min(sq_grid_conc$sq_dist_org, na.rm=TRUE))
w_f <- 600 #To refine. Arrived at by looking at (max(sq_grid_conc$sq_fit_org, na.rm=TRUE)-min(sq_grid_conc$sq_fit_org, na.rm=TRUE)) compared to (max(sq_grid_conc$sq_dist_org, na.rm=TRUE)-min(sq_grid_conc$sq_dist_org, na.rm=TRUE))
  
sq_grid_conc <- sq_grid_conc %>% 
  mutate(cf = sq_dist_org - w_f*sq_fit_org + w_e*sq_elev_org) %>% 
  mutate(cf = case_when(cf<0 ~ 50, cf>=0 ~ cf)) #If cumulative friction is negative -- assign a small positive value. It takes a minimum amount of effort to move. #TODO: adjust amount. Currently chosen at random.


#-------------------------------------------------------------------------------
## Store Output ----
save(sq_grid_conc, file=here('output', 'cumulative_friction.RData'))
     
     
#-------------------------------------------------------------------------------  
# model3 <- nimbleCode({
#   for (i in 1:n_dates)
#   {
#     theta[i] ~ dunif(max=a[id_date_area[i]], min=b[id_date_area[i]]); # Alternatively k=id_area[id_sites[i]] -- but, dynamic indexing not allowed
#     # Calibration
#     mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #C14 age on the relevant calibration curve
#     sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]); #error on the calibration curve
#     error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2); #the samples' C14 error + error on calibration curve
#     cra[i] ~ dnorm(mean=mu[i], sd=error[i]); #observed radiocarbon age of the sample
#   }
#   
#   # Set Prior for Each Region
#   for (k in 1:n_areas){
#     #sigma[k] ~ dcar_normal(adj[L], weights[L], num[k], tau1, zero_mean =0); #spatial random effect modelled using ICAR
#     cf[k] ~ dunif(0, 1000);
#     a[k] <- beta_0 + beta_1*cf[k] #beta_0 + (beta_1+sigma[k])*cf[k] #T(dnorm(beta_0, sigma), min = 4000, max=12000); #dnorm(beta_0 + beta_1*cf[k], sigma); #(beta_1 + sigma[k])*c[k];
#     b[k] ~ dunif(4000, 12000);
#     constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
#   }
#   
#   #priors
#   beta0 ~ dnorm(8000, sd=200);
#   beta1 ~ dexp(1);
#   tau1 <- 1/sigma1^2;
#   sigma ~ dunif(0,100); #dexp(0.01)
#   #sigma1 ~ dunif(0,100)
#   
# })
# 
# 
# # x <- seq(2000, 14000, length=10000)
# # y <- dunif(x, 7000, 9000)
# # z <- dnorm(x, mean=y, sd=200)
# # plot(x, z, type="l", lwd=1)
# 
# 
# 
# # MCMC initialisation ----
# set.seed(1357)
# inits3  <-  list()
# inits3$theta  <- theta_init
# inits3$beta0 <- rnorm(1, 8000, 200)
# inits3$beta1 <- rexp(1, rate=2)
# inits3$a <- init_a
# inits3$b <- init_b
# inits3$cf <- sq_grid_conc$cf
# #inits3$sigma1 <- runif(1,0,100)
# inits3$sigma  <- rexp(1, 0.01)
# #inits3$sigma <- rnorm(f_constants$n_areas, sd=0.001)
#   
# #Initialise parameters ---- 
# constants3 <- list(n_dates = f_constants$n_dates,
#                    n_areas = f_constants$n_areas,
#                    id_date_area = f_constants$id_date_area,
#                    calBP = intcal20$CalBP,
#                    C14BP = intcal20$C14Age,
#                    C14err = intcal20$C14Age.sigma)
# 
# 
# 
# #Run MCMC ----
# mcmc.samples3<- nimbleMCMC(code = model3,
#                            constants = constants3,
#                            data = dat,
#                            niter = 200, 
#                            nchains = 3, 
#                            thin=30, 
#                            nburnin = 100,
#                            monitors=c('a','b','theta'), 
#                            inits=inits3, 
#                            samplesAsCodaMCMC=TRUE)
# 
# #Diagnostics ----
# rhat3  <- gelman.diag(mcmc.samples3, multivariate = FALSE)
# ess3  <- effectiveSize(mcmc.samples3)
# 
# 
