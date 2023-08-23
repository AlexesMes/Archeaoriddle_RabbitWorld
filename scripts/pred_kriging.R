#-------------------------------------------------------------------------------
#EXPLORATORY DATA ANALYSIS AND BASIC KRIGING MODEL
#See: https://github.com/guzmanlopez/IntroductionToGeostatistics

#===============================================================================
# Load Libraries and Data ----
library(sp)
library(sf)
library(gstat)
library(dplyr)
library(ggplot2)
library(here)
library(maptools)
library(tidyr)
library(stringr)
library(raster)
library(terra)
library(stars)


library('geoR')
library('moments')
library('reshape2')
library('forecast')
library('ncf')
library('spdep')
library('MASS')
library('relaimpo')
library('gridExtra')
library('corrgram') 
library('RColorBrewer')


# Load and prepare data ----
load(here('data','c14data_farmers.RData'))

#Read in elevation and resource data ----
fitness_sf <- read_stars("maps/resources.tiff") %>% 
  st_as_sf() %>% 
  rename("fitness_value" = "resources.tiff")

height_sf <- read_stars("maps/east_narnia4x.tif") %>% 
  st_as_sf() %>% 
  rename("height_value" = "east_narnia4x.tif")


#-------------
# glimpse(f_dateInfo)
# glimpse(f_sites)

##Data concentration per site
# f_sites %>% as.data.frame %>% 
#   ggplot(aes(lon, lat)) + geom_point(aes(size=n_dates), color="blue", alpha=3/4) + 
#   ggtitle("Date Concentration per Site") + coord_equal() + theme_bw()

#--------------
#Grid concentrations
f_area_freq  <- plyr::count(f_siteInfo, 'area_id')

sq_site_conc <- sq_grid %>% 
  left_join(f_area_freq, by=join_by(area_id)) %>% 
  mutate(freq = case_when((contains_land==TRUE & is.na(freq)) ~ 0, is.na(contains_land) ~ NA, !is.na(freq) ~ freq)) #0 if there are no sites, NA if the square area is in the sea

#--------------  
#Aggregate fitness and height values for each square area ----
#Fitness
fitness_sf$area_id <- as.integer(st_within(fitness_sf$geometry, sq_grid$geometry))

fitness_sf <- fitness_sf %>% 
  group_by(area_id) %>% 
  summarize(sq_fitness_value = mean(fitness_value))

fitness_df <- as.data.frame(fitness_sf) %>% dplyr::select(area_id, sq_fitness_value)


#Height
height_sf$area_id <- as.integer(st_within(height_sf$geometry, sq_grid$geometry))

height_sf <- height_sf %>% 
  group_by(area_id) %>% 
  summarize(sq_height_value = mean(height_value))

height_df <- as.data.frame(height_sf) %>% dplyr::select(area_id, sq_height_value)

#-----
#Join grid-level data
sq_site_conc <- sq_site_conc %>% 
  left_join(fitness_df, by=join_by(area_id)) %>% 
  left_join(height_df, by=join_by(area_id)) %>% 
  mutate(lat=st_coordinates(area_center)[,1], 
         long=st_coordinates(area_center)[,2])

# ggplot(data = sq_site_conc) +
#   geom_sf(alpha=0.5, aes(fill=freq)) + #sq grid
#   geom_sf(data = as(f_sites, 'sf'), size=3, alpha=0.5, aes(colour="purple")) + #sites
#   geom_sf_label(aes(label = area_ID), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
#   labs(fill="Square grid", x = "Longitude", y = "Latitude") +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")


sq_site_conc_sp <- sq_site_conc %>% as('Spatial')

#===============================================================================
##Frequency distributions----
ggSqData <- ggplot(data = sq_site_conc_sp@data)

#Site concentration
ggHist_siteconc <- ggSqData +
  geom_histogram(aes(x = sq_site_conc_sp$freq),
                 binwidth = 0.5, fill = "orange", col = "white", alpha = 0.8) +
  labs(title = "", x = "Number of sites (per square)", y = "Frequency")

ggBoxPlot_siteconc <- ggSqData +
  geom_boxplot(aes(y = sq_site_conc_sp$freq, 
                   x = rep("RM", length(sq_site_conc_sp$freq))), 
               fill = "orange", col = "black", alpha = 0.8) + 
  labs(title = "Site Concentration", x = "", y = "") + 
  coord_flip()

#Fitness
ggHist_fitness <- ggSqData +
  geom_histogram(aes(x = sq_site_conc_sp$sq_fitness_value), 
                 binwidth = 0.05, fill = "darkblue", col = "white", alpha = 0.7) +
  labs(title = "", x = "Fitness value (per square)", y = "Frequency") 

ggBoxPlot_fitness <- ggSqData +
  geom_boxplot(aes(y = sq_site_conc_sp$sq_fitness_value, 
                   x = rep("PN", length(sq_site_conc_sp$sq_fitness_value))), 
               fill = "darkblue", col = "black", alpha = 0.7) + 
  labs(title = "Fitness", x = "", y = "") + 
  coord_flip()

#Elevation
ggHist_elev <- ggSqData +
  geom_histogram(aes(x = sq_site_conc_sp$sq_height_value), 
                 binwidth = 10, fill = "darkred", col = "white", alpha = 0.7) +
  labs(title = "", x = "Height in meters (per square)", y = "Frequency") 

ggBoxPlot_elev <- ggSqData +
  geom_boxplot(aes(y = sq_site_conc_sp$sq_height_value, 
                   x = rep("PN", length(sq_site_conc_sp$sq_height_value))), 
               fill = "darkred", col = "black", alpha = 0.7) + 
  labs(title = "Elevation", x = "", y = "") + 
  coord_flip()

#Plot
grid.arrange(ggHist_siteconc, 
             ggBoxPlot_siteconc, 
             ggHist_fitness, 
             ggBoxPlot_fitness, 
             ggHist_elev, 
             ggBoxPlot_elev,
             ncol = 2, nrow = 3)
#-------------------------------------------------------------------------------
##Central trend and dispersions----
SummaryStatistics <- function(df) {
  
  ncol = ncol(df)
  
  #Deeclar objects
  max = numeric()
  min = numeric()
  
  mean = numeric()
  median = numeric()
  mode = vector()
  
  variance = numeric()
  sd = numeric()
  cv = numeric()
  skewness = numeric()
  kurtosis = numeric()
  
  NAs = numeric()
  clases = character()
  
  
  #Function to calculate the mode
  mode <- function(x) {
    ux = unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  for(i in 1:ncol) {
    
    if(is.numeric(df[,i])) {
      median[i] = median(df[,i], na.rm = TRUE)
      mean[i] = mean(df[,i], na.rm = TRUE)
      #mode[i] = mode(x = df[,i])
      variance[i] = var(df[,i], na.rm = TRUE)
      sd[i] = sd(df[,i], na.rm = TRUE)
      cv[i] = 100 * (sd(df[,i], na.rm = TRUE) / mean(df[,i], na.rm = TRUE))
      skewness[i] = moments::skewness(df[,i], na.rm = TRUE)
      kurtosis[i] = moments::kurtosis(df[,i], na.rm = TRUE)
      min[i] = min(df[,i], na.rm = TRUE)
      max[i] = max(df[,i], na.rm = TRUE)
      NAs[i] = length(which(is.na(df[,i]) == TRUE))
      clases[i] = class(df[,i])
    } 
    
    else {
      median[i] = NA
      mean[i] = NA
      #mode[i] = mode(x = df[,i])
      variance[i] = NA
      sd[i] = NA
      cv[i] = 100 * NA
      skewness[i] = NA
      kurtosis[i] = NA
      min[i] = NA
      max[i] = NA
      NAs[i] = length(which(is.na(df[,i]) == TRUE))
      clases[i] = class(df[,i])
    }
  } 
  
  summary = data.frame("minimum" = min,
                       "maximum" = max,
                       "mean" = mean,
                       "median" = median,
                       #"mode" = mode,
                       "variance" = variance,
                       "sd" = sd,
                       "cv" = cv,
                       "skewness" = skewness,
                       "kurtosis" = kurtosis,
                       "NAs" = NAs,
                       "class" = clases)
  
  rownames(summary) <- colnames(df)
  
  return(summary)
  
}

sum_variables_df <- as.data.frame(sq_site_conc) %>% 
  dplyr::select(area_id, freq, sq_fitness_value, sq_height_value)

sq_summary_stats <- SummaryStatistics(df = sum_variables_df)

#-------------------------------------------------------------------------------
##Test for normalcy----
shapiro.test(sq_site_conc$sq_height_value) 
shapiro.test(sq_site_conc$sq_fitness_value) 
shapiro.test(sq_site_conc$freq) 

#--------------------
##Transformations----
#attempt to transform the measured values to a new scale on which the distribution is more nearly normal to overcome difficulties arising from departures from normality.

#Site Concentration--
#Sqrt
ConcSqrt <- sqrt(sq_site_conc$freq)
shapiro.test(ConcSqrt)

#Fitness--
#Log
FitLog <- log(sq_site_conc$sq_fitness_value)
shapiro.test(FitLog)

#Sqrt
FitSqrt <- sqrt(sq_site_conc$sq_fitness_value)
shapiro.test(FitSqrt)

# Box-Cox
lambdaFit <- BoxCox.lambda(x = sq_site_conc$sq_fitness_value)
FitBC <- BoxCox(x = sq_site_conc$sq_fitness_value, lambda = lambdaFit) 
shapiro.test(FitBC) 


#Elevation--
#Log
ElevLog <- log(sq_site_conc$sq_height_value)
shapiro.test(ElevLog)

#Sqrt
ElevSqrt <- sqrt(sq_site_conc$sq_height_value)
shapiro.test(ElevSqrt)


#None of the data transformations allowed failure to reject the null hypothesis 
#of normality of the test of Shapiro-Wilk -- i.e. allow us to assume data generated from underlying normal distribution#
#Some transformations offer improvements...

#-------------
#Q-Q Plots----
#Perform a "visual" analysis of normality using Q-Q plots

#Site Concentration--
ggQQ_freq <- ggSqData + 
  geom_qq(aes(sample = freq)) +
  geom_abline(intercept = mean(sq_site_conc_sp@data$freq, na.rm = TRUE), 
              slope = sd(sq_site_conc_sp@data$freq, na.rm = TRUE), col = "red") +
  labs(title = "Site Concentration QQ Plot", x = "Normal value", 
       y = "Data value")

ggQQ_freq_sqrt <- ggSqData + 
  geom_qq(aes(sample = sqrt(freq))) +
  geom_abline(intercept = mean(sqrt(sq_site_conc_sp@data$freq), na.rm = TRUE), 
              slope = sd(sqrt(sq_site_conc_sp@data$freq), na.rm = TRUE), col = "red") +
  labs(title = "Site Concentration QQ Plot", x = "Normal value", 
       y = "Sqrt Transformed Data value")

#Fitness--
ggQQ_fitness <- ggSqData + 
  geom_qq(aes(sample = sq_fitness_value)) +
  geom_abline(intercept = mean(sq_site_conc_sp@data$sq_fitness_value, na.rm = TRUE), 
              slope = sd(sq_site_conc_sp@data$sq_fitness_value, na.rm = TRUE), col = "red") +
  labs(title = "Fitness QQ Plot", x = "Normal value", 
       y = "Data value")

ggQQ_fitness_log <- ggSqData + 
  geom_qq(aes(sample = log(sq_fitness_value))) +
  geom_abline(intercept = mean(log(sq_site_conc_sp@data$sq_fitness_value), na.rm = TRUE), 
              slope = sd(log(sq_site_conc_sp@data$sq_fitness_value), na.rm = TRUE), col = "red") +
  labs(title = "Fitness QQ Plot", x = "Normal value", 
       y = "Log Transformed Data value")

ggQQ_fitness_sqrt <- ggSqData + 
  geom_qq(aes(sample = sqrt(sq_fitness_value))) +
  geom_abline(intercept = mean(sqrt(sq_site_conc_sp@data$sq_fitness_value), na.rm = TRUE), 
              slope = sd(sqrt(sq_site_conc_sp@data$sq_fitness_value), na.rm = TRUE), col = "red") +
  labs(title = "Fitness QQ Plot", x = "Normal value", 
       y = "Sqrt Transformed Data value")

ggQQ_fitness_BC <- ggSqData + 
  geom_qq(aes(sample = FitBC)) +
  geom_abline(intercept = mean(FitBC, na.rm = TRUE), 
              slope = sd(FitBC, na.rm = TRUE), col = "red") +
  labs(title = "Fitness QQ Plot", x = "Normal value", 
       y = "Box-Cox Transformed Data value")


#Elevation--
ggQQ_height <- ggSqData + 
  geom_qq(aes(sample = sq_height_value)) +
  geom_abline(intercept = mean(sq_site_conc_sp@data$sq_height_value, na.rm = TRUE), 
              slope = sd(sq_site_conc_sp@data$sq_height_value, na.rm = TRUE), col = "red") +
  labs(title = "Elevation QQ Plot", x = "Normal value", 
       y = "Data value")

ggQQ_height_log <- ggSqData + 
  geom_qq(aes(sample = log(sq_height_value))) +
  geom_abline(intercept = mean(log(sq_site_conc_sp@data$sq_height_value), na.rm = TRUE), 
              slope = sd(log(sq_site_conc_sp@data$sq_height_value), na.rm = TRUE), col = "red") +
  labs(title = "Elevation QQ Plot", x = "Normal value", 
       y = "Log Transformed Data value")

ggQQ_height_sqrt <- ggSqData + 
  geom_qq(aes(sample = sqrt(sq_height_value))) +
  geom_abline(intercept = mean(sqrt(sq_site_conc_sp@data$sq_height_value), na.rm = TRUE), 
              slope = sd(sqrt(sq_site_conc_sp@data$sq_height_value), na.rm = TRUE), col = "red") +
  labs(title = "Elevation QQ Plot", x = "Normal value", 
       y = "Sqrt Tranformed Data value")
#Sqrt transformation of height data comes the closest to normally distributed data...

#Plot
grid.arrange(ggQQ_freq, ggQQ_freq_sqrt, 
             ggQQ_fitness, ggQQ_fitness_log, ggQQ_fitness_sqrt, ggQQ_fitness_BC, 
             ggQQ_height, ggQQ_height_log, ggQQ_height_sqrt, 
             ncol = 4, nrow = 3,
             layout_matrix = rbind(c(1, 2, NA, NA), c(3, 4, 5, 6), c(7, 8, 9, NA)))

#-------------------------------------------------------------------------------
##Trend Analysis ----
#Discover trends between the response variable and co-variables using scatter plots with smoothed lines and/or correlations and/or correlograms.

# Site Concentration vs Longitude
ggTend_ConcLong <- ggSqData +
  geom_point(aes(y = freq, x = long)) + 
  geom_smooth(aes(y = freq, x = long)) + 
  labs(title = "", x = "Longitude", y = "Site Concentration (per square)")

# Site Concentration vs Latitude
ggTend_ConcLat <- ggSqData +
  geom_point(aes(y = freq, x = lat)) + 
  geom_smooth(aes(y = freq, x = lat)) + 
  labs(title = "", x = "Latitude", y = "Site Concentration (per square)")

# Site Concentration vs Elevation
ggTend_ConcElev <- ggSqData +
  geom_point(aes(y = freq, x = sq_height_value)) + 
  geom_smooth(aes(y = freq, x = sq_height_value)) + 
  labs(title = "", x = "Elevation (m)", y = "Site Concentration (per square)")

# Site Concentration vs Fitness Values
ggTend_ConcFit <- ggSqData +
  geom_point(aes(y = freq, x = sq_fitness_value)) + 
  geom_smooth(aes(y = freq, x = sq_fitness_value)) + 
  labs(title = "", x = "Fitness Values", y = "Site Concentration (per square)")


# Elevation vs Fitness Values
ggTend_FitElev <- ggSqData +
  geom_point(aes(y = sq_height_value, x = sq_fitness_value)) + 
  geom_smooth(aes(y = sq_height_value, x = sq_fitness_value)) + 
  labs(title = "", x = "Fitness Values", y = "Elevation (m)") 
#High fitness values on the coast (elevation lowest in these square areas, since skewed by negative elevation of the ocean parts)

#Plot
grid.arrange(ggTend_ConcLong, 
             ggTend_ConcLat, 
             ggTend_ConcElev, 
             ggTend_ConcFit, 
             ggTend_FitElev)

#-------------------------------------------------------------------------------
##Correlation----
#Pearson's correlation coefficient was calculated to quantify the relationship between the variables for which trend was observed and for which not. To complement the interpretation, correlograms were made to visualize the correlation matrix
#Not particularly helpful here since we don't have variables which are roughly linearly related to one another...

# Site Concentration vs Longitude
cor.test(x = sq_site_conc_sp$long, 
         y = sq_site_conc_sp$freq, 
         method = "pearson")

# Site Concentration vs Latitude
cor.test(x = sq_site_conc_sp$lat, 
         y = sq_site_conc_sp$freq, 
         method = "pearson")

# Site Concentration vs Elevation
cor.test(x = sq_site_conc_sp$sq_height_value, 
         y = sq_site_conc_sp$freq, 
         method = "pearson")

# Site Concentration vs Fitness
cor.test(x = sq_site_conc_sp$sq_fitness_value, 
         y = sq_site_conc_sp$freq, 
         method = "pearson")

# Elevation vs Fitness
cor.test(x = sq_site_conc_sp$sq_fitness_value, 
         y = sq_site_conc_sp$sq_height_value, 
         method = "pearson")


##Correlogram--
#The blue colored cells with a line fill pattern with a positive slope represent positive correlations and the red colored cells with a fill pattern lines with a negative slope represent negative correlations. Color intensity is proportional to the magnitude of the Pearson's correlation coefficient.
corrgram(data.frame("Site_Concentration" = sq_site_conc_sp$freq,
                    "Long" = sq_site_conc_sp$long, 
                    "Lat" = sq_site_conc_sp$lat, 
                    "Elev" = sq_site_conc_sp$sq_height_value,
                    "Fitness" = sq_site_conc_sp$sq_fitness_value),
         cor.method="pearson",
         order = NULL,
         lower.panel = panel.shade,
         upper.panel = NULL, 
         text.panel = panel.txt,
         main = "Correlogram")

#-------------------------------
##Akaike information criterion--
#Evaluate trend models

#Models
freq_M1 <- lm(formula = freq ~ sq_fitness_value + sq_height_value, data = sq_site_conc_sp, na.action=na.exclude)
freq_M2 <- lm(formula = freq ~ sq_fitness_value, data = sq_site_conc_sp, na.action=na.exclude)
freq_M3 <- lm(formula = freq ~ sq_height_value, data = sq_site_conc_sp, na.action=na.exclude)

#AIC calculations
AIC(freq_M1, freq_M2, freq_M3)
#The preferred model is the one with the minimum AIC value
#WAIC is the more preferred metric these days 


##Calculate relative importance of variables--
impRel_M1 <- calc.relimp(freq_M1, 
                         type = c("lmg", "last", "first", "pratt"),
                         rela = TRUE)
#The variance explained by the model is 0.4% for the site concentration. The fitness value is more important than the elevation, but this doesn't really matter with such a poor model.

#-------------------------------------------------------------------------------
##Spatial Autocorrelation----

#Which points values are atypically different from their neighborhood point values, you can use Moran Index to quantify it and make a "posting" plot to explore your data.

sq_site_conc$sq_freq_M1 <- 
  freq_M1$coefficients["sq_fitness_value"] * sq_site_conc$sq_fitness_value +
  freq_M1$coefficients["sq_height_value"] * sq_site_conc$sq_height_value +
  freq_M1$coefficients["(Intercept)"]

sq_site_conc$resfreq_M1 <- residuals(freq_M1)

moran_siteConc <- correlog(x = sq_site_conc$long, 
                            y = sq_site_conc$lat, 
                            z = sq_site_conc$resfreq_M1, 
                            increment = 100, 
                            resamp = 0, 
                            latlon = FALSE, 
                            na.rm = TRUE)

#TODO: finish this section... plot...
#-------------------------------------------------------------------------------
##Empiric variogram---- 
#Characterizing the spatial processes of your data
#geoR Package
dfM1 <- data.frame("Lat" = sq_site_conc$lat,
                   "Long" = sq_site_conc$long, 
                   "sq_height_value" = sq_site_conc$sq_height_value,
                   "sq_fitness_value" = sq_site_conc$sq_fitness_value,
                   "freq" = sq_site_conc$freq) %>% na.omit()

#Object geodata
geodataM1 <- as.geodata(obj = dfM1, 
                        coords.col = 1:2, 
                        data.col = 5, 
                        covar.col = c(3,4))

#Variogram
#TODO: This section isn't working....
varGeoM1 <- variog(geodataM1, 
                   trend = trend.spatial(data~sq_height_value+sq_fitness_value, geodataM1), 
                   max.dist = 1000, 
                   breaks = seq(from = 0, to = 500, by = 20))

varGeoC1dir0 <- variog(geodataM1, trend = trend.spatial(data~sq_height_value+sq_fitness_value, geodataM1), 
                       direction = 0, tolerance = 45/2, unit.angle = "degrees",
                       max.dist = 400, breaks = seq(from = 0, to = 500, by = 20))


#Plot
plot(varGeoM1, pch = 19, col = "#0080FF", ylab = "Semi-variance", 
     xlab = "Variogram of site concentration", main = "M1")

#-----
#gstat Package (alternative method)
varo_M1 <- variogram(freq ~ sq_fitness_value + sq_height_value,
                     data = na.omit(sq_site_conc),
                     width = 20, 
                     cutoff = 1000)

# Plotear variograma
plot(varo_M1, xlab = "Distance", ylab = "Semi-variance",
     main = "Variogram of site concentration", pch = 19)

##Ordinary Kriging Model----
sq_with_info <- sq_site_conc %>% filter(freq!=0)
sq_without_info <- sq_site_conc %>% filter(freq==0)

gstat::krige(formula = freq ~ 1, 
             locations = sq_with_info, 
             newdata = sq_without_info$area_center, 
             model = varo_M1)





















