# Load Libraries and Data ----
library(maptools)
#library(raster)
library(stars)
library(sf)
library(stringr)
library(dplyr)
library(here)
library(ggplot2)
library(ggthemes)
library(viridis)
library(parallel)
library(cowplot)

#-------------------------------------------------------------------------------
## Load and prepare data ----
load(here('data','c14data.RData'))
load(here('data','c14data_farmers.RData'))


#Read in basemaps ----
height_map <- read_stars("maps/east_narnia4x.tif")
fitness_map <- read_stars("maps/resources.tiff")


#-------------------------------------------------------------------------------
## Data preparation ----

#Convert to sf objects
sites_sf <- sf::st_as_sf(siteInfo, 
                         coords = c("lon", "lat"), 
                         remove = F, 
                         crs = 4326, 
                         na.fail = F)


#-------------------------------------------------------------------------------
## Plot Data  ----

#Elevation Map
plt <- ggplot() +
  geom_stars(data=height_map) +
  geom_sf(data = sites_sf,
          aes(colour=economy),
          size = 2,
          alpha=0.8) +
  geom_vline(xintercept = seq(-4, 1, 0.5), colour = 'lightgrey') + #Make grid divisions explicit 
  geom_hline(yintercept = seq(-1, 4, 0.5), colour = 'lightgrey') +
  labs(fill="Elevation", x = "Longitude", y = "Latitude") +
  scale_fill_viridis()
  
#Environmental Fitness Map
ef_plt <- ggplot() +
  geom_stars(data=fitness_map) +
  geom_sf(data = sites_sf,
          aes(colour=economy),
          size = 2,
          alpha=0.8) +
  geom_vline(xintercept = seq(-4, 1, 0.5), colour = 'lightgrey') + #Make grid divisions explicit 
  geom_hline(yintercept = seq(-1, 4, 0.5), colour = 'lightgrey') +
  labs(fill="Environmental fitness probability", x = "Longitude", y = "Latitude") +
  scale_fill_viridis()


#-------------------------------------------------------------------------------
## Save Data  ----

#Plot Elevation Map
pdf(file=here('output', 'figures', 'elevation_map.pdf'), width=8.5, height=7)
cowplot::ggdraw() +
  draw_plot(plt)
dev.off()

#Plot Environmental Fitness Map
pdf(file=here('output', 'figures', 'fitness_map.pdf'), width=8.5, height=7)
cowplot::ggdraw() +
  draw_plot(ef_plt)
dev.off()

#===============================================================================
### Plot GRID areas  ---- FIGURE map_fig_grid

pdf(file=here('output', 'figures', 'map_fig_grid.pdf'), width=8.5, height=7)

ggplot(data = sq_grid) +
  geom_sf(data = as(sample_win_sp, "sf")) + #sampling window with coastal buffer
  geom_sf(alpha=0.1) + #sq grid
  geom_sf(data = as(f_sites, 'sf'), size=3, alpha=0.5, aes(colour="purple")) + #sites
  geom_point(aes(x=f_constants$origin_point[1], y=f_constants$origin_point[2]), colour="darkgreen", size=3) +
  geom_sf_label(aes(label = area_id), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
  #geom_sf(data = sq_grid$area_center, size=2, alpha=1, aes(color = "purple")) + #sq-origins
  labs(fill="Square grid", x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

dev.off()





