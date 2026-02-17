################################################################################
# Script to estimate light availability ########################################
#
# Author: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 17/02/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
# Cite as:
# García-Gómez, G., Sánchez-Hernandez, J., Mas Gutiérrez, J.A., & Arranz, I. (2026). 
# Community size structure of microbial plankton is associated with functional composition along a resource gradient in mountain lakes.
# OSF (DOI 10.17605/OSF.IO/2493Q)
#
# ------------------------------------------------------------------------------
rm(list=ls())# clear the work environment
today <- format(Sys.Date(),"%Y%m%d")# setting the date
# ------------------------------------------------------------------------------
# It needs to be set to Project directory
getwd()# to check
#
# Load libraries ####
#
# Install CRAN packages if needed
#
# required packages:
cran_pkg <- c("dplyr",
              "ggplot2", 
              "terra",
              "suncalc",
              "nasapower",
              "lubridate")

# missing packages:
miss_cran <- cran_pkg[!cran_pkg %in% rownames(installed.packages())]


# install missing packages:
for(pkg in miss_cran){
  install.packages(pkg, dependencies = TRUE)
}

# please check package versions:
# (see example below, but checking all packages is recommended)
#
# library(remotes) # needed to install a specific package version
packageVersion("ggplot2") # should be ‘3.5.2’; otherwise run the next line:
# remotes::install_version("ggplot2", version = "3.5.2", dependencies = TRUE)

packageVersion("terra") # should be ‘1.8.42’; otherwise run the next line:
# remotes::install_version("terra", version = "1.8.42", dependencies = TRUE)

packageVersion("suncalc") # should be ‘0.5.1’; otherwise run the next line:
# remotes::install_version("suncalc", version = "0.5.1", dependencies = TRUE)

packageVersion("nasapower") # should be ‘4.2.5’; otherwise run the next line:
# remotes::install_version("nasapower", version = "4.2.5", dependencies = TRUE)

# NOTE that your R session may have older versions of package dependencies
# that are not updated automatically even if you install the right version
# of required packaged

# load libraries:
library(dplyr)
library(ggplot2)

# handling topography and satellite data:
library(terra)
library(suncalc) # sun position
library(nasapower) # daily irradiance

# handling dates and times:
library(lubridate)

# Content ####
#
## [1] Load processed datasets
## [2] Statistical analysis
## [3] Save merged dataset and results

# [1] Load processed datasets ####
#
## [1.1] Access and load processed data ####

# Load satellite data:
#
# Copernicus GLO-30 Digital Elevation Model (Geotiff)
#
# Obtained from OpenTopography website 
# (available on a free basis for the general public, but creating an account is necessary).
# Check link: https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.3
#
# Overview of file (from OpenTopography website):
# The Copernicus DEM is a Digital Surface Model (DSM) which represents the surface of the Earth including buildings, 
# infrastructure and vegetation. This DSM is derived from an edited DSM named WorldDEM, where flattening of water bodies 
# and consistent flow of rivers has been included. In addition, editing of shore- and coastlines, special features such as 
# airports, and implausible terrain structures has also been applied.
#
# Resolution: 30 m
# Survey Date: 01/01/2011-07/01/2015
# Funder: European Space Agency (ESA)
#
# Citation:
# European Space Agency (2024). 
# Copernicus Global Digital Elevation Model. 
# Distributed by OpenTopography. https://doi.org/10.5069/G9028PQB. Accessed [date of access YYYY-MM-DD]
#
# Note: boundaries in this map were selected to include only
# the coordinates of the lakes investigated in our study.
# Other boundaries and map dimenssions can be selected.
#
# See terms and conditions of GLO-30 DEM here:
# https://docs.sentinel-hub.com/api/latest/static/files/data/dem/resources/license/License-COPDEM-30.pdf
#
# Load .tif file:
dem <- rast("../raw_data/appRasterSelectAPIService1764233714639-1101222210.tif")
#
# Check metadata:
dem # Digital elevation model
crs(dem) # coordinate system
plot(dem) # visualise map (see the two Iberian massifs studied here: 
# Gredos (southern range of map) and Segundera (northern range))

# Project DEM to UTM (meters)
dem_m <- project(dem, "EPSG:25830")  # UTM zone 30N, Spain / in meters

## [1.2] Load lake locations (geographical coordinates and altitude) ####
lakes_location <- read.csv("../raw_data/lakes_location.csv")

lake_df <- 
  
  lakes_location %>%
  
  # get rid of source column to ease data overview
  dplyr::select(-source) %>%
  
  # let's include sampling dates of each lake so we can estimate
  # light intensity for the days prior to sampling collection:
  mutate(date_samp = case_when(
    
    lake %in% c("Trampal 1", "Trampal 2", "Trampal 3") ~ as.Date("2024-09-03"),
    lake == "Negra" ~ as.Date("2024-09-09"),
    lake == "Cuadrada" ~ as.Date("2024-09-10"),
    lake == "Grande" ~ as.Date("2024-10-19"),
    lake == "Peces" ~ as.Date("2024-09-15"),
    lake == "Yeguas" ~ as.Date("2024-09-16"),
    lake == "Payon" ~ as.Date("2024-09-18"),
    lake == "Mancas" ~ as.Date("2024-09-18"),
    TRUE ~ NA)) %>%
  
  # add an id column to trace lake's identity
  mutate(site_id = seq(1, nrow(.), by = 1))

str(lake_df)
# OK

# Plot lake locations on the map (DEM) now

# First, set CRS:
lakes_vect <- vect(
  lake_df,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

# Project in the same CRS as DEM:
lakes_vect_proj <- project(lakes_vect, crs(dem))

# Visualise spatial location of lakes on DEM
plot(dem, col = terrain.colors(50))
points(lakes_vect_proj, pch = 16, cex = 0.75)
# all seems correct!

# [2] Estimate light intensity ####

## [2.1] Function to estimate horizon angles in a lake ####

# First, we define a function to computs horizon angle for one azimuth for each lake

# Note that we check different values of distance (maxdist parameter)
# lake by lake to estimate shading by surronding topology

horizon_angle <- function(dem_m, # DEM
                          x, y, # lake coordinates 
                          azimuth_deg, # azimuth degrees
                          maxdist, # max. distance that surrounding topology may project shades onto lakes
                          step = 30 # number of points sampled along azimuth line
) {
  
  # Convert azimuth to radians:
  az <- azimuth_deg * pi/180
  
  # Sample points along the azimuth line
  distances <- seq(30, maxdist, by = step)
  xs <- x + distances * sin(az)
  ys <- y + distances * cos(az)
  
  pts <- cbind(xs, ys) # points along azimuth line
  
  # Extract elevation of surroundings:
  z <- terra::extract(dem_m, pts)[,1]
  
  if (all(is.na(z))) return(0) # if surrounding is flat
  
  # Height of lake point location:
  lake_z <- terra::extract(dem_m, cbind(x, y))[,1]
  
  # Angle to terrain
  angles <- atan((z - lake_z) / distances)
  
  # Maximum obstruction angle (shading) in degrees
  return(max(angles, na.rm = TRUE) * 180/pi)
}

## [2.2] Conversion factor: values ####
#
# Photosynthetically active radiation (PAR)
# fraction of daily irradiance that is PAR: photosynthetically active waveband 400-700 nm
PAR_fraction <- 0.48  # typically between 0.45-0.50 [see ref. 1 below] 
umol_per_J   <- 4.6   # µmol photons per J (depending on light spectrum, 4.6 is at green light, at 550 nm)[see ref. 2 and 3 below]

# to convert from MJ/day to J/seg:
MJ_to_J = 1e6
day_to_seg <- 24*60*60

conv_factor <- (MJ_to_J/day_to_seg) * PAR_fraction * umol_per_J
# conv. factor ~ MJ / m2 / day x 25.56 = µmol photons / m2 / seg

## [2.3] Calculate PAR across lakes ####

# Let's loop this function over the lake locations dataset:

# Create an empty data frame to store results:
lake_light_int.df <- data.frame()

# Loop over lakes now:
for(i in seq_along(unique(lake_df$lake))) {
  
  # select a lake:
  lake_coord_i <- lake_df[i,]
  
  message("Processing site: ", lake_coord_i$lake, " ~ ",
          
          lake_coord_i$site_id, "/", nrow(lake_df))
  
  # Set coordinate system of the lake
  lake_pt_i <- vect(data.frame(lon = lake_coord_i$longitude, 
                               lat = lake_coord_i$latitude), 
                    
                    crs="EPSG:4326")
  
  lake_pt.crs_i <- project(lake_pt_i, crs(dem_m)) # project to UTM meters (as in DEM)
  
  xy_i <- geom(lake_pt.crs_i)[, c("x","y")] # lake coordinates in UTM
  
  alt_est <- terra::extract(dem_m, cbind(xy_i[1], xy_i[2]))[[1]] # lake altitude (m a.s.l.)
  
  # Compute horizon for 360 azimuths:
  azimuths <- 0:359
  
  # Estimate shading from surrounding topology using different distances around the lake:
  maxdist <- c(100, 500, 1000, 3000, 5000, 8000)
  
  # Let's loop again within each lake (several values for a lake following maxdist values):
  
  # Create empty dataset to store results from this intra-loop:
  lake_d_df <- data.frame()
  
  for(d in seq_along(maxdist)){
    
    maxdist_d <- maxdist[[d]] # pick a distance value (in m)
    
    # calculate horizon angles for every maxdist value:
    horiz_i <- sapply(azimuths, function(a) horizon_angle(dem_m, xy_i[1], xy_i[2], a, maxdist_d))
    
    # Store it into a dataframe:
    horizon_df_i <- data.frame(
      azimuth = azimuths,
      horizon_angle = horiz_i
    )
    
    # Compute sun elevation & azimuth through time:
    date_t.samp <- lake_coord_i$date_samp # sampling date
    date_t.init <- date_t.samp - 30 # select 30 days before the sampling date
    
    # Obtain irradiance throughout the day:
    times_t.samp <- paste(date_t.init, "00:00:00")
    times_t.init <- paste(date_t.samp, "23:59:59")
    
    # Create 10-min intervals:
    times <- seq(
      as.POSIXct(times_t.samp, tz = "UTC"),
      as.POSIXct(times_t.init, tz = "UTC"),
      by = "10 min"
    )
    
    # Get sun position throughout the selected period:
    sun <- getSunlightPosition(
      date = times,
      lat = lake_coord_i$latitude,
      lon = lake_coord_i$longitude
    )
    
    # Create dataset including times, azimuths and corresponding sun elevation:
    sun_df <- data.frame(
      time = times,
      azimuth = (sun$azimuth * 180/pi + 180) %% 360, # convert from radians (suncalc default) to decimal angles
      elevation = sun$altitude * 180/pi
    )
    
    # Compare sun elevation to terrain (surrounding) horizon:
    sun_df.c <- sun_df %>%
      rowwise() %>%
      mutate(
        horizon_block = horizon_df_i$horizon_angle[which.min(abs(horizon_df_i$azimuth - azimuth))], # potential horizon block
        sun_above_terrain = elevation > horizon_block # is sunlight elevation above terrain elevation surrounding the lake?
      ) %>%
      
      data.frame()
    
    # Obtain daily irradiance throughout this period:
    irr <- get_power(
      community = "AG",
      lonlat = c(lake_coord_i$longitude, lake_coord_i$latitude),
      pars = c("ALLSKY_SFC_SW_DWN"),  # surface shortwave radiation (MJ/m^2/day)
      dates = c(date_t.init, date_t.samp)
    ) 
    
    # mean daily irradiance (30-day period prior to samples collection)
    irr.df <- data.frame(irr) # convert to data frame
    
    # daily light exposure and intensity data:
    sun_df.30day <- 
      
      sun_df.c %>%
      
      mutate(day = day(time),
             month = month(time),
             year = year(time)) %>%
      
      # estimate daily means:
      group_by(year, month, day) %>%
      
      summarise(
        # Astronomical daylight (i.e., sun above 0°)
        # Terrain-corrected “actual sunlight hours”:
        tc_daylight_hours = sum(sun_above_terrain)/6, # divided by 6 to convert 10-min intervals in 1 hour intervals
        
        # Compute terrain-corrected daylight fraction:
        terrain_fraction =  sum(sun_above_terrain) / sum(elevation > 0)) %>%  # fraction of sunlight in the lake vs total sunlight 
      
      data.frame()
    
    # Merge daylight and irradiance datasets
    sun_irr_df.30day <- 
      
      sun_df.30day %>%
      
      left_join(., irr.df %>% select(-YYYYMMDD, -DOY), by = c("day" = "DD", "month" = "MM", "year" = "YEAR")) %>%
      
      # We now calculate the proportion of daily irradiance that the 
      # lake may have received after accounting for surrounding shading:
      mutate(tc_irr = ALLSKY_SFC_SW_DWN * terrain_fraction) # terrain-corrected irradiance (or light intensity)
    
    # final light availability dataset (following 30-day average):
    lake_d_df <- rbind(lake_d_df,
                       
                       data.frame(massif = lake_coord_i$massif,
                                  lake = lake_coord_i$lake,
                                  lon = lake_coord_i$longitude,
                                  lat = lake_coord_i$latitude,
                                  
                                  altitude = lake_coord_i$altitude_m,
                                  altitude.est = alt_est,
                                  
                                  maxdist = maxdist_d,
                                  
                                  light_hrs = mean(sun_irr_df.30day$tc_daylight_hours), # daylight hours
                                  light_irr_MJ.m2.d = mean(sun_irr_df.30day$ALLSKY_SFC_SW_DWN), # total mean daily irradiance
                                  light_irr.tc_MJ.m2.d = mean(sun_irr_df.30day$tc_irr)) %>% # mean terrain corrected irradiance
                         
                         # Conversion:
                         # 1 MJ/m2/day ~ 25.56 µmol photons/m2/s
                         mutate(light_irr.tc_PAR_umol.m2.s = light_irr.tc_MJ.m2.d * conv_factor))
    
    
  }
  
  # Store results from all lakes in a global dataset:
  lake_light_int.df <- rbind(lake_light_int.df, lake_d_df)
  
}

# Check dataset:
head(lake_light_int.df)

## [2.4] Compare PAR among lakes ####

# Check accuracy of DEM comparing:
# altitude information measured in the field and 
# vs
# altitude estimated from point location in DEM (altitude.est):
plot(altitude.est ~ altitude, data = lake_light_int.df)
abline(0, 1, lty = 2)
# It does seem accurate enough.

# Visualise results and decide at which maxdist value comparisons should be made:
ggplot(lake_light_int.df, aes(x = maxdist, y = light_irr.tc_PAR_umol.m2.s)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_log10() +
  facet_wrap(~lake, scales = "free_y")

# It looks like a (max) distance of 3000 m is the best
# compromise between too close and too far, as light values become 
# stable in all lakes for distances ≥ 3000 m:

ggplot(lake_light_int.df %>%
         dplyr::filter(maxdist >= 3000), aes(x = maxdist, y = light_irr.tc_PAR_umol.m2.s)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_log10() +
  facet_wrap(~lake, scales = "free_y")

# OK! Comparing estimated light availability among lakes (using maxdist ≥ 3000 m):

ggplot(lake_light_int.df %>%
         dplyr::filter(maxdist == 3000), aes(x = lake, y = light_irr.tc_PAR_umol.m2.s)) +
  geom_point(size = 3) +
  geom_hline(yintercept = mean(subset(lake_light_int.df, maxdist == 3000)$light_irr.tc_PAR_umol.m2.s),
             lty = 2)

# Check average light availability (in µmol/m^2/s within PAR (photosynthetic active radiation))
summary(subset(lake_light_int.df, maxdist == 3000)$light_irr.tc_PAR_umol.m2.s)
# min.: 198.7
# median: 491.8
# mean: 464.2
# max.: 548.6

# [3] References ####

# (1)
# Yu, X., Wu, Z., Jiang, W., & Guo, X. (2015). 
# Predicting daily photosynthetically active radiation from global solar radiation in the Contiguous United States. 
# Energy Conversion and management, 89, 71-82.
# 10.1016/j.enconman.2014.09.038
#
# https://www.sciencedirect.com/science/article/pii/S0196890414008395?casa_token=oL-I0u-Z4tAAAAAA:n4Wz05Kok3_YyCajGLiG-6Jg3j5sclnMIl_8EDw6MvZJ71xdDgxb1mlvUKCTuI_VpCXl1b_v2iY
#
# PAR/solar irradiance is ca. 0.45-0.50 in temperate zones

# (2)
# Jones, H. G., Archer, N., Rotenberg, E., & Casa, R. (2003). 
# Radiation measurement for plant ecophysiology. 
# Journal of Experimental Botany, 54(384), 879-889.
# https://doi.org/10.1093/jxb/erg116
#
# https://academic.oup.com/jxb/article/54/384/879/631226?login=false
#
# 1 J = 4.6 µmol photons at 550 nm (green light)

# (3)
# Wang, W., Cai, S., Huang, J., Ding, R., & Chen, L. (2024). 
# Variation in the Quanta-to-Energy Ratio of Photosynthetically Active Radiation under the Cloudless Atmosphere. 
# Atmosphere, 15(10), 1166.
# https://doi.org/10.3390/atmos15101166
#
# https://www.mdpi.com/2073-4433/15/10/1166
# 
# Review: 1 J = 4.43-4.97 µmol photons
#
#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/light_availability.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################
