# Part 1 of 4

library(maptools)
library(gstat)
library(rgdal)
library(sp)
library(spacetime)
library(raster)
library(xts)
library(reshape2)

setwd("~/Projects/spatiotemporal_kriging")

# Set environment timezone to UTC to avoid timezone conflicts
Sys.setenv(TZ="UTC")

# http://spatialreference.org/ref/esri/usa-contiguous-lambert-conformal-conic/
crs_lcc <- CRS(paste0("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0",
                      "+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"))

# Read in temperature observations and transform year column to POSIXct
df_tmin <- read.csv('data/ann_anoms_tmin.csv',stringsAsFactors = FALSE)
df_tmax <- read.csv('data/ann_anoms_tmax.csv',stringsAsFactors = FALSE)

# Transform the dates in POSIXct dates
df_tmin$year <- as.POSIXct(paste(df_tmin$year, "01", "01", sep = "-"))
df_tmax$year <- as.POSIXct(paste(df_tmax$year, "01", "01", sep = "-"))

# Extract out unique station locations and build SpatialPointsDataFrame of 
# station locations and metadata. There should be 1218 total stations
stns <- unique(df_tmin[,c('station_id','station_name','longitude','latitude','elevation')])
# Change stns to SpatialPointsDataFrame
coordinates(stns) <- ~longitude+latitude
proj4string(stns) <- CRS("+proj=longlat +datum=WGS84")

# Project station locations to LCC
stns <- spTransform(stns, crs_lcc)

# Set the row names as the station ids. This will let us reference a station point
# by id in the STFDF
row.names(stns) <- stns$station_id

# Change obs_tmin_df and obs_tmax_df from long tables to a space-wide tables using dcast
# from the reshape2 package. In a space-wide table, each column is a station
# time series
obs_tmin_spacewide <- dcast(df_tmin[,c('year','station_id','tmin')],year~station_id)
obs_tmax_spacewide <- dcast(df_tmax[,c('year','station_id','tmax')],year~station_id)

# Separate out time index from the dataframe and save in this format for later with the CV
time_index <- unique(df_tmin$year)
obs_tmin_spacewide <- obs_tmin_spacewide[,-1]
obs_tmax_spacewide <- obs_tmax_spacewide[,-1]

# Contruct the stfdf objects for tmin and tmax
stidf_tmin <- stConstruct(x=df_tmin,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stidf_tmax <- stConstruct(x=df_tmax,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmin <- as(stidf_tmin, "STFDF")
stfdf_tmax <- as(stidf_tmax, "STFDF")
stfdf_tmin <- spTransform(stfdf_tmin, crs_lcc)
stfdf_tmax <- spTransform(stfdf_tmax, crs_lcc)
