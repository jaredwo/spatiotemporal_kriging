# Part 1 of 4

library(maptools)
library(gstat)
library(rgdal)
library(sp)
library(spacetime)
library(raster)
library(xts)
library(reshape2)

# Set environment timezone to UTC to avoid timezone conflicts
Sys.setenv(TZ="UTC")

# Definition of CONUS Lambert Conformal Conic Projection by which input stations
# and interpolation grid will be transformed
# http://spatialreference.org/ref/esri/usa-contiguous-lambert-conformal-conic/
crs_lcc <- CRS(paste0("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0",
                      "+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"))

# Read in temperature observations and transform year column to POSIXct
obs_tmin_df <- read.csv('data/ann_anoms_tmin.csv',stringsAsFactors = FALSE)
obs_tmin_df$year <- as.POSIXct(paste(obs_tmin_df$year, "01", "01", sep = "-"))

obs_tmax_df <- read.csv('data/ann_anoms_tmax.csv',stringsAsFactors = FALSE)
obs_tmax_df$year <- as.POSIXct(paste(obs_tmax_df$year, "01", "01", sep = "-"))

# Extract out unique station locations and build SpatialPointsDataFrame of 
# station locations and metadata. There should be 1218 total stations
stns <- unique(obs_tmin_df[,c('station_id','station_name','longitude','latitude','elevation')])
# Change stns to SpatialPointsDataFrame
coordinates(stns) <- ~longitude+latitude
proj4string(stns) <- CRS("+proj=longlat +datum=WGS84")

# Transform stations from lon/lat WGS84 to LCC
stns <- spTransform(stns, crs_lcc)
# Set the row names as the station ids. This will let us reference a station point
# by id in the STFDF
row.names(stns) <- stns$station_id
plot(stns)

# Change obs_tmin_df and obs_tmax_df from long tables to a space-wide tables using dcast
# from the reshape2 package. In a space-wide table, each column is a station
# time series
obs_tmin_df <- dcast(obs_tmin_df[,c('year','station_id','tmin')],year~station_id)
obs_tmax_df <- dcast(obs_tmax_df[,c('year','station_id','tmax')],year~station_id)

# Separate out time index from the dataframe
time_index <- obs_tmin_df$year
obs_tmin_df <- obs_tmin_df[,-1]
obs_tmax_df <- obs_tmax_df[,-1]

# Save the spacewide data frame for later in the cross validation
obs_tmin_spacewide <- obs_tmin_df
obs_tmax_spacewide <- obs_tmax_df
 
# Make sure stns is in the same order as columns of obs_tmin_df
stns <- stns[colnames(obs_tmin_df),]

# For input to STFDF,obs_tmin_df needs n*m rows with the spatial index moving
# the fastest
obs_tmin_df <- data.frame('tmin'=as.numeric(t(obs_tmin_df)))
obs_tmax_df <- data.frame('tmax'=as.numeric(t(obs_tmax_df)))

# Create STFDF objects
obs_tmin_spdf <- STFDF(sp=stns, time=time_index, data=obs_tmin_df)
obs_tmax_spdf <- STFDF(sp=stns, time=time_index, data=obs_tmax_df)

# Create dummy dates so that each year is considered a day in the variogram
# modeling. This allows for the ST anisitropy to be calculated correctly 
dummy_dates <- as.POSIXct(seq(as.Date('2015-01-01'), by='day',
                              length.out=length(obs_tmin_spdf@time)))
# Replace dates with dummy dates in the obs_tmin STFDF object
index(obs_tmin_spdf@time) <- dummy_dates
index(obs_tmax_spdf@time) <- dummy_dates

# Set dummy end dates on the obs_tmin STFDF object (times + one day)
obs_tmin_spdf@endTime <- (dummy_dates + (60*60*24))
obs_tmax_spdf@endTime <- (dummy_dates + (60*60*24))