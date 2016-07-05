# Example of creating a "pooled" spatial variogram from a spatiotemporal dataset

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

df_tmin$year <- as.POSIXct(paste(df_tmin$year, "01", "01", sep = "-"))
df_tmax$year <- as.POSIXct(paste(df_tmax$year, "01", "01", sep = "-"))

# Extract out unique station locations and build SpatialPointsDataFrame of 
# station locations and metadata. There should be 1218 total stations
stns <- df_tmin
# Change stns to SpatialPointsDataFrame
coordinates(stns) <- ~longitude+latitude
proj4string(stns) <- CRS("+proj=longlat +datum=WGS84")
# Transform stations from lon/lat WGS84 to LCC
stns <- spTransform(stns, crs_lcc)
# Set the row names as the station ids. This will let us reference a station point
# by id in the STFDF
row.names(stns) <- stns$station_id
plot(stns)

stidf_tmin <- stConstruct(x=df_tmin,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmin <- as(stidf_tmin, "STFDF")
stidf_tmax <- stConstruct(x=df_tmax,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmax <- as(stidf_tmax, "STFDF")
stfdf_tmin <- spTransform(stfdf_tmin, crs_lcc)
stfdf_tmax <- spTransform(stfdf_tmax, crs_lcc)

# Calculate a "pooled" spatial only variogram by setting the tlags param = 0
# Set cutoff to 3500 km to make sure we get the range
emp_varst_tmin <- variogramST(tmin~1, data=stfdf_tmin, tlags=0,
                              progress=TRUE, assumeRegular=TRUE,
                              cutoff=3500, width=3500/30)
emp_varst_tmax <- variogramST(tmax~1, data=stfdf_tmax, tlags=0,
                              progress=TRUE, assumeRegular=TRUE,
                              cutoff=3500, width=3500/30)
# Plot pooled spatial sample variogram
plot(emp_varst_tmin, map=F)
plot(emp_varst_tmax, map=F)

# To use this sample variogram with pure spatial gstat functions, we must convert it
# from a StVariogram object to a gstatVariogram object

# Function to convert spatiotemporal variogram to spatial variogram
stvario_to_svario <- function(stvario)
{
  svario <- stvario[stvario$timelag == 0,]
  class(svario) <- c("gstatVariogram", "data.frame")
  svario <- svario[-1,1:3]
  svario$dir.hor <- 0
  svario$dir.ver <- 0
  return(svario)
}

emp_var_tmin <- stvario_to_svario(emp_varst_tmin)
emp_var_tmax <- stvario_to_svario(emp_varst_tmax)

# Can now use emp_var_tmin/tmax with any pure spatial gstat functions
#Tmin model to begin with
vgmModel_tmin <- vgm(psill = .7, model = "Exp", range = 2000, nugget = .15)

#Does this fit?
plot(emp_var_tmin, vgmModel_tmin, all=T)

#Try Fitting
fitModel_tmin <- fit.variogram(emp_var_tmin, vgmModel_tmin)
plot(emp_var_tmin, fitModel_tmin, all=T)

#Tmax model to begin with
vgmModel_tmax <- vgm(psill = 1.1, model = "Exp", range = 2000, nugget = .11)

#Does this fit?
plot(emp_var_tmax, vgmModel_tmax, all=T)

#Try Fitting
fitModel_tmax <- fit.variogram(emp_var_tmax, vgmModel_tmax)
plot(emp_var_tmax, fitModel_tmax, all=T)

#Pure Spatial Data Frame
coordinates(df_tmin) <- ~longitude+latitude
proj4string(df_tmin) <- CRS("+proj=longlat +datum=WGS84")
coordinates(df_tmax) <- ~longitude+latitude
proj4string(df_tmax) <- CRS("+proj=longlat +datum=WGS84")

#Change to LCC projection
df_tmin <- spTransform(df_tmin, crs_lcc)
df_tmax <- spTransform(df_tmax, crs_lcc)

#Grid
spGrid <- raster('data/mask_conus_25deg.tif')

# Set up spatial grid
proj4string(spGrid) <- proj4string(df_tmin)
spGrid_spdf <- as(spGrid,'SpatialPixelsDataFrame')
spplot(spGrid_spdf)

#Get just 1995 data
df_tmin_1995 <- stfdf_tmin[,"1995", 'tmin']
df_tmax_1995 <- stfdf_tmax[,"1995", 'tmax']

#Kriging
tminPredict_pooled <- krige(tmin~1,df_tmin_1995, spGrid_spdf, fitModel_tmin) 
spplot(tminPredict_pooled, main = "T-Min")

tmaxPredict_pooled <- krige(tmax~1, df_tmax_1995, spGrid_spdf, fitModel_tmax) 
spplot(tmaxPredict_pooled, main = "T-Max")

#Crossvalidation
# Make years in a vector
years <- as.character(seq(1900,2015,1))

#Set up Null data frames
pooledPrediction_tmin <- NULL
pooledPrediction_tmax <- NULL

bMask_tmin <- NULL
bMask_tmax <- NULL

  # Loop through years
for (y in 1:length(years)){
  # Get just that year's data for main grid
  df_tmin_year <- stfdf_tmin[,years[y], 'tmin']
  df_tmax_year <- stfdf_tmax[,years[y], 'tmax']
  
  # Get rid of NA stations
  bMask_tmin <- cbind(bMask_tmin, is.finite(df_tmin_year$tmin))
  bMask_tmax <- cbind(bMask_tmax, is.finite(df_tmax_year$tmax))
  df_tmin_year <- df_tmin_year[is.finite(df_tmin_year$tmin),]
  df_tmax_year <- df_tmax_year[is.finite(df_tmax_year$tmax),]

  # Do  leave one out CV
  tmin.cv.temp <- krige.cv(tmin~1,df_tmin_year, fitModel_tmin)
  tmax.cv.temp <- krige.cv(tmax~1,df_tmax_year, fitModel_tmax)
  
  # Create a vector with all predicted values and NA for non applicable stations
  stationsVector_tmin <- NULL
  stationsVector_tmax <- NULL
  
  #Set index variables to the start of the vector
  station_tmin = 1
  station_tmax = 1
    
  #Tmin
  for (i in 1:length(bMask_tmin)){ # For all stations if there is an observation store the predicted vale and increment the station index
    if (bMask_tmin[i] == TRUE) {
      stationsVector_tmin[i] = tmin.cv.temp$var1.pred[station_tmin]
      station_tmin = station_tmin+1
    }
    #Otherwise, set the predicted value as NA
  }
  
  #Repeat for Tmax
  for (i in 1:bMask_tmax){
    if (bMask_tmax[i] == TRUE) {
      stationsVector_tmax[i] = tmax.cv.temp$var1.pred[station_tmax]
      station_tmax = station_tmax+1
    }
  }
  
  # Store CV Prediction vectors in a data frame by row (1 row = 1 year and all stations)
  pooledPrediction_tmin <- rbind(pooledPrediction_tmin, stationsVector_tmin)
  pooledPrediction_tmax <- rbind(pooledPrediction_tmax, stationsVector_tmax)
}