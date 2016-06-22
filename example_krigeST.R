# Example of using gstat krigeST for spatiotemporal interpolation.

library(spacetime)
library(sp)
library(gstat)
library(xts)
library(raster)

# Set environment timezone to UTC to avoid timezone conflicts
Sys.setenv(TZ="UTC")

# Definition of CONUS Lambert Conformal Conic Projection by which input stations
# and interpolation grid will be transformed
# http://spatialreference.org/ref/esri/usa-contiguous-lambert-conformal-conic/
crs_lcc <- CRS(paste0("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0",
                      "+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"))

# Read in Tmin observations and transform year column to POSIXct
obs_tmin <- read.csv('data/ann_anoms_tmin.csv',stringsAsFactors = FALSE)
obs_tmin$year <- as.POSIXct(paste(obs_tmin$year, "01", "01", sep = "-"))

# Change to SpatialPointsDataFrame
coordinates(obs_tmin) <- ~longitude+latitude
proj4string(obs_tmin) <- CRS("+proj=longlat +datum=WGS84")

# Transfrom points from lon/lat WGS84 to LCC
obs_tmin <- spTransform(obs_tmin, crs_lcc)

# Convert back to dataframe for input to stConstruct
obs_tmin <- as.data.frame(obs_tmin)

# Convert to spacetime object. Note: spatial coordinates are now x and y, not
# longitude and latitude
obs_tmin <- stConstruct(x=obs_tmin, c('x','y'), time="year", crs = crs_lcc)
obs_tmin <- as(obs_tmin, "STFDF")

# Check out data: plot years 2012-2015
stplot(obs_tmin[,'2012/2015','tmin'])
       
# Create dummy dates so that each year is considered a day in the variogram
# modeling. Better way to do this?
dummy_dates <- as.POSIXct(seq(as.Date('2015-01-01'), by='day',
                              length.out=length(obs_tmin@time)))
# Set dummy dates on the obs_tmin STFDF object
index(obs_tmin@time) <- dummy_dates
# Set dummy end dates on the obs_tmin STFDF object (times + one day)
obs_tmin@endTime <- (dummy_dates + (60*60*24))

# Build empirical variogram going out to 5 time lags
emp_varst_tmin <- variogramST(tmin~1, data=obs_tmin, tlags=0:5,
                              progress=TRUE, assumeRegular=TRUE,
                              width=3500/30,
                              cutoff=3500)

# For this example, just fit a separable model
vgm_sep <- vgmST("separable",
                 space = vgm(.5, "Sph", 3000, 0.18),
                 time = vgm(.5, "Sph", 3, 0.19),
                 sill = 1)
vgm_sep <- fit.StVariogram(emp_varst_tmin, vgm_sep, fit.method=6,
                           method = "L-BFGS-B")

# Load interpolation raster grid
agrid <- raster('data/mask_conus_25deg.tif')
# Transform to LCC projection going from 0.25 deg resolution to 30 km resolution
agrid <- projectRaster(agrid, res=30, crs=crs_lcc, method='ngb')
# Convert to SpatialPixelsDataFrame to define grid
# for interpolations
agrid <- as(agrid,'SpatialPixelsDataFrame')

# Create spatialtemporal object to which we will interpolate
# In this case, we'll run interpolations for the first 3 dates
astf <- STF(agrid, dummy_dates[1:3])

# Estimate spatiotemproal anisotropy parameter
stani <- estiStAni(emp_varst_tmin, interval=c(0,500), method="linear")
vals <- krigeST(tmin~1, data=obs_tmin, newdata=astf, modelList=vgm_sep, nmax=100, stAni=stani)

# Plot results
stplot(vals)

