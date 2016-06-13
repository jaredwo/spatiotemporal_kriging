# Basic example of using gstat to interpolate point observations

# Load required libraries
# If using miniconda version of R, need to set environmental variable pointing 
# rgdal to proj4 share location
# Sys.setenv(PROJ_LIB=paste0(Sys.getenv("HOME"),"/miniconda2/envs/r/share/proj"))
library(rgdal)
library(raster)
library(gstat)
library(sp)

# Load station data
# Contains monthly Tmax normals (1981-2010)
stns <- read.csv("data/tmax_stns.csv", stringsAsFactors=FALSE)

# Inspect data.frame structure
str(stns)

# Change to SpatialPointsDataFrame
coordinates(stns) <- ~lon+lat
proj4string(stns) <- CRS("+proj=longlat +datum=WGS84")

# Plot station points
plot(stns)

# Read in domain extent for Pennsylvania (PA)
bbox <- extent(readOGR('data/bbox_pa/bbox_fnl.shp','bbox_fnl'))
# Add domain bbox to station plot
plot(bbox, col='orange', add=TRUE)

# Subset stations to PA domain
stns <- crop(stns, bbox)

# Re-plot PA stations in red
plot(stns, col='red', add=TRUE)

# Load PRISM DEM raster (4km resolution)
dem <- raster('data/PRISM_us_dem_4km_tif.tif')
# Crop to domain bbox
dem <- crop(dem, bbox)
# Make sure projection is the same as stns.
# The raster is technically NAD83, but we can consider
# it equal to WGS84 for our purposes
proj4string(dem) <- proj4string(stns)
# Convert to SpatialPixelsDataFrame to define grid
# for interpolations
dem_spdf <- as(dem,'SpatialPixelsDataFrame')
# Plot the DEM
spplot(dem_spdf)

##############################
# Examples for interpolating tmax08 variable:
# 1981-2010 August normal Tmax
##############################

# Basic inverse distance weighting
tmax_idw <- idw(tmax08~1,stns,dem_spdf)

# Plot results
spplot(tmax_idw['var1.pred'])

# Ordinary Kriging
# Create sample variogram
a_vgm_okrig <- variogram(tmax08~1, stns)
# Fit empirical variogram
a_vgm_fit_okrig <- fit.variogram(a_vgm_okrig, vgm("Sph"))
# Plot sample and fit variogram
plot(a_vgm_okrig,a_vgm_fit_okrig)
# Perform ordinary kriging
tmax_okrig <- krige(tmax08~1, stns, dem_spdf, model=a_vgm_fit_okrig)
# Plot result predictions and associated variance
spplot(tmax_okrig['var1.pred'])
spplot(tmax_okrig['var1.var'])

# Universal/KED kriging (tmax08~lon+lat+elev)
# Setup interpolation grid with required predictors
names(dem_spdf) <- c('elev')
dem_spdf@data['lon'] <- dem_spdf@coords[,1]
dem_spdf@data['lat'] <- dem_spdf@coords[,2]

# Create sample variogram
a_vgm_kedkrig <- variogram(tmax08~elev+lon+lat, stns)
# Fit empirical variogram
a_vgm_fit_kedkrig <- fit.variogram(a_vgm_kedkrig, vgm("Sph"))
# Plot sample and fit variogram
plot(a_vgm_kedkrig, a_vgm_fit_kedkrig)

# Perform universal/KED kriging
tmax_kedkrig <- krige(tmax08~elev+lon+lat, stns, dem_spdf, model=a_vgm_fit_kedkrig)
# Plot result predictions and associated variance
spplot(tmax_kedkrig['var1.pred'])
spplot(tmax_kedkrig['var1.var'])

# Compare interpolation grids from IDW, oridinary kriging, and KED
names(tmax_kedkrig) <- c("ked.pred","ked.var")
tmax_kedkrig@data['idw.pred'] <- tmax_idw$var1.pred
tmax_kedkrig@data['okrig.pred'] <- tmax_okrig$var1.pred
spplot(tmax_kedkrig[c('idw.pred', 'okrig.pred', 'ked.pred')])

# Run leave-one-out cross validation for ordinary kriging and KED
# These with both take several minutes
cv_kedkrig <- krige.cv(tmax08~elev+lon+lat, stns, a_vgm_fit_kedkrig)
cv_okrig <- krige.cv(tmax08~1, stns, a_vgm_fit_okrig)

# Compare mean absolute error
mean(abs(cv_kedkrig$residuals))
mean(abs(cv_okrig$residuals))
