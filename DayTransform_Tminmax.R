library(maptools)
library(gstat)
library(rgdal)
library(sp)
library(spacetime)
library(raster)
library(xts)

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

obs_tmax <- read.csv('data/ann_anoms_tmax.csv', stringsAsFactors = FALSE)
obs_tmax$year <- as.POSIXct(paste(obs_tmax$year, "01", "01", sep = "-"))

# Change to SpatialPointsDataFrame
coordinates(obs_tmin) <- ~longitude+latitude
proj4string(obs_tmin) <- CRS("+proj=longlat +datum=WGS84")

coordinates(obs_tmax) <- ~longitude+latitude
proj4string(obs_tmax) <- CRS("+proj=longlat +datum=WGS84")

# Transfrom points from lon/lat WGS84 to LCC
obs_tmin <- spTransform(obs_tmin, crs_lcc)
obs_tmax <- spTransform(obs_tmax, crs_lcc)

# Convert back to dataframe for input to stConstruct
obs_tmin <- as.data.frame(obs_tmin)
obs_tmax <- as.data.frame(obs_tmax)

# Convert to spacetime object. Note: spatial coordinates are now x and y, not
# longitude and latitude
obs_tmin <- stConstruct(x=obs_tmin, c('x','y'), time="year", crs = crs_lcc)
obs_tmin <- as(obs_tmin, "STFDF")

obs_tmax <- stConstruct(x=obs_tmax, c('x','y'), time="year", crs = crs_lcc)
obs_tmax <- as(obs_tmax, "STFDF")

# Create dummy dates so that each year is considered a day in the variogram
# modeling. Better way to do this?
dummy_dates <- as.POSIXct(seq(as.Date('2015-01-01'), by='day',
                              length.out=length(obs_tmin@time)))
# Set dummy dates on the obs_tmin STFDF object
index(obs_tmin@time) <- dummy_dates
index(obs_tmax@time) <- dummy_dates

# Set dummy end dates on the obs_tmin STFDF object (times + one day)
obs_tmin@endTime <- (dummy_dates + (60*60*24))
obs_tmax@endTime <- (dummy_dates + (60*60*24))

#Rename for convenience
stfdf_tmin_day <- obs_tmin
stfdf_tmax_day <- obs_tmax

# Compute sample variogram
tmin_var_day <- variogramST(tmin~1,data=stfdf_tmin_day,tunit="days",assumeRegular=T,na.omit=T, progress = TRUE, cutoff=3500) 
tmax_var_day <- variogramST(tmax~1,data=stfdf_tmax_day,tunit="days",assumeRegular=T,na.omit=T, progress = TRUE, cutoff=3500) 

#Inspect the sample variograms in two forms in order to find estimates for fitting the variogram
#Tmin
# Spatial: Nugget ~ .2 range ~ 3000km sill ~.7
# Temporal: Nugget ~ .42 range ~ 3 years sill ~ .55
# Joint: sill = .7
plot(tmin_var_day,map=F) 
plot(tmin_var_day,wireframe=T)
#Tmax
# Spatial: Nugget ~ .15 range ~ 2000km sill ~.8
# Temporal: Nugget ~ .18 range ~ 3 years sill ~ .7
# Joint: sill = .9
plot(tmax_var_day,map=F)
plot(tmax_var_day,wireframe=T) 

#Estimate ST Anisotropy
tmin_stAni_day <- estiStAni(tmin_var_day, interval = c(10, 1000), method = "linear")
tmax_stAni_day <- estiStAni(tmax_var_day, interval = c(10, 1000), method = "linear")

#Estimate lower and upper bounds
bound.l.tmin_day <- c(sill.s = .4, range.s = 1300, nugget.s = .1,
                  sill.t = .4, range.t = 1, nugget.t = .1,
                  sill.st = .55, range.st = 1300, nugget.st = 0.1,
                  anis = 120)
bound.u.tmin_day <- c(sill.s = .7, range.s = 2000, nugget.s = .2,
                  sill.t = .8, range.t = 3, nugget.t = .2,
                  sill.st = 200, range.st = 2000, nugget.st = .2,
                  anis = 140)  

bound.l.tmax_day <- c(sill.s = .7, range.s = 1500, nugget.s = .1,
                  sill.t = .6, range.t = 1, nugget.t = .1,
                  sill.st = .9, range.st = 1000, nugget.st = .1,
                  anis = 140)
bound.u.tmax_day <- c(sill.s = .9, range.s = 2500, nugget.s = .2,
                  sill.t = .8, range.t = 3, nugget.t = .2,
                  sill.st = 1.1, range.st = 2000, nugget.st = .2,
                  anis = 150)   

# Sum Metric Model
#Tmin
#Using parameters from earlier model estimations
sumMetricModel_tmin1_day <- vgmST("sumMetric",
                              space = vgm(.45, "Sph", 3000, 0.11),
                              time = vgm(0.1335337, "Sph", 1.95, 0.2495566),
                              joint = vgm(.01,"Mat", 1500, .1),
                              stAni = tmin_stAni_day)

# Does this fit?
plot(tmin_var_day, sumMetricModel_tmin1_day,map=F, all=T)

# Try Fitting
fitSumMetric_tmin1_day <- fit.StVariogram(tmin_var_day, sumMetricModel_tmin1_day, fit.method = 8, #Fit method 8 seems like a good fit. (Weighting with distance and bin amount)
                                      stAni = tmin_stAni_day, method="L-BFGS-B")

# Try the various fit methods and choose the one with the lowest MSE
for (variable in c(1:13)){
  temp.alternate <- fit.StVariogram(tmin_var_day, fitSumMetric_tmin1_day, fit.method=variable, method="BFGS", stAni=tmin_stAni_day)
  temp.mse = attr(temp.alternate, "MSE")
  main.mse = attr(fitSumMetric_tmin1_day, "MSE")
  if(temp.mse<main.mse){
    fitSumMetric_tmin1_day <- temp.alternate
    print("New best " + variable)
  }
}

# Does this fit better?
plot(tmin_var_day, fitSumMetric_tmin1_day,map=F, all=T) #No

#Best Model
sumMetric_tmin_day <- sumMetricModel_tmin1_day

#Tmax
sumMetricModel_tmax1_day <- vgmST("sumMetric",
                              space = vgm(0.65, "Exp", 1400, 0.14),
                              time = vgm(.15, "Sph", 3, 0.15),
                              joint = vgm(.3,"Mat", 2500, 0),
                              stAni = tmax_stAni_day)

#Does this fit?
plot(tmax_var_day, sumMetricModel_tmax1_day,map=F, all=T)

#Try Fitting
fitSumMetric_tmax1_day <- fit.StVariogram(tmax_var_day, sumMetricModel_tmax1_day, fit.method = 8,
                                      stAni = tmax_stAni_day, method="L-BFGS-B")

# Try the various fit methods and choose the one with the lowest MSE
for (variable in c(1:13)){
  temp.alternate <- fit.StVariogram(tmax_var_day, fitSumMetric_tmax1_day, fit.method=variable, method="BFGS", stAni=tmax_stAni_day)
  temp.mse = attr(temp.alternate, "MSE")
  main.mse = attr(fitSumMetric_tmax1_day, "MSE")
  if(temp.mse<main.mse){
    fitSumMetric_tmax1_day <- temp.alternate
    print("New best " + variable)
  }
}

#Does this fit better?
plot(tmax_var_day, fitSumMetric_tmax1_day,map=F, all=T) #No

#Best
sumMetric_tmax_day <- sumMetricModel_tmax1_day

# Plot empirical vs. fitted variograms
plot(tmin_var_day, sumMetric_tmin_day, map=F,main = "SumMetric Tmin", all=T)
plot(tmax_var_day, sumMetric_tmax_day, map=F,main = "SumMetric Tmax", all=T)

# Differences
plot(tmin_var_day, sumMetric_tmin_day, wireframe=T,diff = T, main = "SumMetric Difference Tmin")
plot(tmax_var_day, sumMetric_tmax_day, wireframe=T,diff = T, main = "SumMetric Difference Tmax")

# MSE Measures
attr(sumMetric_tmin_day, "MSE")
attr(sumMetric_tmax_day, "MSE")

#Create prediction grid
tm.grid_day <- dummy_dates[96:97]
spGrid <- raster('data/mask_conus_25deg.tif')
spGrid <- projectRaster(spGrid, res=30, crs=crs_lcc, method='ngb')

# Convert to SpatialPixelsDataFrame to define grid
# for interpolations
spGrid_spdf <- as(spGrid,'SpatialPixelsDataFrame')

#Create the spatiotemporal grid object to interpolate to
grid.ST_day <- STF(spGrid_spdf,tm.grid_day)

#ST Kriging
precipPredict_tmin_day <- krigeST(tmin~1, data=stfdf_tmin_day, modelList=sumMetric_tmin_day, newdata=grid.ST_day, nmax = 100,  stAni=tmin_stAni_day) 
stplot(precipPredict_tmin_day)

precipPredict_tmax_day <- krigeST(tmax~1, data=stfdf_tmax_day, modelList=sumMetric_tmax_day, newdata=grid.ST_day, nmax = 100,  stAni=tmax_stAni_day) 
stplot(precipPredict_tmax_day)

#Cross Validation



# Bibliography
#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/