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

# Make sure stns is in the same order as columns of obs_tmin_df
stns <- stns[colnames(obs_tmin_df),]

# For input to STFDF,obs_tmin_df needs n*m rows with the spatial index moving
# the fastest
obs_tmin_df <- data.frame('tmin'=as.numeric(t(obs_tmin_df)))
obs_tmax_df <- data.frame('tmax'=as.numeric(t(obs_tmax_df)))

# Create STFDFs
obs_tmin_spdf <- STFDF(sp=stns, time=time_index, data=obs_tmin_df)
obs_tmax_spdf <- STFDF(sp=stns, time=time_index, data=obs_tmax_df)

# Create dummy dates so that each year is considered a day in the variogram
# modeling. Better way to do this?
dummy_dates <- as.POSIXct(seq(as.Date('2015-01-01'), by='day',
                              length.out=length(obs_tmin_spdf@time)))
# Set dummy dates on the obs_tmin STFDF object
index(obs_tmin_spdf@time) <- dummy_dates
index(obs_tmax_spdf@time) <- dummy_dates

# Set dummy end dates on the obs_tmin STFDF object (times + one day)
obs_tmin_spdf@endTime <- (dummy_dates + (60*60*24))
obs_tmax_spdf@endTime <- (dummy_dates + (60*60*24))
                          
# Set environment timezone to UTC to avoid timezone conflicts
Sys.setenv(TZ="UTC")

#Rename for convenience
stfdf_tmin_day <- obs_tmin_spdf
stfdf_tmax_day <- obs_tmax_spdf

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
plot(tmin_var_day[tmin_var_day$spacelag==0,'timelag'],tmin_var_day[tmin_var_day$spacelag==0,'gamma'])
#Tmax
# Spatial: Nugget ~ .15 range ~ 2000km sill ~.8
# Temporal: Nugget ~ .18 range ~ 3 years sill ~ .7
# Joint: sill = .9
plot(tmax_var_day,map=F)
plot(tmax_var_day,wireframe=T) 
plot(tmax_var_day[tmax_var_day$spacelag==0,'timelag'],tmax_var_day[tmax_var_day$spacelag==0,'gamma'])

#Estimate ST Anisotropy
tmin_stAni_day <- estiStAni(tmin_var_day, interval = c(10, 1000), method = "linear")
tmax_stAni_day <- estiStAni(tmax_var_day, interval = c(10, 1000), method = "linear")

# Sum Metric Model
#Tmin
#Using parameters from earlier model estimations
sumMetricModel_tmin1_day <- vgmST("sumMetric",
                                  space = vgm(psill=0, "Exp", range=1500, nugget=0),
                                  time = vgm(psill=0,"Exp", range=1, nugget=0),
                                  joint = vgm(psill=0.53, "Sph", range=3000, nugget=0.17 ),
                                  stAni = tmin_stAni_day)

# Does this fit?
plot(tmin_var_day, sumMetricModel_tmin1_day,map=F, all=T)

# Try Fitting
fitSumMetric_tmin1_day <- fit.StVariogram(tmin_var_day, sumMetricModel_tmin1_day,
                                          control = list(parscale = c(1,10000,1,1,1,1,1,10000,1,10000), maxit=2e4))
# parscale: sill.s, range.s, nugget.s, sill.t ,  range.t ,   nugget.t , sill.st, range.st , nugget.st, ani

# Try the various fit methods and choose the one with the lowest MSE
for (variable in 1:13){
  temp.alternate <- fit.StVariogram(tmin_var_day, sumMetricModel_tmin1_day, fit.method=variable, method="L-BFGS-B", stAni=tmin_stAni_day)
  plot(tmin_var_day, temp.alternate, map=F, all=T, main="Fit method = " + variable)
  temp.mse = attr(temp.alternate, "MSE")
  main.mse = attr(fitSumMetric_tmin1_day, "MSE")
#   if(temp.mse<main.mse){
#     fitSumMetric_tmin1_day <- temp.alternate
#     print("New best " + variable)
#   }
}

# Does this fit better?
plot(tmin_var_day, fitSumMetric_tmin1_day,map=F, all=T) #yes

#Best Model
sumMetric_tmin_day <- fitSumMetric_tmin1_day

#Tmax
sumMetricModel_tmax1_day <- vgmST("sumMetric",
                                  space = vgm(psill=0, "Exp", range=1500, nugget=0),
                                  time = vgm(psill=0,"Exp", range=1, nugget=0),
                                  joint = vgm(psill=0.68, "Sph", range=3000, nugget=0.16),
                                  stAni = tmax_stAni_day)

# Does this fit?
plot(tmax_var_day, sumMetricModel_tmax1_day,map=F, all=T)

# Try Fitting
fitSumMetric_tmax1_day <- fit.StVariogram(tmax_var_day, sumMetricModel_tmax1_day,
                                          control = list(parscale = c(1,10000,1,1,1,1,1,10000,1,10000), maxit=2e4))

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
plot(tmax_var_day, fitSumMetric_tmax1_day,map=F, all=T) #Yes

#Best
sumMetric_tmax_day <- fitSumMetric_tmax1_day

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
# stfdf_tmin[station, time, attribute] how to do???
# SUbset out one station, predict time series, repeat for all stations and combine predictions in data frame
# Do same for pooled (taking out on station), but predict each year for all stations for all years

### Create grids for storing the predicted time series
prediction_tmin <- data.frame(time=unique(index(obs_tmin_spdf)))
prediction_tmax <- data.frame(time=unique(index(obs_tmin_spdf)))

for (i in 1:length(stns)){
  #The station to be kriged as an STF object
  stnid <- stns$station_id[i]
  stn_stf <- as(obs_tmin_spdf[stnid,,'tmin',drop=F],'STF')
  
  # The stations to be used in the kriging (all except the station to be kriged)
  stnids_keep <- stns$station_id[-i]
  
  # Create a temporary SPDF of all stations except the one to be kriged
  obs_temp_tmin <- obs_tmin_spdf[stnids_keep]
  obs_temp_tmax <- obs_tmax_spdf[stnids_keep]
  
  # The kriged time series
  temp.krige.tmin <- krigeST(tmin~1, data=obs_temp_tmin, modelList=sumMetric_tmin_day, newdata=stn_stf, nmax = 100,  stAni=tmin_stAni_day) 
  temp.krige.tmax <- krigeST(tmax~1, data=obs_temp_tmax, modelList=sumMetric_tmax_day, newdata=stn_stf, nmax = 100,  stAni=tmax_stAni_day)
  
  # Store for Leave One Out CV
  prediction_tmin[stnid] <- temp.krige.tmin$var1.pred
  prediction_tmax[stnid] <- temp.krige.tmax$var1.pred
}
# Bibliography
#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/