# Part 4 of 4
load("~/Projects/spatiotemporal_kriging/PooledWorkspace.RData")
library(maptools)
library(gstat)
library(rgdal)
library(sp)
library(spacetime)
library(raster)
library(xts)
library(reshape2)

#Crossvalidation
# Make years in a vector of characters
years <- as.character(seq(1900,2015,1))

#Set up Null data frames for the storage for the CV predicted values
pooledPrediction_tmin <- NULL
pooledPrediction_tmax <- NULL

#Do the same for a boolean mask for which stations have an observation in each specific year
bMask_tmin <- NULL
bMask_tmax <- NULL

# Loop through each year
for (y in 1:length(years)){
  # Extract just that year's data
  df_tmin_year <- stfdf_tmin[,years[y], 'tmin']
  df_tmax_year <- stfdf_tmax[,years[y], 'tmax']
  
  # Get rid of NA stations but keep the boolean mask for which station were removed
  bMask_tmin <- is.finite(df_tmin_year$tmin)
  bMask_tmax <- is.finite(df_tmax_year$tmax)
  df_tmin_year <- df_tmin_year[bMask_tmin,]
  df_tmax_year <- df_tmax_year[bMask_tmax,]
  
  # Do  leave one out CV with all stations with observation from that year
  tmin.cv.temp <- krige.cv(tmin~1,df_tmin_year, fitModel_tmin)
  tmax.cv.temp <- krige.cv(tmax~1,df_tmax_year, fitModel_tmax)
  
  # Create a vector with all predicted values for that year and NA for stations without observations
  stationsVector_tmin <- rep(NA, length(bMask_tmin))
  stationsVector_tmax <- rep(NA, length(bMask_tmax))
  stationsVector_tmin[bMask_tmin] <- tmin.cv.temp$var1.pred
  stationsVector_tmax[bMask_tmax] <- tmax.cv.temp$var1.pred
  
  # Store CV Prediction vectors in a data frame by row (1 row = 1 year and all stations)
  pooledPrediction_tmin <- rbind(pooledPrediction_tmin, stationsVector_tmin)
  pooledPrediction_tmax <- rbind(pooledPrediction_tmax, stationsVector_tmax)
  print(y)
}

# Create a data frame so that each station's predictions can be accesses via the station ID
rownames(pooledPrediction_tmin) <- years
colnames(pooledPrediction_tmin) <- stns$station_id
rownames(pooledPrediction_tmax) <- years
colnames(pooledPrediction_tmax) <- stns$station_id

# Create a matrix of differences between predicted and observed values
differences_tmin <- pooledPrediction_tmin - obs_tmin_spacewide
differences_tmax <- pooledPrediction_tmax - obs_tmax_spacewide

# Get the absolute value of the differences
differences_tmin <- abs(differences_tmin)
differences_tmax <- abs(differences_tmax)

# Apply mean for MAE
MAE_tmin <- apply(differences_tmin, 2, mean, na.rm = TRUE)
MAE_tmax <- apply(differences_tmax, 2, mean, na.rm = TRUE)
#Get a global average
avgMAE_tmin <- mean(MAE_tmin)
avgMAE_tmax <- mean(MAE_tmax)

# For MSE
# Square differences_tmin^2)
squareDiff_tmin <- differences_tmin^2
squareDiff_tmax <- differences_tmax^2
# Sums
MSE_tmin <- apply(squareDiff_tmin, 2, sum, na.rm = TRUE)
MSE_tmax <- apply(squareDiff_tmax, 2, sum, na.rm = TRUE)
# Divide by number of non NA
numNA_tmin <- apply(differences_tmin, 2, function(x) sum(is.finite(x)))
numNA_tmax <- apply(differences_tmax, 2, function(x) sum(is.finite(x)))
MSE_tmin <- MSE_tmin/numNA_tmin
MSE_tmax <- MSE_tmax/numNA_tmax

# Global average
avgMSE_tmin <- mean(MSE_tmin)
avgMSE_tmax <- mean(MSE_tmax)

# For RMSE
# Square of MSE
RMSE_tmin <- MSE_tmin^(1/2)
RMSE_tmax <- MSE_tmax^(1/2)

# Global average
avgRMSE_tmin <- mean(RMSE_tmin)
avgRMSE_tmax <- mean(RMSE_tmax)

# Create data frames for spatial trends of error
differences_tmin_spdf <- STFDF(sp=stns, time=time_index, data=data.frame('differences'=as.numeric(t(differences_tmin))))
differences_tmax_spdf <- STFDF(sp=stns, time=time_index, data=data.frame('differences'=as.numeric(t(differences_tmax))))
tmin_results <- stns
tmax_results <- stns
tmin_results$MSE <- MSE_tmin
tmax_results$MSE <- MSE_tmax
tmin_results$MAE <- MAE_tmin
tmax_results$MAE <- MAE_tmax
tmin_results$RMSE <- RMSE_tmin
tmax_results$RMSE <- RMSE_tmax

# Plot spatial trends of errors
spplot(tmin_results, "MSE", main = "Tmin - MSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmax_results, "MSE", main = "Tmax - MSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmin_results, "MAE", main = "Tmin - MAE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmax_results, "MAE", main = "Tmax - MAE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmin_results, "RMSE", main = "Tmin - RMSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmax_results, "RMSE", main = "Tmax - RMSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))

# Create multi-panel plot of station observations for years 2012-2015
stplot(differences_tmin_spdf[,'2012/2015','differences'])
stplot(differences_tmax_spdf[,'2012/2015','differences'])

# Look at a single station's time series of error (chose 2 stations with large tmin error - Laketown, UT and Odessa, WA)
plot(biasDifferences_tmin$USH00424856~time_index, main = "Laketown, UT Tmin prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)
plot(biasDifferences_tmin$USH00456039~time_index, main = "Odessa, WA Tmin prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)

plot(biasDifferences_tmax$USH00424856~time_index, main = "Laketown, UT Tmax prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)
plot(biasDifferences_tmax$USH00456039~time_index, main = "Odessa, WA Tmax prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit

# Looke at summary statistics     
avgMAE_tmin
avgMAE_tmax
avgMSE_tmin
avgMSE_tmax
avgRMSE_tmin
avgRMSE_tmax

save.image(file = "~/Projects/spatiotemporal_kriging/Spatial/CV_S.RData")