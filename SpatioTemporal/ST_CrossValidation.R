# Part 4 of 4

#Cross Validation
load("~/Projects/spatiotemporal_kriging/Spatial/ST_Data.RData")

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
  temp.krige.tmin <- krigeST(tmin~1, data=obs_temp_tmin, modelList=sumMetric_tmin, newdata=stn_stf, nmax = 100,  stAni=tmin_stAni_day) 
  temp.krige.tmax <- krigeST(tmax~1, data=obs_temp_tmax, modelList=sumMetric_tmax, newdata=stn_stf, nmax = 100,  stAni=tmax_stAni_day)
  
  # Store for Leave One Out CV
  #One row = all stations for one year
  prediction_tmin[stnid] <- temp.krige.tmin$var1.pred
  prediction_tmax[stnid] <- temp.krige.tmax$var1.pred
}

predictionTimes <- prediction_tmin[1]

prediction_tmin <- prediction_tmin[,-1]
prediction_tmax <- prediction_tmax[,-1]

# Create a matrix of differences
differences_tmin <- prediction_tmin - obs_tmin_spacewide
differences_tmax <- prediction_tmax - obs_tmax_spacewide

biasDifferences_tmin <- differences_tmin
biasDifferences_tmax <- differences_tmax

# Get the absolute value for each station
differences_tmin <- abs(differences_tmin)
differences_tmax <- abs(differences_tmax)

# Apply mean for MAE
MAE_tmin <- apply(differences_tmin, 2, mean, na.rm = TRUE)
MAE_tmax <- apply(differences_tmax, 2, mean, na.rm = TRUE)

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

avgMSE_tmin <- mean(MSE_tmin)
avgMSE_tmax <- mean(MSE_tmax)

# For RMSE
# Square of MSE
RMSE_tmin <- MSE_tmin^(1/2)
RMSE_tmax <- MSE_tmax^(1/2)

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

# Plot spatial trends of errors
spplot(tmin_results, "MSE", main = "Tmin - MSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmax_results, "MSE", main = "Tmax - MSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmin_results, "MAE", main = "Tmin - MAE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmax_results, "MAE", main = "Tmax - MAE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmin_results, "RMSE", main = "Tmin - RMSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))
spplot(tmax_results, "RMSE", main = "Tmax - RMSE", cuts = c(0,.25,.5,1,1.75, 2.75,4))

# Create multi-panel plot of station observations for years 2012-2015
stplot(differences_tmin_spdf[,'2012/2015','tmin'])
stplot(differences_tmax_spdf[,'2012/2015','tmin'])

# Look at a single station's time series of error (chose 2 stations with large tmin error - Laketown, UT and Odessa, WA)
plot(biasDifferences_tmin$USH00424856~time_index, main = "Laketown, UT Tmin prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)
plot(biasDifferences_tmin$USH00456039~time_index, main = "Odessa, WA Tmin prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)

plot(biasDifferences_tmax$USH00424856~time_index, main = "Laketown, UT Tmax prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)
plot(biasDifferences_tmax$USH00456039~time_index, main = "Odessa, WA Tmax prediction error", xlab = "Year", ylab = Error in degrees Fahrenheit)

#Look at summary statistics
avgMAE_tmin
avgMAE_tmax
avgMSE_tmin
avgMSE_tmax
avgRMSE_tmin
avgRMSE_tmax

save.image(file = "~/Projects/spatiotemporal_kriging/SpatioTemporal/CV_ST.RData")