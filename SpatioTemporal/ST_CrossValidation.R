# Part 4 of 4

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