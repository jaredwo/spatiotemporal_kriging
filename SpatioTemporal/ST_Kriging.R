# Part 3 of 4

#Create prediction grid for 1995 1996
tm.grid_day <- dummy_dates[96:97]
spGrid <- raster('data/mask_conus_25deg.tif')
spGrid <- projectRaster(spGrid, res=30, crs=crs_lcc, method='ngb')

# Convert to SpatialPixelsDataFrame to define grid
# for interpolations
spGrid_spdf <- as(spGrid,'SpatialPixelsDataFrame')

#Create the spatiotemporal grid object to interpolate to
grid.ST_day <- STF(spGrid_spdf,tm.grid_day)

#Perform spatio temporal kriging and plot results
precipPredict_tmin_day <- krigeST(tmin~1, data=stfdf_tmin_day, modelList=sumMetric_tmin, newdata=grid.ST_day, nmax = 100,  stAni=tmin_stAni_day) 
stplot(precipPredict_tmin_day, main = "Tmin: Predicted Variances for 1995 1996")

precipPredict_tmax_day <- krigeST(tmax~1, data=stfdf_tmax_day, modelList=sumMetric_tmax, newdata=grid.ST_day, nmax = 100,  stAni=tmax_stAni_day) 
stplot(precipPredict_tmax_day,  main = "Tmax: Predicted Variances for 1995 1996")