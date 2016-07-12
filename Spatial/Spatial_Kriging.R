# Part 3 of 4

#Create a spatial object from the tmin and tmax data tables using long lat
coordinates(df_tmin) <- ~longitude+latitude
proj4string(df_tmin) <- CRS("+proj=longlat +datum=WGS84")
coordinates(df_tmax) <- ~longitude+latitude
proj4string(df_tmax) <- CRS("+proj=longlat +datum=WGS84")

#Project to LCC
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
spplot(tminPredict_pooled, "var1.pred", main = "Tmin: Predicted Annomalies 1995")

tmaxPredict_pooled <- krige(tmax~1, df_tmax_1995, spGrid_spdf, fitModel_tmax) 
spplot(tmaxPredict_pooled, "var1.pred", main = "Tmax: Predicted Annomalies 1995")