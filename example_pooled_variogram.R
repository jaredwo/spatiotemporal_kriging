# Example of creating a "pooled" spatial variogram from a spatiotemporal dataset

library(spacetime)
library(sp)
library(gstat)
library(raster)
library(rgdal)

setwd("~/Projects/spatiotemporal_kriging")
df_tmin <- read.csv('data/ann_anoms_tmin.csv',stringsAsFactors = FALSE)
df_tmax <- read.csv('data/ann_anoms_tmax.csv',stringsAsFactors = FALSE)
df_tmin$year <- as.POSIXct(paste(df_tmin$year, "01", "01", sep = "-"), tz='GMT')
df_tmax$year <- as.POSIXct(paste(df_tmax$year, "01", "01", sep = "-"), tz='GMT')
stidf_tmin <- stConstruct(x=df_tmin,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmin <- as(stidf_tmin, "STFDF")
stidf_tmax <- stConstruct(x=df_tmax,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmax <- as(stidf_tmax, "STFDF")

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

# Can now use emp_var_tmin with any pure spatial gstat functions
vgmModel_tmin <- vgm(psill = .7, model = "Exp", range = 2000, nugget = .15)

#Does this fit?
plot(emp_var_tmin, vgmModel_tmin, all=T)

#Try Fitting
fitModel_tmin <- fit.variogram(emp_var_tmin, vgmModel_tmin)
plot(emp_var_tmin, fitModel_tmin, all=T)

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
precipPredict_tmin <- krige(tmin~1,df_tmin_1995, spGrid_spdf, fitModel_tmin) 
spplot(precipPredict_tmin, main = "T-Min")

precipPredict_tmax <- krige(tmax~1, df_tmax_1995, spGrid_spdf, fitModel_tmax) 
spplot(precipPredict_tmax, main = "T-Max")
