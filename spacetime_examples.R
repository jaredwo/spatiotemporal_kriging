# Example of loading data into spacetime data structures

library(spacetime)
library(sp)

df_tmin <- read.csv('data/ann_anoms_tmin.csv',stringsAsFactors = FALSE)
df_tmax <- read.csv('data/ann_anoms_tmax.csv',stringsAsFactors = FALSE)

#coordinates(df_tmin) <- ~longitude+latitude
#proj4string(df_tmin)=CRS("+proj=longlat +datum=WGS84")
df_tmin$year <- as.Date(paste(df_tmin$year, "01", "01", sep = "-"))
df_tmax$year <- as.Date(paste(df_tmax$year, "01", "01", sep = "-"))

stidf_tmin <- stConstruct(x=df_tmin,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stidf_tmax <- stConstruct(x=df_tmax,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmin <- as(stidf_tmin, "STFDF")
stfdf_tmax <- as(stidf_tmax, "STFDF")

# Plot data for years 2010-2015
stplot(stfdf_tmax[,"2010/2015",'tmax'])