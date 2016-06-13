# Example of loading precipitation station observations for a netcdf file

library(RNetCDF)
library(ncdf4)
library(sp)

# Open prcp netcdf file
ds <- nc_open('data/prcp_tobs_adj_19480101_20151231.nc')

# Print out header of netcdf dataset showing variables and dimensions
print(ds)

#Load dataframe of station metadata
stn_id <- ncvar_get(ds,varid='station_id')
stn_name <- ncvar_get(ds,varid='station_name')
stn_lon <- ncvar_get(ds,varid='longitude')
stn_lat <- ncvar_get(ds,varid='latitude')
stn_elev <- ncvar_get(ds,varid='elevation')
stns <- data.frame(stn_id=stn_id,stn_name=stn_name,stn_lon=stn_lon,
                   stn_lat=stn_lat,stn_elev=stn_elev)

#Turn station dataframe into spatial dataframe
coordinates(stns) <- ~stn_lon+stn_lat
proj4string(stns)=CRS("+proj=longlat +datum=WGS84")

# Get units for the time dimension
time_units <- ncatt_get(ds,'time','units')$value
# Convert time values to date objects
dates <- as.Date(utcal.nc(time_units, ncvar_get(ds,'time'), type="c"))

# Get indices for 2011
time_i <- which((dates >= "2011-01-01") & (dates <= "2011-12-31"))
start_i <- time_i[1]
count_i <- length(time_i)

# Read prcp observations for 2011
obs_prcp <- ncvar_get(ds,varid='prcp',start=c(1,start_i),count=c(-1,count_i))

# Transpose so that matrix is: time X station (i.e.--each column is a station time series)
obs_prcp <- t(obs_prcp)

# Get number of missing observations for each station
nmiss <- colSums(is.na(obs_prcp))

# We want stations that have at least one observation over this time period
# We can drop the rest
mask_keep <- nmiss < 365

# Subset stns SpatialPointsDataFrame and obs_prcp using the keep mask
stns <- stns[mask_keep,]
obs_prcp <- obs_prcp[,mask_keep]

# Convert to spacetime object...
# obs_prcp is a "space-wide" table

