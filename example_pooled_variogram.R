# Example of creating a "pooled" spatial variogram from a spatiotemporal dataset

library(spacetime)
library(sp)
library(gstat)

df_tmin <- read.csv('data/ann_anoms_tmin.csv',stringsAsFactors = FALSE)
df_tmin$year <- as.POSIXct(paste(df_tmin$year, "01", "01", sep = "-"), tz='GMT')
stidf_tmin <- stConstruct(x=df_tmin,c('longitude','latitude'),
                          time="year", crs = CRS("+proj=longlat +datum=WGS84"))
stfdf_tmin <- as(stidf_tmin, "STFDF")

# Calculate a "pooled" spatial only variogram by setting the tlags param = 0
# Set cutoff to 3500 km to make sure we get the range
emp_varst_tmin <- variogramST(tmin~1, data=stfdf_tmin, tlags=0,
                              progress=TRUE, assumeRegular=TRUE,
                              cutoff=3500)
# Plot pooled spatial sample variogram
plot(emp_varst_tmin, map=F)

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
# Can now use emp_var_tmin with any pure spatial gstat functions
