Sys.setenv(PROJ_LIB=paste0(Sys.getenv("HOME"),"/miniconda2/envs/r/share/proj"))

library(maptools)
library(gstat)
library(rgdal)
library(sp)
library(spacetime)
library(raster)

setwd('/mizuna/s0/jwo118/summer_scholars/spacetime_example')

df_tmin <- read.csv('ann_anoms_tmin.csv',stringsAsFactors = FALSE)
df_tmax <- read.csv('ann_anoms_tmax.csv',stringsAsFactors = FALSE)

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

# Compute sample variogram
tmin_var <- variogramST(tmin~1,data=stfdf_tmin,tunit="years",assumeRegular=T,na.omit=T, progress = TRUE) 
tmax_var <- variogramST(tmax~1,data=stfdf_tmax,tunit="years",assumeRegular=T,na.omit=T, progress = TRUE) 

#Check it out
plot(tmin_var,wireframe=T) 
plot(tmax_var,wireframe=T) 

#Estimate ST Anisotropy
tmin_stAni <- estiStAni(tmin_var, interval = c(.01, 2), method = "linear")
tmax_stAni <- estiStAni(tmax_var, interval = c(.01, 2), method = "linear")

#Seperable Model
separableModel_tmin <- vgmST("separable",
                        space = vgm(.9, "Exp", 400, 0.1),
                        time = vgm(.95, "Sph", 3, 0.1),
                        sill = 124)

separableModel_tmax <- vgmST("separable",
                        space = vgm(0.9, "Exp", 400, 0.1),
                        time = vgm(1, "Sph", 3, 0.1),
                        sill = 124)
fitSeparable_tmin <- fit.StVariogram(tmin_var, separableModel_tmin, fit.method=7, method="BFGS",
                                     stAni=tmin_stAni)
fitSeparable_tmax <- fit.StVariogram(tmax_var, separableModel_tmax, fit.method=7, method="BFGS",
                                stAni=tmax_stAni)

#Product-Sum Model
prodSumModel <- vgmST("productSum",space = vgm(1, "Exp", 150, 0.5),
                      time = vgm(1, "Exp", 5, 0.5),k = 50) 
fitProdSum_tmin <- fit.StVariogram(tmin_var, prodSumModel,method = "BFGS")
fitProdSum_tmax <- fit.StVariogram(tmax_var, prodSumModel,method = "BFGS")

# Metric Model
metric_tmin <- vgmST("metric", joint = vgm(50,"Mat", 500, 0), stAni=tmin_stAni) 
metric_tmax <- vgmST("metric", joint = vgm(50,"Mat", 500, 0), stAni=tmax_stAni) 
fitMetric_tmin <- fit.StVariogram(tmin_var, metric_tmin, method="BFGS")
fitMetric_tmax <- fit.StVariogram(tmax_var, metric_tmax, method="BFGS")


# Sum Metric Model
sumMetricModel_tmin <- vgmST("sumMetric",space = vgm(20, "Sph", 150, 1),time = vgm(10, "Exp", 2, 0.5),joint = vgm(80, "Sph", 1500, 2.5),
                        stAni = tmin_stAni)
sumMetricModel_tmax <- vgmST("sumMetric",space = vgm(20, "Sph", 150, 1),time = vgm(10, "Exp", 2, 0.5),joint = vgm(80, "Sph", 1500, 2.5),
                             stAni = tmax_stAni)
fitSumMetric_tmin <- fit.StVariogram(tmin_var, sumMetricModel_tmin, fit.method = 7,
                                stAni = tmin_stAni, method = "BFGS",
                                control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                               maxit=1e4),
                                lower = c(sill.s = 0, range.s = 10, nugget.s = 0,
                                          sill.t = 0, range.t = 0.1, nugget.t = 0,
                                          sill.st= 0, range.st = 10, nugget.st = 0,
                                          anis = tmin_stAni-.1),
                                upper = c(sill.s = 200, range.s = 1E3, nugget.s = 20,
                                          sill.t = 200, range.t = 75, nugget.t = 20,
                                          sill.st= 200, range.st = 5E3, nugget.st = 20,
                                          anis = tmin_stAni+.1))

fitSumMetric_tmax <- fit.StVariogram(tmax_var, sumMetricModel_tmin, fit.method = 7,
                                     stAni = tmax_stAni, method = "BFGS",
                                     control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                                    maxit=1e4),
                                     lower = c(sill.s = 0, range.s = 10, nugget.s = 0,
                                               sill.t = 0, range.t = 0.1, nugget.t = 0,
                                               sill.st= 0, range.st = 10, nugget.st = 0,
                                               anis = tmax_stAni-.1),
                                     upper = c(sill.s = 200, range.s = 1E3, nugget.s = 20,
                                               sill.t = 200, range.t = 75, nugget.t = 20,
                                               sill.st= 200, range.st = 5E3, nugget.st = 20,
                                               anis = tmax_stAni+.1))

# Simple Sum Metric Model
simpleSumMetricModel_tmin <- vgmST("simpleSumMetric", space=vgm(120,"Sph", 150),
                                   time =vgm(120,"Exp", 10), 
                                   joint=vgm(120,"Sph", 150), 
                                   nugget=10, stAni=tmin_stAni)
simpleSumMetricModel_tmax <- vgmST("simpleSumMetric", space=vgm(120,"Sph", 150),
                                   time =vgm(120,"Exp", 10), 
                                   joint=vgm(120,"Sph", 150), 
                                   nugget=10, stAni=tmax_stAni)
#Not sure about upper/lower parameters
fitSimpleSumMetric_tmin <- fit.StVariogram(tmin_var, simpleSumMetricModel_tmin, fit.method = 7, stAni = tmin_stAni, method = "L-BFGS-B",
                                           control = list(parscale = c(1,10,1,1,1,100,1,10)),
                lower = c(sill.s = 0, range.s = 10, sill.t = 0, range.t = 0.1, sill.st= 0, range.st= 10, nugget=0, anis = tmin_stAni-.1), 
                upper = c(sill.s = 200, range.s = 500, sill.t = 200, range.t = 20, sill.st= 200, range.st = 5000, nugget = 100, anis = tmin_stAni+.1))

fitSimpleSumMetric_tmax <- fit.StVariogram(tmax_var, simpleSumMetricModel_tmax, fit.method = 7, stAni = tmax_stAni, method = "L-BFGS-B",
                                           control = list(parscale = c(1,10,1,1,1,100,1,10)),
                                           lower = c(sill.s = 0, range.s = 10, sill.t = 0, range.t = 0.1, sill.st= 0, range.st= 10, nugget=0, anis = tmax_stAni-.1), 
                                           upper = c(sill.s = 200, range.s = 500, sill.t = 200, range.t = 20, sill.st= 200, range.st = 5000, nugget = 100, anis = tmax_stAni+.1))

# Plot empirical vs. fitted variograms
par(mfrow=c(3,3))
plot(tmin_var, wireframe=T, main = "Sample Tmin")
plot(tmin_var, fitSeparable_tmin, wireframe=T,main = "Separable Tmin")
plot(tmin_var, fitProdSum_tmin, wireframe=T,main = "ProdSum Tmin")
plot(tmin_var, fitMetric_tmin,wireframe=T,main = "Metric Tmin")
plot(tmin_var, fitSumMetric_tmin, wireframe=T,main = "SumMetric Tmin")
plot(tmin_var, fitSimpleSumMetric_tmin, wireframe=T,main = "SimSumMetric Tmin")

plot(tmax_var, wireframe=T, main = "Sample Tmax")
plot(tmax_var, fitSeparable_tmax, wireframe=T,main = "Separable Tmax")
plot(tmax_var, fitProdSum_tmax, wireframe=T,main = "ProdSum Tmax")
plot(tmax_var, fitMetric_tmax,wireframe=T,main = "Metric Tmax")
plot(tmax_var, fitSumMetric_tmax, wireframe=T,main = "SumMetric Tmax")
plot(tmax_var, fitSimpleSumMetric_tmax, wireframe=T,main = "SimSumMetric Tmax")

# Differences
plot(tmin_var, wireframe=T, main = "Sample")
plot(tmin_var, fitSeparable_tmin, wireframe=T, diff = T, main = "Separable Difference Tmin")
plot(tmin_var, fitProdSum_tmin, wireframe=T,diff = T, main = "ProdSum Difference Tmin")
plot(tmin_var, fitMetric_tmin,wireframe=T,diff = T, main = "Metric Difference Tmin")
plot(tmin_var, fitSumMetric_tmin, wireframe=T,diff = T, main = "SumMetric Difference Tmin")
plot(tmin_var, fitSimpleSumMetric_tmin, wireframe=T,diff = T, main = "SimSumMetric TDifference min")

plot(tmax_var, wireframe=T, main = "Sample Tmax")
plot(tmax_var, fitSeparable_tmax, wireframe=T,diff = T, main = "Separable Difference Tmax")
plot(tmax_var, fitProdSum_tmax, wireframe=T,diff = T, main = "ProdSum Difference Tmax")
plot(tmax_var, fitMetric_tmax,wireframe=T,diff = T, main = "Metric Difference Tmax")
plot(tmax_var, fitSumMetric_tmax, wireframe=T,diff = T, main = "SumMetric Difference Tmax")
plot(tmax_var, fitSimpleSumMetric_tmax, wireframe=T,diff = T, main = "SimSumMetric Difference Tmax")


plot(tmax_var, wireframe=T, main = "Sample Tmax")
plot(tmax_var, fitSeparable_tmax, wireframe=T,diff = T, main = "Separable Difference Tmax")
plot(tmax_var, fitProdSum_tmax, wireframe=T,diff = T, main = "ProdSum Difference Tmax")
plot(tmax_var, fitMetric_tmax,wireframe=T,diff = T, main = "Metric Difference Tmax")
plot(tmax_var, fitSumMetric_tmax, wireframe=T,diff = T, main = "SumMetric Difference Tmax")
plot(tmax_var, fitSimpleSumMetric_tmax, wireframe=T,diff = T, main = "SimSumMetric Difference Tmax")

# MSE Measures
attr(fitSeparable_tmin, "MSE")
attr(fitProdSum_tmin, "MSE")
attr(fitMetric_tmin, "MSE")
attr(fitSumMetric_tmin, "MSE")
attr(fitSimpleSumMetric_tmin, "MSE")

attr(fitProdSum_tmax, "MSE")
attr(fitMetric_tmax, "MSE")
attr(fitSumMetric_tmax, "MSE")
attr(fitSimpleSumMetric_tmax, "MSE")


#Create prediction grid
tm.grid <- seq(as.POSIXct('2011-12-12 06:00 CET'),as.POSIXct('2011-12-14 09:00 CET'),
               length.out=5) 
grid.ST <- STF(sp.grid.UTM,tm.grid) 

#Kriging
precipPredict <- krigeST(PPB~1, data=timeDF, modelList=fitSumMetric, newdata=grid.ST) 
stplot(precipPredict)

#Cross Validation



# Bibliography
#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/
