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

#Inspect the sample variograms in two forms in order to find estimates for fitting the variogram
#Tmin
# Spatial: Nugget ~ .18 range ~ 2000km sill ~.8
# Temporal: Nugget ~ .19 range ~ 3 years sill ~ .5
# Joint: sill = .7
plot(tmin_var,map=F) 
plot(tmin_var,wireframe=T)
#Tmax
# Spatial: Nugget ~ .15 range ~ 2000km sill ~.8
# Temporal: Nugget ~ .18 range ~ 3 years sill ~ .7
# Joint: sill = .9
plot(tmax_var,map=F)
plot(tmax_var,wireframe=T) 

#Estimate ST Anisotropy
tmin_stAni <- estiStAni(tmin_var, interval = c(.01, 2), method = "linear")
tmax_stAni <- estiStAni(tmax_var, interval = c(.01, 2), method = "linear")

#Seperable Models
  #Tmin (try different vgm models too - Exp vs. Sph)
separableModel_tmin1 <- vgmST("separable",
                        space = vgm(.6, "Exp", 2000, 0.18),
                        time = vgm(.5, "Exp", 3, 0.19),
                        sill = .8)
separableModel_tmin2 <- vgmST("separable",
                             space = vgm(.6, "Sph", 2000, 0.18),
                             time = vgm(.5, "Exp", 3, 0.19),
                             sill = .75)
separableModel_tmin3 <- vgmST("separable",
                             space = vgm(.6, "Exp", 2000, 0.18),
                             time = vgm(.5, "Sph", 3, 0.19),
                             sill = .8)
separableModel_tmin4 <- vgmST("separable",
                             space = vgm(.6, "Sph", 2000, 0.18),
                             time = vgm(.5, "Sph", 3, 0.19),
                             sill = .75)
    #Does this fit?
plot(tmin_var, separableModel_tmin1,map=F, all=T, main = "Exp Exp")
plot(tmin_var, separableModel_tmin2,map=F, all=T, main = "Sph Exp")   #2 Choice
plot(tmin_var, separableModel_tmin3,map=F, all=T, main = "Exp Sph")
plot(tmin_var, separableModel_tmin4,map=F, all=T, main = "Sph Sph")   #1 Choice
    #Try fitting best 2 (how about fit method 8?)
fitSeparable_tmin4 <- fit.StVariogram(tmin_var, separableModel_tmin4, fit.method=6, method="BFGS",
                                     stAni=tmin_stAni)
fitSeparable_tmin2 <- fit.StVariogram(tmin_var, separableModel_tmin2, fit.method=6, method="BFGS",
                                     stAni=tmin_stAni)
    #Does this fit better?
plot(tmin_var, fitSeparable_tmin4,map=F, all=T, main = "Sph Sph Fitted") # No
plot(tmin_var, fitSeparable_tmin2,map=F, all=T, main = "Sph Exp Fitted") # No
#Best
separableModel_tmin <- separableModel_tmin4
  #Tmax (try different vgm models too - Exp vs. Sph)
separableModel_tmax1 <- vgmST("separable",
                        space = vgm(0.8, "Exp", 2000, 0.17),
                        time = vgm(.6, "Exp", 3, 0.18),
                        sill = 1)
separableModel_tmax2 <- vgmST("separable",
                             space = vgm(0.75, "Sph", 2000, 0.17),
                             time = vgm(.7, "Exp", 3, 0.18),
                             sill = .9)
separableModel_tmax3 <- vgmST("separable",
                             space = vgm(0.8, "Exp", 2000, 0.17),
                             time = vgm(.7, "Sph", 3, 0.18),
                             sill = 1)
separableModel_tmax4 <- vgmST("separable",
                             space = vgm(0.75, "Sph", 2000, 0.17),
                             time = vgm(.6, "Sph", 3, 0.1),
                             sill = .9)
    #Does this fit?
plot(tmax_var, separableModel_tmax1,map=F, all=T, main = "Exp Exp")
plot(tmax_var, separableModel_tmax2,map=F, all=T, main = "Sph Exp") #2 Choice
plot(tmax_var, separableModel_tmax3,map=F, all=T, main = "Exp Sph")
plot(tmax_var, separableModel_tmax4,map=F, all=T, main = "Sph Sph") #1 Choice
#Try fitting best 2 (how about fit method 8? or others)
fitSeparable_tmax2 <- fit.StVariogram(tmax_var, separableModel_tmax2, fit.method=6, method="BFGS",
                                     stAni=tmax_stAni)
fitSeparable_tmax4 <- fit.StVariogram(tmax_var, separableModel_tmax4, fit.method=6, method="BFGS",
                                      stAni=tmax_stAni)
    #Does this fit better?
plot(tmax_var, fitSeparable_tmax2,map=F, all=T, main = "Sph Exp Fitted") #No
plot(tmax_var, fitSeparable_tmax4,map=F, all=T, main = "Sph Sph Fitted") #No
#Best
separableModel_tmax <- separableModel_tmax4
#Product-Sum Model
  #Tmin
prodSumModel_tmin1 <- vgmST("productSum", 
                      space = vgm(.6, "Exp", 2000, 0.18),
                      time = vgm(.2, "Exp", 3, 0.19),
                      k = .1) 
prodSumModel_tmin2 <- vgmST("productSum", 
                            space = vgm(.4, "Sph", 2000, 0.18),
                            time = vgm(.1, "Exp", 3, 0.19),
                            k = .1) 
prodSumModel_tmin3 <- vgmST("productSum", 
                            space = vgm(.65, "Exp", 2000, 0.18),
                            time = vgm(.01, "Sph", 3, 0.3),
                            k = .05) 
prodSumModel_tmin4 <- vgmST("productSum", 
                            space = vgm(.4, "Sph", 2000, 0.18),
                            time = vgm(.01, "Sph", 3, 0.3),
                            k = .05) 
#Does this fit?
plot(tmin_var, prodSumModel_tmin1,map=F, all=T, main = "Exp Exp")
plot(tmin_var, prodSumModel_tmin2,map=F, all=T, main = "Sph Exp")
plot(tmin_var, prodSumModel_tmin3,map=F, all=T, main = "Exp Sph") #2 Choice
plot(tmin_var, prodSumModel_tmin4,map=F, all=T, main = "Sph Sph") #1 Choice
#Try fitting best 2 (how about fit method 8?)
fitProdSum_tmin4 <- fit.StVariogram(tmin_var, prodSumModel_tmin4,fit.method=6,method = "BFGS")
fitProdSum_tmin3 <- fit.StVariogram(tmin_var, prodSumModel_tmin3,fit.method=6,method = "BFGS")
#Does this fit better?
plot(tmin_var, fitProdSum_tmin4,map=F, all=T, main = "Sph Sph Fitted") #Yes
plot(tmin_var, fitProdSum_tmin3,map=F, all=T, main = "Exp Sph Fitted") #Yes
#Best
prodSum_tmin <- fitProdSum_tmin4
  #Tmax
prodSumModel_tmax1 <- vgmST("productSum", 
                            space = vgm(0.9, "Exp", 2000, 0.17),
                            time = vgm(.25, "Exp", 3, 0.1),
                            k = .05) 
prodSumModel_tmax2 <- vgmST("productSum", 
                            space = vgm(0.75, "Sph", 2000, 0.17),
                            time = vgm(.1, "Exp", 3, 0.18),
                            k = .05)  
prodSumModel_tmax3 <- vgmST("productSum", 
                            space = vgm(0.9, "Exp", 2000, 0.17),
                            time = vgm(.25, "Sph", 3, 0.25),
                            k = .05)   
prodSumModel_tmax4 <- vgmST("productSum", 
                            space = vgm(0.6, "Sph", 2000, 0.17),
                            time = vgm(.15, "Sph", 3, .3),
                            k = .05)  
#Does this fit?
plot(tmax_var, prodSumModel_tmax1,map=F, all=T, main = "Exp Exp")
plot(tmax_var, prodSumModel_tmax2,map=F, all=T, main = "Sph Exp")
plot(tmax_var, prodSumModel_tmax3,map=F, all=T, main = "Exp Sph") #1 Choice
plot(tmax_var, prodSumModel_tmax4,map=F, all=T, main = "Sph Sph") #2 Choice
#Try fitting best 2 (how about fit method 8?)
fitProdSum_tmax3 <- fit.StVariogram(tmax_var, prodSumModel_tmax3,fit.method=6,method = "BFGS")
fitProdSum_tmax4 <- fit.StVariogram(tmax_var, prodSumModel_tmax4,fit.method=6,method = "BFGS")
#Does this fit better?
plot(tmax_var, fitProdSum_tmax3,map=F, all=T, main = "Sph Sph Fitted") #Yes
plot(tmax_var, fitProdSum_tmax4,map=F, all=T, main = "Exp Sph Fitted") #Yes
#Best
prodSum_tmax <- fitProdSum_tmax3

# Metric Model
  #Tmin
metric_tmin1 <- vgmST("metric", joint = vgm(.6,"Exp", 1500, .18), stAni=10) # Est stAni seems to be way too low
metric_tmin2 <- vgmST("metric", joint = vgm(.45,"Sph", 2000, .16), stAni=10)
metric_tmin3 <- vgmST("metric", joint = vgm(.6,"Mat", 1500, .16), stAni=10)
#Does this fit?
plot(tmin_var, metric_tmin1,map=F, all=T, main = "Exp") #2 Choice
plot(tmin_var, metric_tmin2,map=F, all=T, main = "Sph")
plot(tmin_var, metric_tmin3,map=F, all=T, main = "Mat") #1 Choice
#Try fitting best 2 (how about fit method 8?)
fitMetric_tmin3 <- fit.StVariogram(tmin_var,metric_tmin3,fit.method=6,method = "BFGS")
fitMetric_tmin1 <- fit.StVariogram(tmin_var,metric_tmin1,fit.method=6,method = "BFGS")
#Does this fit better?
plot(tmin_var, fitMetric_tmin3,map=F, all=T, main = "Mat Fitted") #Maybe?
plot(tmin_var, fitMetric_tmin1,map=F, all=T, main = "Exp Fitted") #Maybe?
#Best?
metric_tmin <- fitMetric_tmin3
  #Tmax
metric_tmax1 <- vgmST("metric", joint = vgm(.9,"Exp", 1500, .14), stAni=10)
metric_tmax2 <- vgmST("metric", joint = vgm(.65,"Sph", 2000, .16), stAni=7)
metric_tmax3 <- vgmST("metric", joint = vgm(.85,"Mat", 1500, .14), stAni=10)
#Does this fit?
plot(tmax_var, metric_tmax1,map=F, all=T, main = "Exp") #2 Choice
plot(tmax_var, metric_tmax2,map=F, all=T, main = "Sph")
plot(tmax_var, metric_tmax3,map=F, all=T, main = "Mat") #1 Choice
#Try fitting best 2 (how about fit method 8?)
fitMetric_tmax3 <- fit.StVariogram(tmax_var,metric_tmax3,fit.method=6,method = "BFGS")
fitMetric_tmax1 <- fit.StVariogram(tmax_var,metric_tmax1,fit.method=6,method = "BFGS")
#Does this fit better?
plot(tmax_var, fitMetric_tmax3,map=F, all=T, main = "Mat Fitted") #No
plot(tmax_var, fitMetric_tmax1,map=F, all=T, main = "Exp Fitted") #No
#Best
metric_tmax <- metric_tmax3

# Sum Metric Model
  #Tmin
#Using parameters from earlier model estimations
sumMetricModel_tmin1 <- vgmST("sumMetric",
                             space = vgm(.3, "Sph", 2000, 0.18),
                             time = vgm(.01, "Sph", 3, 0.3),
                             joint = vgm(.25,"Mat", 1500, 0),
                             stAni = 10)
# Does this fit?
plot(tmin_var, sumMetricModel_tmin,map=F, all=T)
# Try Fitting1
fitSumMetric_tmin1 <- fit.StVariogram(tmin_var, sumMetricModel_tmin1, fit.method = 6,
                                     stAni = 10, method = "BFGS")
# Does this fit better?
plot(tmin_var, fitSumMetric_tmin1,map=F, all=T) #YES
#Best
sumMetric_tmin <- fitSumMetric_tmin1
  #Tmax
sumMetricModel_tmax1 <- vgmST("sumMetric",
                             space = vgm(0.8, "Exp", 2000, 0.17),
                             time = vgm(.25, "Sph", 3, 0.25),
                             joint = vgm(.3,"Mat", 1500, 0),
                             stAni = 10)
#Does this fit?
plot(tmax_var, sumMetricModel_tmax1,map=F, all=T)
#Try Fitting
fitSumMetric_tmax1 <- fit.StVariogram(tmax_var, sumMetricModel_tmax1, fit.method = 6,
                                     stAni = 10, method = "BFGS")
#Does this fit better?
plot(tmax_var, fitSumMetric_tmax1,map=F, all=T) #YES 
#Best
sumMetric_tmax <- fitSumMetric_tmax1
# Simple Sum Metric Model
  #Tmin (using parameters from sum metric model to start)
simpleSumMetricModel_tmin1 <- vgmST("simpleSumMetric", 
                                   space=vgm(.12,"Sph",2000,.12),
                                   time =vgm(0,"Sph", 3,.03), 
                                   joint=vgm(.47,"Mat", 1500,.001), 
                                   nugget=0, stAni=10)
#Does this fit?
plot(tmin_var, simpleSumMetricModel_tmin1,map=F, all=T)
#Try Fitting
fitSimpleSumMetric_tmin1 <- fit.StVariogram(tmin_var, simpleSumMetricModel_tmin1, fit.method = 6,
                                     stAni = 10, method = "BFGS")
#Does this fit better?
plot(tmin_var, fitSimpleSumMetric_tmin1,map=F, all=T) #No
#Best
simpleSumMetric_tmin <- simpleSumMetricModel_tmin1
  #Tmax (using parameters from sum metric model to start)
simpleSumMetricModel_tmax1 <- vgmST("simpleSumMetric", 
                                   space=vgm(.35,"Exp",2000,.08),
                                   time =vgm(.01,"Sph", 3,.01), 
                                   joint=vgm(.66,"Mat", 1500,.003), 
                                   nugget=0, stAni=10)
#Does this fit?
plot(tmax_var, simpleSumMetricModel_tmax1,map=F, all=T)
#Try Fitting
fitSimpleSumMetric_tmax1 <- fit.StVariogram(tmax_var, simpleSumMetricModel_tmax1, fit.method = 6,
                                           stAni = 10, method = "BFGS")
#Does this fit better?
plot(tmax_var, fitSimpleSumMetric_tmax1,map=F, all=T) #No
#Best
simpleSumMetric_tmax <- simpleSumMetricModel_tmax1
#Not sure about upper/lower parameters
# fitSimpleSumMetric_tmin <- fit.StVariogram(tmin_var, simpleSumMetricModel_tmin, fit.method = 7, stAni = tmin_stAni, method = "L-BFGS-B",
#                                            control = list(parscale = c(1,10,1,1,1,100,1,10)),
#                 lower = c(sill.s = 0, range.s = 10, sill.t = 0, range.t = 0.1, sill.st= 0, range.st= 10, nugget=0, anis = tmin_stAni-.1), 
#                 upper = c(sill.s = 200, range.s = 500, sill.t = 200, range.t = 20, sill.st= 200, range.st = 5000, nugget = 100, anis = tmin_stAni+.1))
# 
# fitSimpleSumMetric_tmax <- fit.StVariogram(tmax_var, simpleSumMetricModel_tmax, fit.method = 7, stAni = tmax_stAni, method = "L-BFGS-B",
#                                            control = list(parscale = c(1,10,1,1,1,100,1,10)),
#                                            lower = c(sill.s = 0, range.s = 10, sill.t = 0, range.t = 0.1, sill.st= 0, range.st= 10, nugget=0, anis = tmax_stAni-.1), 
#                                            upper = c(sill.s = 200, range.s = 500, sill.t = 200, range.t = 20, sill.st= 200, range.st = 5000, nugget = 100, anis = tmax_stAni+.1))

# Plot empirical vs. fitted variograms
plot(tmin_var, separableModel_tmin, map=F,main = "Separable Tmin", all=T)
plot(tmin_var, prodSum_tmin, map=F,main = "ProdSum Tmin", all=T)
plot(tmin_var, metric_tmin,map=F,main = "Metric Tmin", all=T)
plot(tmin_var, sumMetric_tmin, map=F,main = "SumMetric Tmin", all=T)
plot(tmin_var, simpleSumMetric_tmin, map=F,main = "SimSumMetric Tmin", all=T)

plot(tmax_var, separableModel_tmax, map=F,main = "Separable Tmax", all=T)
plot(tmax_var, prodSum_tmax, map=F,main = "ProdSum Tmax", all=T)
plot(tmax_var, metric_tmax,map=F,main = "Metric Tmax", all=T)
plot(tmax_var, sumMetric_tmax, map=F,main = "SumMetric Tmax", all=T)
plot(tmax_var, simpleSumMetric_tmax, map=F,main = "SimSumMetric Tmax", all=T)

# Differences
plot(tmin_var, separableModel_tmin, wireframe=T, diff = T, main = "Separable Difference Tmin")
plot(tmin_var, prodSum_tmin, wireframe=T,diff = T, main = "ProdSum Difference Tmin")
plot(tmin_var, metric_tmin,wireframe=T,diff = T, main = "Metric Difference Tmin")
plot(tmin_var, sumMetric_tmin, wireframe=T,diff = T, main = "SumMetric Difference Tmin")
plot(tmin_var, simpleSumMetric_tmin, wireframe=T,diff = T, main = "SimSumMetric TDifference min")

plot(tmax_var, separableModel_tmax, wireframe=T,diff = T, main = "Separable Difference Tmax")
plot(tmax_var, prodSum_tmax, wireframe=T,diff = T, main = "ProdSum Difference Tmax")
plot(tmax_var, metric_tmax,wireframe=T,diff = T, main = "Metric Difference Tmax")
plot(tmax_var, sumMetric_tmax, wireframe=T,diff = T, main = "SumMetric Difference Tmax")
plot(tmax_var, simpleSumMetric_tmax, wireframe=T,diff = T, main = "SimSumMetric Difference Tmax")

# MSE Measures
attr(separableModel_tmin, "MSE")
attr(prodSum_tmin, "MSE")
attr(metric_tmin, "MSE")
attr(sumMetric_tmin, "MSE")
attr(simpleSumMetric_tmin, "MSE")

attr(separableModel_tmax, "MSE")
attr(prodSum_tmax, "MSE")
attr(metric_tmax, "MSE")
attr(sumMetric_tmax, "MSE")
attr(simpleSumMetric_tmax, "MSE")


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

separable <- vgmST("separable", space = vgm(-60,"Sph", 500, 1),time = vgm(35,"Sph", 500, 1), sill=0.56)
plot(separable, map=F)
