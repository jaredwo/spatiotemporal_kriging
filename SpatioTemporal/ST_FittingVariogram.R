# Part 2 of 4

# Compute sample variogram
tmin_var <- variogramST(tmin~1,data=obs_tmin_spdf,tunit="days",assumeRegular=T,na.omit=T, progress = TRUE, cutoff=3500) 
tmax_var <- variogramST(tmax~1,data=obs_tmax_spdf,tunit="days",assumeRegular=T,na.omit=T, progress = TRUE, cutoff=3500) 

#Inspect the sample variograms in two forms in order to find estimates for fitting the variogram

#Tmin
# Spatial: Nugget ~ .2 range ~ 3000km sill ~.7
# Temporal: Nugget ~ .42 range ~ 3 years sill ~ .55
# Joint: sill = .7
# Plot to view spatial component
plot(tmin_var,map=F) 
# 3D Plot
plot(tmin_var,wireframe=T)
# Plot to view temporal component
plot(tmin_var[tmin_var$spacelag==0,'timelag'],tmin_var[tmin_var$spacelag==0,'gamma'])

#Repeat for Tmax
# Spatial: Nugget ~ .15 range ~ 2000km sill ~.8
# Temporal: Nugget ~ .18 range ~ 3 years sill ~ .7
# Joint: sill = .9
plot(tmax_var,map=F)
plot(tmax_var,wireframe=T) 
plot(tmax_var[tmax_var$spacelag==0,'timelag'],tmax_var[tmax_var$spacelag==0,'gamma'])

#Estimate ST Anisotropy
tmin_stAni_day <- estiStAni(tmin_var, interval = c(10, 1000), method = "linear")
tmax_stAni_day <- estiStAni(tmax_var, interval = c(10, 1000), method = "linear")

# Fit using Sum Metric covariance model
#Tmin
#Using parameters from earlier model estimations
sumMetricModel_tmin <- vgmST("sumMetric",
                             space = vgm(psill=0, "Exp", range=1500, nugget=0),
                             time = vgm(psill=0,"Exp", range=1, nugget=0),
                             joint = vgm(psill=0.53, "Sph", range=3000, nugget=0.17 ),
                             stAni = tmin_stAni_day)

# Check the fit of the estimated model and make neccesary adjustments
plot(tmin_var, sumMetricModel_tmin,map=F, all=T)

# Try fitting using optim
# Use parscale to optimize the fitting of all parameters
fitSumMetric_tmin <- fit.StVariogram(tmin_var, sumMetricModel_tmin, fit.method=6, 
                                     method = "L-BFGS-B", stAni=tmin_stAni_day,
                                     control = list(parscale = c(1,10000,1,1,1,1,1,10000,1,10000), maxit=2e4))
# parscale: sill.s, range.s, nugget.s, sill.t ,  range.t ,   nugget.t , sill.st, range.st , nugget.st, ani

# Check the result of the fitting
plot(tmin_var, fitSumMetric_tmin, map=F, all=T)

# Save the best model
sumMetric_tmin <- fitSumMetric_tmin

#Repeat for tmax
sumMetricModel_tmax <- vgmST("sumMetric",
                             space = vgm(psill=0, "Exp", range=1500, nugget=0),
                             time = vgm(psill=0,"Exp", range=1, nugget=0),
                             joint = vgm(psill=0.68, "Sph", range=3000, nugget=0.16),
                             stAni = tmax_stAni_day)

# Check the estimated fit
plot(tmax_var, sumMetricModel_tmax, map=F, all=T)

# Try fitting with optim
fitSumMetric_tmax <- fit.StVariogram(tmax_var, sumMetricModel_tmax, fit.method=6, 
                                     method = "L-BFGS-B", stAni=tmin_stAni_day,
                                     control = list(parscale = c(1,10000,1,1,1,1,1,10000,1,10000), maxit=2e4))

#Check the results of the optim fit
plot(tmax_var, fitSumMetric_tmax,map=F, all=T)

#Save the best model
sumMetric_tmax <- fitSumMetric_tmax

# Look at both tmin and tmax results
# Plot empirical vs. fitted variograms
plot(tmin_var, sumMetric_tmin, map=F,main = "SumMetric Tmin", all=T)
plot(tmax_var, sumMetric_tmax, map=F,main = "SumMetric Tmax", all=T)

# Plot 3D differences
plot(tmin_var, sumMetric_tmin, wireframe=T,diff = T, main = "SumMetric Difference Tmin")
plot(tmax_var_day, sumMetric_tmax, wireframe=T,diff = T, main = "SumMetric Difference Tmax")

# MSE Measures
attr(sumMetric_tmin, "MSE")
attr(sumMetric_tmax, "MSE")