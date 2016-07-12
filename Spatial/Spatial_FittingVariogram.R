# Part 2 of 4

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

#Convert to gstatVariogram
emp_var_tmin <- stvario_to_svario(emp_varst_tmin)
emp_var_tmax <- stvario_to_svario(emp_varst_tmax)

# Can now use emp_var_tmin/tmax with any pure spatial gstat functions
#Estimate the tmin model to begin with
vgmModel_tmin <- vgm(psill = .7, model = "Exp", range = 2000, nugget = .15)

#Check how well this estimate fits. Make adjustments as neccesary
plot(emp_var_tmin, vgmModel_tmin, all=T)

#Try fitting using optim and check the result
fitModel_tmin <- fit.variogram(emp_var_tmin, vgmModel_tmin)
plot(emp_var_tmin, fitModel_tmin, all=T, main = "Tmin: Fit Pooled Spatial Variogram")

#Repeat for tmax
vgmModel_tmax <- vgm(psill = 1.1, model = "Exp", range = 2000, nugget = .11)

plot(emp_var_tmax, vgmModel_tmax, all=T)

fitModel_tmax <- fit.variogram(emp_var_tmax, vgmModel_tmax)
plot(emp_var_tmax, fitModel_tmax, all=T,  main = "Tmax: Fit Pooled Spatial Variogram")