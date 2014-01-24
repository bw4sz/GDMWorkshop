####################################################################################
## Generalized Dissimilarity modelling - Estimating turnover due to climate change ##

## 4 January 2014 - Dan Rosauer, Australian National University (dan.rosauer@anu.edu.au) ##

## This script requires:
## (1) transformed grids from a GDM model;
## (2) transformed grids from the same model, projected to a different climate;
## (3) the y intercept value from the model.

###################################################################################

library(raster)

# define the file locations
work.dir          <- "C:/GDM_workshop/out/"   # EDIT THIS TO YOUR LOCATION
output.filename   <- "climate_turnover.asc"

# define the input data.  The grids from the two time periods must be listed in the same order
tgrids_2010 <- c("isothermality_2010_1kmTran.flt","jan_rad_2010_1kmTran.flt","mintempjuly_2010_1kmTran.flt","precip_pet_ratio_2010_1kmTran.flt")
tgrids_2100 <- c("isothermality_2100_1kmTran.flt","jan_rad_2100_1kmTran.flt","mintempjuly_2100_1kmTran.flt","precip_pet_ratio_2100_1kmTran.flt")
model_intercept <- 0 ## use zero, though model intercept = 0.745343; if modeling through space and not time you might use the intercept and not set it to 0

setwd(work.dir)

# count the predictor grids
grid_count <- length(tgrids_2010)

# loop through the predictors, and for each predictor get the absolute difference in the transformed value for each pixel,
#   summed across all predictors
for (i in 1:grid_count) {
  current_trans.ras <- raster(tgrids_2010[i])
  future_trans.ras  <- raster(tgrids_2100[i])
  if (exists("ecol_dist.ras")) {
    ecol_dist.ras <- ecol_dist.ras + (abs(current_trans.ras - future_trans.ras))
  } else {
    ecol_dist.ras <- abs(current_trans.ras - future_trans.ras)
  }
}

# convert the ecological distance to a value for potential turnover
dist.ras = 1 - exp(-1 * (model_intercept + ecol_dist.ras))  

# plot and save the turnover result
plot(dist.ras)
writeRaster(dist.ras, output.filename, overwrite=TRUE)

