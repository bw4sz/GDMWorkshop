#####################################################################################################
## Generalized Dissimilarity modelling - Estimating habitat and species representation in reserves ##

## 3 January 2014 - Dan Rosauer, Australian National University (dan.rosauer@anu.edu.au)           ##

## This script estimates for each pixel:
##    the effective habitat area of similar habitat reserved as a proportion of the original habitat
##    the estimated proportion of the pixel's species that would be retained in the reserved areas

## The methods follow Allnutt et al. (2008) in Conservation Letters, pp173-181

## It could easily be adapted to estimate habitat and species loss under different land-use change
##    scenarios by replacing the reserve density grid.

## It requires:
## (1) an unmasked (pristine) density grid from a GDM model;
## (2) a density grid from the same GDM model, using a domain grid for habitat condition or reserve status (reserve).
## Both density grids must be produced using the sampling strategy (eg same number of training samples)

###################################################################################

library(raster)

# define the file locations
work.dir <- "C:/GDM_workshop/out"   # EDIT THIS TO YOUR LOCATION
pristine_density.filename  <- "density_presitine.flt"
reserve_density.filename   <- "density_reserve.flt"
output_prop.filename       <- "reserve_proportion.asc"
output_sar.filename        <- "reserve_sar.asc"

setwd(work.dir)

# load the density grids produced in GDM Modeller
pristine_dens.ras <- raster(pristine_density.filename)
reserve_dens.ras  <- raster(reserve_density.filename)

# calculate the effective habitat area ratio
output_prop.ras <- reserve_dens.ras / pristine_dens.ras

# calculate the species proportion using the species area curve function
output_sar.ras  <- output_prop.ras^0.25

# save the results
writeRaster(output_prop.ras,output_prop.filename,format="ascii",overwrite=TRUE)
writeRaster(output_sar.ras,output_sar.filename,format="ascii",overwrite=TRUE)

cat("\nHabitat loss estimates written to",output_prop.filename,"and",output_sar.filename,"\n")

# optional plot of the results
plot(output_sar.ras, main="Estimated proportion of species protected by reserves alone")

