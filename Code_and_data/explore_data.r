###################################################################################################
## Generalized Dissimilarity modelling - Display the gridded environment data & samplE locations ##

## 2 January 2014 - Dan Rosauer, Australian National University (dan.rosauer@anu.edu.au)         ##

## This script maps each of the 2010 envrironment rasters and the sampled flora sites.           ##

###################################################################################################

library("raster")

# define the file locations
tas.dir <- "C:/GDM_workshop/Code_and_data/Tas_Plant_GDM_Files/"   # EDIT THIS TO YOUR LOCATION

grid.dir <- paste(tas.dir,"ASCII grids",sep="")

# change directories
setwd(grid.dir)

# load the environment rasters
jan_rad_2010 <- raster("jan_rad_2010_1km.asc")
mintempjuly_2010 <- raster("mintempjuly_2010_1km.asc")
precip_pet_ratio_2010 <- raster("precip_pet_ratio_2010_1km.asc")
isothermality_2010 <- raster("isothermality_2010_1km.asc")
bulkdensity <- raster("bulkdensity_1km.asc")
reserves <- raster("reserves1km.asc")
pristine <- raster("pristine_1mask_1km.asc")

# examine one of the rasters
mintempjuly_2010

# plot the environment rasters - all on the one graphics window
windows(12,12)
par(mfrow=c(3,3))
plot(jan_rad_2010,main="Jan rad 2010")
plot(mintempjuly_2010,main="min temp July 2010")
plot(precip_pet_ratio_2010,main="precip pet ratio 2010")
plot(isothermality_2010,main="isothermality 2010")
plot(bulkdensity,main="bulk density")
plot(reserves,main="reserves")
plot(pristine,main="pristine")

# load the sites
setwd(tas.dir)
composition.data<-read.table("Tas_Plant_Composition_5Dec13.csv",header=T,sep=",")
sites <- composition.data[,1:2]

# examine the composition data (just the first 20 columns)
head(composition.data[,1:25])

# plot the sites geographically
windows(12,12)
plot(jan_rad_2010, main="January radiation & vegetation sites")
points(sites, pch=20)

