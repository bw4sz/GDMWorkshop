####################################################################################
## Generalized Dissimilarity modelling - Coloured map from an unsupervised classification ##

## 3 January 2014 - Dan Rosauer, Australian National University (dan.rosauer@anu.edu.au) ##

## This script requires:
## (1) transformed grids from a GDM model;
## (2) transformed grids from the same model, projected to a different climate;
## (3) the y intercept value from the model.

###################################################################################

library(raster)

# define the file locations
class.dir <- "C:/GDM_workshop/out/Class_01/"   # EDIT THIS TO YOUR LOCATION
input.raster.filename  <- "Class_01.flt"
input.colour.filename  <- "Class_01_RGB.csv"
output.image.name <-  "Classes_01_rgb"

setwd(class.dir)

# read the colours file
col <- read.csv(input.colour.filename)

# create 3 copies of the classes
classes.ras <- raster(input.raster.filename)
classes1.ras <- classes.ras
classes2.ras <- classes.ras
classes3.ras <- classes.ras

# apply colour
class.count <- nrow(col)
for (i in 1:class.count) {
  classes1.ras[classes.ras==i] <- col[i,1]
  classes2.ras[classes.ras==i] <- col[i,2]
  classes3.ras[classes.ras==i] <- col[i,3]
  cat(i,col[i,1],"\t",col[i,2],"\t",col[i,3],"\n")
}

rm(classes.ras)

# create a stack of the three layers
classes.stack <- stack(classes1.ras,classes2.ras,classes3.ras)

# plot the colour class map
windows()
plotRGB(classes.stack,r=1,g=2,b=3)

# show all of the different allocations of red, green and blue
windows(12,12)
par(mfrow=c(2,3))
plotRGB(classes.stack,r=1,g=2,b=3)
plotRGB(classes.stack,r=1,g=3,b=2)
plotRGB(classes.stack,r=2,g=1,b=3)
plotRGB(classes.stack,r=2,g=3,b=1)
plotRGB(classes.stack,r=3,g=1,b=2)
plotRGB(classes.stack,r=3,g=2,b=1)

# replot the colour class map comparing colours stretch options
windows()
par(mfrow=c(1,3))
plotRGB(classes.stack,r=1,g=2,b=3)
plotRGB(classes.stack,r=1,g=2,b=3, stretch="lin")
plotRGB(classes.stack,r=1,g=2,b=3, stretch="hist")

# show all of the different allocations of red, green and blue
windows(12,12)
par(mfrow=c(2,3))
plotRGB(classes.stack,r=1,g=2,b=3, stretch="hist")
plotRGB(classes.stack,r=1,g=3,b=2, stretch="hist")
plotRGB(classes.stack,r=2,g=1,b=3, stretch="hist")
plotRGB(classes.stack,r=2,g=3,b=1, stretch="hist")
plotRGB(classes.stack,r=3,g=1,b=2, stretch="hist")
plotRGB(classes.stack,r=3,g=2,b=1, stretch="hist")

# and a final preferred image to save
windows()
plotRGB(classes.stack,r=2,g=1,b=3, stretch="hist")
