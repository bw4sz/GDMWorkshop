######################################################################################################
## Generalized Dissimilarity modelling - Generate Site Pairs                                        ##

## 3 January 2014 - Dan Rosauer, Australian National University (dan.rosauer@anu.edu.au)            ##

## This script generates a table with compositional dssimilarity between all pairs of sites         ##
##    along with their location, in the format required by GDM for Tables in R and by GDM Modeller  ##

## It requires:                                                                                     ##
## (1) a sites by species matrix, with a row for each site, and columns                             ##
##    for X, Y, then for each species.                                                              ##

######################################################################################################

library(ecodist)

# define the file locations
work.dir <- "Tas_Plant_GDM_Files/"   # EDIT THIS TO YOUR LOCATION
input.filename <- "Tas_Plant_Composition_5Dec13.csv"
output.filename <- "site_pairs_formatted.csv"

# load the sites
setwd(work.dir)
composition <- read.table(input.filename,header=T,sep=",")

# examine the composition data (just the first 20 columns)
head(composition[,1:25])

# remove the X, Y columns
composition_no_xy <- composition[,-(1:2)] 

# generate a Bray-Curtis (Sorenson) distance matrix for the first 10 sites and view it
dist.matrix <- bcdist(composition_no_xy[1:10,],rmzero=TRUE)
dist.matrix <- fixdmat(dist.matrix) # convert to regular matrix format
dist.matrix

# now generate the full matrix
dist.matrix <- bcdist(composition_no_xy,rmzero=TRUE)
dist.matrix <- fixdmat(dist.matrix) # convert to regular matrix format

# reformat in the 6 column format for GDM
site_count <- nrow(composition)

# Create an empty data frame
GDM_input <- data.frame(dist=0,weights=0,X0=0,Y0=0,X1=0,Y1=0,stringsAsFactors=FALSE)
m <- 1 # a row counter

cat("\nSite: ")

for (i in 1:(site_count-1)) {
  cat(i,"\t")
  
  #add each site pair
  for (j in (i+1):site_count) {
    dist  <- dist.matrix[i,j]
    XY0   <- composition[i,1:2]
    XY1   <- composition[j,1:2]
    row <- as.numeric(c(dist,1,XY0,XY1))
    GDM_input[m,] <- row
    m <- m+1
  }
}

rm(XY0,XY1,composition_no_xy,dist.matrix,dist,i,j,m,row)

write.csv(GDM_input,output.filename,row.names=FALSE)

cat("\nGDM response data created as",output.filename)