
####################################################################################
## Generalized Dissimilarity modelling - Significance test for models & variables ##

## Version 3 - 16 December 2013 - Karel Mokany, CSIRO Ecosystem Sciences (Karel.Mokany@csiro.au)##

## The following code is broken into three sections:
## (1) the first section demonstrates how to prepare the data & run a single GDM;
## (2) the second section is the significance test function;
## (3) the third section explains and demonstrates the application of the 
## significance test function.


###################################################################################
## (1) DATA PREPARATION & SIMPLE MODEL FITTING ##
## Now load the base GDM functions and your data into R ## 

## The following loads the GDM functions
## Users are to specify their own path to the relevant file.
source(paste(gitpath,"GDM_workspace\\Rfunctions\\GDM_Table_Funcs64_1_0.R",sep=""))
mydllpath<-"GDM_workspace\\Rfunctions\\GDMBinLib64.dll"

setwd("Tas_Plant_GDM_Files/")

## The following loads the composition data (sites = rows, species = columns) as "composition.data"
## and the independent variables data (sites = rows, independent variables = columns) as "ind.variables"
## Users are to specify their own path to the relevant files.

## upload the composition data
composition.data <- read.table("Tas_Plant_Composition_5Dec13.csv",header=T,sep=",")
composition.data<-as.matrix(composition.data[,-c(1,2)])

## upload the environment data
ind.variables<-read.table("Tas_Env_Predictors_5Dec13.csv",header=T,sep=",")
ind.variables<-as.matrix(ind.variables)

## convert the composition data to a dissimilarity matrix
dissimilarity.data<-matrix(NA,nrow(composition.data),nrow(composition.data))
for(i_site in 1:nrow(composition.data))
{
  richness_i_site<-0
  for(i_spp in 1:ncol(composition.data))
  {
    richness_i_site<-richness_i_site + composition.data[i_site,i_spp] 
  } ## end for i_spp
  for(j_site in 1:nrow(composition.data))
  {
    richness_j_site<-0    
    sharedspecies<-0
    for(i_spp in 1:ncol(composition.data))
    {
      richness_j_site<-richness_j_site + composition.data[j_site,i_spp]
      if(composition.data[i_site,i_spp]*composition.data[j_site,i_spp] == 1)
      {      
        sharedspecies<-sharedspecies + 1
      } ## end if  
    } ## end for i_spp    
    dissimilarity.data[i_site,j_site]<-(1 - ((2 * sharedspecies)/(richness_i_site + richness_j_site)))
  } ## end for j_site
} ## end for i_site


## Create a J-Table for application in GDM ##
nrows_JTable<-((nrow(dissimilarity.data)^2)-nrow(dissimilarity.data))/2
ncols_JTable<-2+(ncol(ind.variables)*2)
nEvar<-ncol(ind.variables)-2
JTable<-matrix(NA,nrows_JTable,ncols_JTable)
colnames(JTable)<-c("Dissimilarity","weight",colnames(ind.variables[,c(1,2)]),colnames(ind.variables[,c(1,2)]),colnames(ind.variables[,-c(1,2)]),colnames(ind.variables[,-c(1,2)]))
pair<-1
for(i_site in 1:(ncol(dissimilarity.data)-1))
{
  for(j_site in (i_site+1):nrow(dissimilarity.data))
  {
    JTable[pair,1]<-dissimilarity.data[i_site,j_site]
    JTable[pair,2]<-1
    JTable[pair,3]<-ind.variables[i_site,1]
    JTable[pair,4]<-ind.variables[i_site,2]
    JTable[pair,5]<-ind.variables[j_site,1]
    JTable[pair,6]<-ind.variables[j_site,2]
    for(var in 1:nEvar)
    {
      JTable[pair,(6+var)]<-ind.variables[i_site,(2+var)]
      JTable[pair,(6+nEvar+var)]<-ind.variables[j_site,(2+var)]
    } ## end for nvar
    pair<-pair+1 
  } ## end for j_site
} ## end for i_site

## Write the JTable to file
write.table(JTable, file = "TasPlant_JTable.csv", append = FALSE, quote = FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)

## convert matrix to data frame
JTable.F<-as.data.frame(JTable)
## Run GDM on the JTable
mod.1<-gdm.fitfromtable(JTable.F, geo=FALSE)
## show model summary
gdm.summaryfromtable(mod.1)
## plot model

windows()
gdm.plotfromtable(mod.1)


####################################################################################
## (2) THE SIGNIFICANCE TEST FUNCTION ##
## load in the following function (i.e. copy and paste into R) ##
## START FUNCTION
gdm.model.variable.test<-function(ind.variables, dissimilarity.data, geo=TRUE, splines=NULL, quantiles=NULL, single.model=TRUE ,permutations=10)
{
  nsim<-permutations
  geo.still.in<-geo
  n.variables<-if(geo.still.in==TRUE) length(ind.variables[1,])-1 else length(ind.variables[1,])-2
  var.names<-if(geo.still.in==TRUE) c("geo",colnames(ind.variables[,-c(1,2)])) else colnames(ind.variables[,-c(1,2)])
  Model.test.values<-matrix(NA,3,n.variables,dimnames = list(c("Model deviance", "Percent deviance explained", "Model p-value"),c(paste("Model.", 1:n.variables, sep = ""))))
  Deviance.reduction.variables<-matrix(NA,n.variables,((n.variables-1)*2),dimnames = list(c(1:n.variables),c(rep(c("variable","deviance reduction"),(n.variables-1)))))
  Pvalues.variables<-matrix(NA,n.variables,((n.variables-1)*2),dimnames = list(c(1:n.variables),c(rep(c("variable","p-value"),(n.variables-1)))))
  manip.ind.variables<-ind.variables
  permuted.deviance.list<-vector(mode="numeric",length=nsim)
  
  ## create the J-Table
  nrows_JTable<-((nrow(dissimilarity.data)^2)-nrow(dissimilarity.data))/2
  ncols_JTable<-2+(ncol(ind.variables)*2)
  nEvar<-ncol(ind.variables)-2
  nEvar.still.in<-nEvar
  n.var.still.in<-if(geo.still.in==TRUE) (nEvar.still.in + 1) else nEvar.still.in
  JTable<-matrix(NA,nrows_JTable,ncols_JTable)
  colnames(JTable)<-c("Dissimilarity","weight",colnames(ind.variables[,c(1,2)]),colnames(ind.variables[,c(1,2)]),colnames(ind.variables[,-c(1,2)]),colnames(ind.variables[,-c(1,2)]))
  pair<-1
  for(i_site in 1:(ncol(dissimilarity.data)-1))
  {
    for(j_site in (i_site+1):nrow(dissimilarity.data))
    {
      JTable[pair,1]<-dissimilarity.data[i_site,j_site]
      JTable[pair,2]<-1
      JTable[pair,3]<-ind.variables[i_site,1]
      JTable[pair,4]<-ind.variables[i_site,2]
      JTable[pair,5]<-ind.variables[j_site,1]
      JTable[pair,6]<-ind.variables[j_site,2]
      for(var in 1:nEvar)
      {
        JTable[pair,(6+var)]<-ind.variables[i_site,(2+var)]
        JTable[pair,(6+nEvar+var)]<-ind.variables[j_site,(2+var)]
      } ## end for nvar
      pair<-pair+1 
    } ## end for j_site
  } ## end for i_site
  
  ## the big loop over all v variables
  for(v in 1:n.variables)  
  {
    if(v==n.variables)
    {
      n.var.still.in<-1
    } ## end if v==n.variables
    else
    {
      var.names<-if(geo.still.in==TRUE) c("geo",colnames(manip.ind.variables[,-c(1,2)])) else colnames(manip.ind.variables[,-c(1,2)])
      Deviance.reduction.variables[c(1:n.var.still.in),((v*2)-1)]<-var.names
      Pvalues.variables[c(1:n.var.still.in),((v*2)-1)]<-var.names
    } ## end else v==n.variables
    ## run gdm including all variables
    JTable.F<-as.data.frame(JTable)
    all.variables.gdm<-gdm.fitfromtable(JTable.F, geo=geo.still.in, splines=NULL, quantiles=NULL)
    ## Test model significance by permuting all data together
    for(i in 1:nsim)
    {
      perm.ind.variables<-manip.ind.variables[sample(length(manip.ind.variables[,1]),length(manip.ind.variables[,1]),),] ## permute the rows of the independant variables
      ## Create permuted J-Table
      perm.JTable<-JTable
      pair<-1
      for(i_site in 1:(ncol(dissimilarity.data)-1))
      {
        for(j_site in (i_site+1):nrow(dissimilarity.data))
        {
          perm.JTable[pair,3]<-perm.ind.variables[i_site,1]
          perm.JTable[pair,4]<-perm.ind.variables[i_site,2]
          perm.JTable[pair,5]<-perm.ind.variables[j_site,1]
          perm.JTable[pair,6]<-perm.ind.variables[j_site,2]
          for(var in 1:nEvar.still.in)
          {
            perm.JTable[pair,(6+var)]<-perm.ind.variables[i_site,(2+var)]
            perm.JTable[pair,(6+nEvar.still.in+var)]<-perm.ind.variables[j_site,(2+var)]
          } ## end for nvar
          pair<-pair+1 
        } ## end for j_site
      } ## end for i_site
      ## run gdm using permuted variables
      JTable.F<-as.data.frame(perm.JTable)
      perm.variables.gdm<-gdm.fitfromtable(JTable.F, geo=geo.still.in, splines=NULL, quantiles=NULL)  
      permuted.deviance.list[i]<-perm.variables.gdm$gdmdeviance
    } ## end for i
    Model.test.values[1,v]<-all.variables.gdm$gdmdeviance
    Model.test.values[2,v]<-all.variables.gdm$explained
    Model.test.values[3,v]<-sum(permuted.deviance.list<=all.variables.gdm$gdmdeviance)/nsim
    if(n.var.still.in<2)
    {
      break
    } ## end if n.var.still.in<2 
    ## Now test each variable. Run this if geographic distance is included
    if(geo.still.in==TRUE) 
    {
      ## run gdm without geographic distance
      JTable.F<-as.data.frame(JTable)
      gdm.without.geo<-gdm.fitfromtable(JTable.F, geo=FALSE, splines=NULL, quantiles=NULL)
      permuted.deviance.reduction.list<-vector(mode="numeric",length=nsim)
      perm.ind.variables<-manip.ind.variables
      for(k in 1:nsim)
      {
        perm.ind.variables[,1]<-sample(perm.ind.variables[,1],length(perm.ind.variables[,1]))
        perm.ind.variables[,2]<-sample(perm.ind.variables[,2],length(perm.ind.variables[,2])) 
        perm.JTable<-JTable
        pair<-1
        for(i_site in 1:(ncol(dissimilarity.data)-1))
        {
          for(j_site in (i_site+1):nrow(dissimilarity.data))
          {
            perm.JTable[pair,3]<-perm.ind.variables[i_site,1]
            perm.JTable[pair,4]<-perm.ind.variables[i_site,2]
            perm.JTable[pair,5]<-perm.ind.variables[j_site,1]
            perm.JTable[pair,6]<-perm.ind.variables[j_site,2]
            pair<-pair+1
          } ## end for j_site
        } ## end for i_site        
        JTable.F<-as.data.frame(perm.JTable)
        perm.j.gdm<-gdm.fitfromtable(JTable.F,geo=geo.still.in, splines=NULL, quantiles=NULL) 
        permuted.deviance.reduction.list[k]<-gdm.without.geo$gdmdeviance-perm.j.gdm$gdmdeviance
      } ## end for k
      Deviance.reduction.variables[1,(v*2)]<-gdm.without.geo$gdmdeviance-all.variables.gdm$gdmdeviance
      Pvalues.variables[1,(v*2)]<-sum(permuted.deviance.reduction.list>=(gdm.without.geo$gdmdeviance-all.variables.gdm$gdmdeviance))/nsim
      for(j in 1:nEvar.still.in) ## now test all other variables
      {
        ## run gdm without variable j        
        perm.JTable<-JTable[,-c((6+j),(6+j+nEvar.still.in))]
        JTable.F<-as.data.frame(perm.JTable)
        gdm.without.j<-gdm.fitfromtable(JTable.F,geo=geo.still.in, splines=NULL, quantiles=NULL)
        permuted.deviance.reduction.list<-vector(mode="numeric",length=nsim)
        perm.JTable<-JTable
        perm.ind.variables<-manip.ind.variables
        for (k in 1:nsim)
        {
          perm.ind.variables[,(j+2)]<-sample(perm.ind.variables[,(j+2)],length(perm.ind.variables[,(j+2)]))
          pair<-1
          for(i_site in 1:(ncol(dissimilarity.data)-1))
          {
            for(j_site in (i_site+1):nrow(dissimilarity.data))
            {
              perm.JTable[pair,(6+j)]<-perm.ind.variables[i_site,(2+j)]
              perm.JTable[pair,(6+nEvar.still.in+j)]<-perm.ind.variables[j_site,(2+j)]
              pair<-pair+1 
            } ## end for j_site
          } ## end for i_site
          JTable.F<-as.data.frame(perm.JTable)
          perm.j.gdm<-gdm.fitfromtable(JTable.F,geo=geo.still.in, splines=NULL, quantiles=NULL) 
          permuted.deviance.reduction.list[k]<-gdm.without.j$gdmdeviance-perm.j.gdm$gdmdeviance
        } ## end for k 
        Deviance.reduction.variables[(j+1),(v*2)]<-gdm.without.j$gdmdeviance-all.variables.gdm$gdmdeviance
        Pvalues.variables[(j+1),(v*2)]<-sum(permuted.deviance.reduction.list>=(gdm.without.j$gdmdeviance-all.variables.gdm$gdmdeviance))/nsim
      } ## end for j loop 
    } ## end if geo.still.in==TRUE
    else ## run this if geographic distance is not included
    { 
      for(j in 1:nEvar.still.in)
      { 
        perm.JTable<-JTable[,-c((6+j),(6+j+nEvar.still.in))]
        JTable.F<-as.data.frame(perm.JTable)
        gdm.without.j<-gdm.fitfromtable(JTable.F,geo=geo.still.in, splines=NULL, quantiles=NULL)
        permuted.deviance.reduction.list<-vector(mode="numeric",length=nsim)
        perm.JTable<-JTable
        perm.ind.variables<-manip.ind.variables
        for(k in 1:nsim)
        {
          perm.ind.variables[,(j+2)]<-sample(perm.ind.variables[,(j+2)],length(perm.ind.variables[,(j+2)]))
          pair<-1
          for(i_site in 1:(ncol(dissimilarity.data)-1))
          {
            for(j_site in (i_site+1):nrow(dissimilarity.data))
            {
              perm.JTable[pair,(6+j)]<-perm.ind.variables[i_site,(2+j)]
              perm.JTable[pair,(6+nEvar.still.in+j)]<-perm.ind.variables[j_site,(2+j)]
              pair<-pair+1 
            } ## end for j_site
          } ## end for i_site
          JTable.F<-as.data.frame(perm.JTable)
          perm.j.gdm<-gdm.fitfromtable(JTable.F,geo=geo.still.in, splines=NULL, quantiles=NULL) 
          permuted.deviance.reduction.list[k]<-gdm.without.j$gdmdeviance-perm.j.gdm$gdmdeviance          
        } ## end for k
        Deviance.reduction.variables[j,(v*2)]<-gdm.without.j$gdmdeviance-all.variables.gdm$gdmdeviance
        Pvalues.variables[j,(v*2)]<-sum(permuted.deviance.reduction.list>=(gdm.without.j$gdmdeviance-all.variables.gdm$gdmdeviance))/nsim
      } ## end for j 
    } ## end else geo.still.in==TRUE
    if(single.model)
    {
      break
    } ## end if single.model
    ## This is where we pick the variable to omit, based on P-value, then deviance reduction
    Temp.Ps<-as.numeric(Pvalues.variables[c(1:n.var.still.in),(v*2)])
    omit.var<-which.max(Temp.Ps)
    Temp.Devs<-as.numeric(Deviance.reduction.variables[c(1:n.var.still.in),(v*2)])
    for(i.check in 1:n.var.still.in)
    {
      if(Temp.Ps[i.check] == Temp.Ps[omit.var])
      {
        if(Temp.Devs[i.check] < Temp.Devs[omit.var]) omit.var<-i.check
      } ## end if Temp.Ps[i.check] == Temp.Ps[omit.var]
    } ## end for i.check
    if(geo.still.in==TRUE)
    {
      if(omit.var==1)
      {
        geo.still.in<-FALSE 
        n.var.still.in<-n.var.still.in - 1
      } ## end if omit.var==1
      else
      {
        manip.ind.variables<-manip.ind.variables[,-(omit.var+1)] 
        JTable<-JTable[,-c((omit.var+5),(omit.var+5+nEvar.still.in))]  					
        nEvar.still.in<-nEvar.still.in - 1
        n.var.still.in<-n.var.still.in - 1
      } ## end else omit.var==1
    } ## end if geo.still.in==TRUE
    else
    { 
      manip.ind.variables<-manip.ind.variables[,-(omit.var+2)]
      JTable<-JTable[,-c((omit.var+6),(omit.var+6+nEvar.still.in))]						
      nEvar.still.in<-nEvar.still.in - 1
      n.var.still.in<-n.var.still.in - 1
    } ## end else geo.still.in==TRUE
  } ## end for v (n.variables loop)
  if(single.model)
  {
    Model.test.values<-Model.test.values[,1]
    Deviance.reduction.variables<-Deviance.reduction.variables[,c(1,2)]
    Pvalues.variables<-Pvalues.variables[,c(1,2)]
  } ## end if single.model
  list(Model.test.values=Model.test.values,Deviance.reduction.variables=Deviance.reduction.variables,Pvalues.variables=Pvalues.variables)
} ## END FUNCTION


#################################################################################
## (3) APPLY THE SIGNIFICANCE TEST FUNCTION ##
## Now run the model and variable significance test function for GDM ##

## Set geo=TRUE to specify whether geographic distance is included as an independent variable
## in "ind.variables". If so, this must be as x and y coordinates in columns 1 & 2 (respectfully)
## Set single.model=TRUE to test only one model, or FALSE to run a backward elimination variable
## selection procedure. Set the number of permutations to apply in the significance testing (as a
## default this has been set to 1000).

## This function will generate three sets of data. The first are the "model.test.values", which are the 
## model deviance and model significance (as assessed by permuting sites within "ind.variables".
## The second set of data are the deviance reduction values for each variable. These values for each
## variable are the amount of model deviance reduction achieved by adding that variable to a model 
## that does not contain that variable. The third set of data are the significance values of the
## deviance reduction values reported above (i.e. the proportion of random permutations of variable i
## that generate a greater reduction in model deviance than is achieved by including variable i in the model.
## the lower the p-value, the better that variable is at reducing model deviance, in comparison to random data.

output<-gdm.model.variable.test(ind.variables, dissimilarity.data, geo=TRUE, splines=NULL, quantiles=NULL, single.model=FALSE ,permutations=5)
output

## You can also write the results out so they can be copied and pasted into Excel ##

write.csv(output$Model.test.values,"model_test_values.csv")
write.csv(output$Deviance.reduction.variables,"Deviance.reduction.variables.csv")
write.csv(output$Pvalues.variables,"Pvalues.variables.csv")

#####################################################################################
