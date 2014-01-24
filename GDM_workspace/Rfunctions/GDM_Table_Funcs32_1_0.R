##
## edit this path to reference the GDM table dll if it is in a different location on your machine
##
"mydllpath" <- "D:\\GDM_2012\\GDM_Table_version_For_R\\GDMBinLib32.dll"


"gdm.fitfromtable" <- 
function (data, geo=FALSE, splines=NULL, quantiles=NULL) 
{
	options(warn.FPU = FALSE)
##
##	sanity check on the data table
##
	if ( ncol(data) < 6 )
		stop("Not enough columns in data")
	if ( nrow(data) < 1 )
		stop("Not enough rows in data")


##
##	check that the response data is [0..1]
##
	rtmp <- data[,1]
	if (length(rtmp[rtmp<0]) > 0)
		stop("Response data has negative values")
	if (length(rtmp[rtmp>1]) > 0)
		stop("Response data has values greater than 1")


##
##	current data format is response,weights,X0,Y0,X1,Y1 before any predictor data (thus 6 leading columns)
##
	LEADING_COLUMNS <- 6
	if ( geo ) {
	      nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2 + 1		
	}
	else {
	     	nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2
	}

	if ( nPreds < 1 )
		stop("Data has NO predictors")


##
## 	setup the predictor name list
##
	if ( geo ) {
		if ( nPreds > 1 )
			predlist <- c("Geographic", names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds-1)])
		else
			predlist <- c("Geographic")
	}
	else {
		predlist <- names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds)]
	}


##
##	deal with the splines and quantiles
##
	if (is.null(quantiles))
	{
		##
		## generate quantiles internally from the data
		##
		if ( is.null(splines) ) {
			nSplines <- 3
			quantvec <- rep(0, nPreds * nSplines)
			splinvec <- rep(nSplines, nPreds)

			if ( geo ) {
      		      ## get quantiles for the geographic distance
            		v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
		            quantvec[1] <- min(v)
      		      quantvec[2] <- median(v)
            		quantvec[3] <- max(v)

				if ( nPreds > 1 ) {
		            	## get quantiles for the environmental predictors
				      for (i in seq(from = 1, to = nPreds-1, by = 1)) {
            				v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])                 
      	            		index = i * nSplines
		                  	quantvec[index+1] <- min(v)
		      	            quantvec[index+2] <- median(v)
      	      	      	quantvec[index+3] <- max(v)
	      	      	}
				}
    			}
		
			else {
				## get quantiles for the environmental predictors after skipping geographic preds
	      	      for (i in seq(from = 1, to = nPreds, by = 1)) {
					v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])                   
            	      	index = (i-1) * nSplines
	            	      quantvec[index+1] <- min(v)
	      	            quantvec[index+2] <- median(v)
            	      	quantvec[index+3] <- max(v)
		            }
			}
		}
		
		else {
			##
			##	otherwise check that the supplied splines vector has enough data and minumum spline values of 3
			##
			if ( length(splines) != nPreds ) {
				stop(paste("splines argument has", length(splines), "items but needs", nPreds, "items."))
			}

			## count the total number of user defined splines to dimension the quantiles vector
			quantvec <- rep(0, sum(splines))
			splinvec <- splines

			if ( geo ) {
				if ( splines[1] < 3 )
					stop("Must greater than or equal to 3 splines per predictor")

				## get quantiles for the geographic distance
      	      			v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
		      	      	quantvec[1] <- min(v)			## 0% quantile
      	      			quantvec[splines[1]] <- max(v)	## 100% quantile

				quant_increment <- 1.0 / (splines[1]-1)
				this_increment <- 1
				for (i in seq(from = 2, to = (splines[1]-1), by = 1)) {  
					## mid % quantiles
					quantvec[i] <- quantile(v,quant_increment*this_increment)
					this_increment <- this_increment + 1
				}

				if ( nPreds > 1 ) {
					## get quantiles for the environmental predictors
					current_quant_index <- splines[1]
					for (i in seq(from = 1, to = nPreds-1, by = 1)) {
      	      		
						num_splines <- splines[i+1]
						if ( num_splines < 3 )
							stop("Must greater than or equal to 3 splines per predictor")

						v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])                 
						quantvec[current_quant_index+1] <- min(v)	            ## 0% quantile
	      			      		quantvec[current_quant_index+num_splines] <- max(v)	## 100% quantile

						quant_increment <- 1.0 / (num_splines-1)
						this_increment <- 1
						for (i in seq(from = 2, to = (num_splines-1), by = 1)) {  
							## mid % quantiles
							quantvec[current_quant_index+i] <- quantile(v,quant_increment*this_increment)
							this_increment <- this_increment + 1
						}
						current_quant_index <- current_quant_index + num_splines
					}
				}
			}

			else {
				## get quantiles for the environmental predictors
				current_quant_index <- 0
				for (i in seq(from = 1, to = nPreds, by = 1)) {
      	      	
					num_splines <- splines[i]
					if ( num_splines < 3 )
						stop("Must greater than or equal to 3 splines per predictor")

					v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])          
					quantvec[current_quant_index+1] <- min(v)	            ## 0% quantile
	      		      		quantvec[current_quant_index+num_splines] <- max(v)	## 100% quantile

					quant_increment <- 1.0 / (num_splines-1)
					this_increment <- 1
					for (i in seq(from = 2, to = (num_splines-1), by = 1)) {  
						## mid % quantiles
						quantvec[current_quant_index+i] <- quantile(v,quant_increment*this_increment)
						this_increment <- this_increment + 1
					}
					current_quant_index <- current_quant_index + num_splines
				}
			}
		}
	}

	else
	{
		##
		## user defined quantiles supplied as an argument
		##
		if ( is.null(splines) ) {
			## check that there are nPreds * 3 quantiles in the user defined vector
			if ( length(quantiles) != (nPreds * 3) ) {
				stop(paste("There should be", (nPreds * 3), "items in the quantiles argument, not", length(quantiles), "items."))
			}

			## now check that each of the three quantiles for each predictor are in ascending order
			for (i in seq(from = 1, to = nPreds, by = 1)) {
                        index = i * 3
				if ((quantiles[index-1] < quantiles[index-2] ) || 
                            (quantiles[index] < quantiles[index-2]) || 
                            (quantiles[index] < quantiles[index-1])) {
						stop(paste("Wrong quantiles for predictor:", predlist[i]))
				}
			}

			nSplines <- 3
			quantvec <- quantiles
			splinvec <- rep(nSplines, nPreds)

		}
		
		else {
			## check that there are sum(splines) quantiles in the user defined vector
			if ( length(quantiles) != sum(splines) ) {
				stop(paste("There should be", sum(splines), "items in the quantiles argument, not", length(quantiles), "items."))
			}
				
			## now check that each of the quantiles for each predictor are in ascending order
			index = 0
			for (i in seq(from = 1, to = nPreds, by = 1)) {
				for (j in seq(from = 2, to = splines[i], by = 1)) {
                              if (quantiles[index+j] < quantiles[index+j-1]) {
						stop(paste("Wrong quantiles for predictor:", predlist[i]))
					}
				}
				index <- index + splines[i]
			}

			quantvec <- quantiles
			splinvec <- splines
		}
	}


      	p1 <- 0
	p2 <- 0
    	p3 <- 0
    	p4 <- 0
    	p5 <- rep(0,times=length(quantvec))
      	p6 <- rep(0,times=nrow(data))
      	p7 <- rep(0,times=nrow(data))
      	p8 <- rep(0,times=nrow(data))
	

##
##  Call the dll
##
	dyn.load(mydllpath)

	

	z <- .C( "GDM_FitFromTable",
               paste(getwd()),
               as.matrix(data),
               as.integer(geo),
               as.integer(nPreds), 
               as.integer(nrow(data)), 
               as.integer(ncol(data)),
               as.integer(splinvec),
               as.double(quantvec),                 
               gdmdev = as.double(p1),
               nulldev = as.double(p2),
               expdev = as.double(p3),
               intercept = as.double(p4),         
               coeffs = as.double(p5),
               response = as.double(p6),
               preddata = as.double(p7),
               ecodist = as.double(p8)
               )

        dyn.unload(mydllpath)

	call <- match.call()
        m <- match.call(expand = F)
        
    	return(structure(list(dataname = m[[2]],
                       geo = geo,
                       sample = nrow(data),
                       gdmdeviance = z$gdmdev,
                       nulldeviance = z$nulldev,
                       explained = z$expdev,
                       intercept = z$intercept,
                       predictors = predlist,
                       coefficients = z$coeffs,
                       quantiles = quantvec,
                       splines = splinvec,
                       creationdate = date(),
                       observed = z$response,
                       predicted = z$preddata,
                       ecological = z$ecodist)))         
}





"gdm.plotfromtable" <- 
function (model, plot.layout = c(2,2), plot.color = rgb(0.0,0.0,1.0), plot.linewidth=2.0) 
{
      	options(warn.FPU = FALSE)
      	PSAMPLE <- 200
      	preddata <- rep(0,times=PSAMPLE)
			
	##
	## establish what plot layout to use
	##  
	thisplot <- 0
	one_page_per_plot <- FALSE
	if ((plot.layout[1]==1) && (plot.layout[2]==1)) 
		one_page_per_plot <- TRUE
	else
	  par(mfrow=plot.layout)	


      	##
      	## apply the link function and plot.....
      	##
      	plot( model$ecological, 
              model$observed, 
              xlab="Predicted Ecological Distance", 
              ylab="Observed Compositional Dissimilarity", type="n" )

        points( model$ecological, model$observed, pch=20, cex=0.25, col=plot.color )

        overlayX <- seq( from=min(model$ecological), to=max(model$ecological), length=PSAMPLE )
        overlayY <- 1 - exp( - overlayX )
        lines( overlayX, overlayY, lwd=plot.linewidth ) 
	thisplot <- thisplot + 1


        ##
        ## use the raw data and plot.....
        ##
	if (one_page_per_plot){
      		x11()
        	dev.next()
	}
        plot( model$predicted, 
              model$observed, 
              xlab="Predicted Compositional Dissimilarity", 
              ylab="Observed Compositional Dissimilarity", type="n" )

        points( model$predicted, model$observed, pch=20, cex=0.25, col=plot.color )

        overlayX <- overlayY <- seq( from=min(model$predicted), to=max(model$predicted), length=PSAMPLE )
        lines( overlayX, overlayY, lwd=plot.linewidth ) 
	thisplot <- thisplot + 1


        ##
        ## determine the max of all the predictor data
        ##
        dyn.load(mydllpath)
        preds <- length(model$predictors)
        predmax <- 0
        splineindex <- 1
        for ( i in 1:preds ) {  

       		## only if the sum of the coefficients associated with this predictor is > 0.....
                numsplines <- model$splines[i]
                if ( sum(model$coefficients[splineindex:(splineindex+numsplines-1)]) > 0 ) {
			
			## get predictor plot Y-data                            
                        z <- .C( "GetPredictorPlotData", 
                                 pdata = as.double(preddata),
                                 as.integer(PSAMPLE),
                                 as.double(model$coefficients[splineindex:(splineindex+numsplines-1)]),
                                 as.double(model$quantiles[splineindex:(splineindex+numsplines-1)]),
                                 as.integer(numsplines)
                                 )
                        

                        v <- max(z$pdata);                      
                        if (v > predmax ) predmax <- v
        	}
                splineindex <- splineindex + numsplines
	}


	
        ##
        ## plot the predictors with non-zero sum of coefficients
        ##      
        splineindex <- 1
	for ( i in 1:preds ) {  

      		## only if the sum of the coefficients associated with this predictor is > 0.....
                numsplines <- model$splines[i]
                if ( sum(model$coefficients[splineindex:(splineindex+numsplines-1)]) > 0 ) {
						
			if (one_page_per_plot){
      				x11()
		        	dev.next()
			}
			else {
				thisplot <- thisplot + 1
				if (thisplot > (plot.layout[1] * plot.layout[2])) {	
					x11()                              				
                        		dev.next()					
					thisplot <- 1
                        		par(mfrow=plot.layout)	
				}
			}
                  
			## get predictor plot Y-data    
                  z <- .C( "GetPredictorPlotData", 
                           pdata = as.double(preddata),
                           as.integer(PSAMPLE),
                           as.double(model$coefficients[splineindex:(splineindex+numsplines-1)]),
                           as.double(model$quantiles[splineindex:(splineindex+numsplines-1)]),
                           as.integer(numsplines)
                           )

                        
                  plot( seq(from=model$quantiles[[(i*3)-2]],to=model$quantiles[[(i*3)]], length=PSAMPLE),
                  	z$pdata, 
                        xlab=model$predictors[i], 
                        ylab=paste("f(", model$predictors[i], ")", sep="" ), 
            	  	ylim=c(0,predmax), type="l" )
      		}
	      	splineindex <- splineindex + numsplines
      	}       
	dyn.unload(mydllpath)
}






"gdm.predictfromtable" <- 
function (model, data) 
{
    options(warn.FPU = FALSE)

#
#  Call the dll
#
        dyn.load(mydllpath)

        predicted <- rep(0,times=nrow(data))
        z <- .C( "GDM_PredictFromTable",
                        as.matrix(data),
                        as.integer(model$geo),
                        as.integer(length(model$predictors)), 
                        as.integer(nrow(data)), 
                        as.double(model$quantiles),
                        as.integer(model$splines),
                        as.double(c(model$intercept,model$coefficients)),
                        preddata = as.double(predicted)
               )

        dyn.unload(mydllpath)

    return(z$preddata)
}




"gdm.summaryfromtable" <- 
function (model) 
{
        print( "", quote=F )    
        print( "", quote=F )    
        print( "GDM Modelling Summary", quote=F );
        print( paste( "Creation Date: ", model$creationdate ), quote=F );
        print( "", quote=F )    
        call <- match.call()
        m <- match.call(expand = F)
        print( paste( "Name: ", m[[2]] ), quote=F )
        print( "", quote=F )    
        print( paste( "Data: ", model$dataname ), quote=F )
        print( "", quote=F )    
        print( paste( "Samples: ", model$sample ), quote=F )
        print( "", quote=F )    
        print( paste( "Use Geographical Distance: ", model$geo ), quote=F )
        print( "", quote=F )    
        print( paste( "NULL Deviance: ", model$nulldeviance ), quote=F )
        print( paste( "GDM Deviance: ", model$gdmdeviance ), quote=F )  
        print( paste( "Deviance Explained: ", model$explained ), quote=F )
        print( "", quote=F )    
        print( paste( "Intercept: ", model$intercept ), quote=F )
        print( "", quote=F )    
        thiscoeff <- 1
        thisquant <- 1
        for ( i in 1:length(model$predictors) ) {
        print( paste( "Predictor ",i,": ",model$predictors[[i]], sep="" ), quote=F )            
        print( paste( "Splines: ",model$splines[[i]], sep="" ), quote=F )
        numsplines <- model$splines[[i]]
        for ( j in 1:numsplines ) {
                if ( j == 1 ) print( paste( "Min Quantile: ",model$quantiles[[thisquant]], sep="" ), quote=F )          
                else if ( j == numsplines ) print( paste( "Max Quantile: ",model$quantiles[[thisquant]], sep="" ), quote=F )
                else print( paste( round(100/(numsplines-1),digits=2),"% Quantile: ",model$quantiles[[thisquant]], sep="" ), quote=F )
                thisquant <- thisquant + 1
        }
        for ( j in 1:numsplines ) {
            print( paste( "Coefficient[",j,"]: ",model$coefficients[[thiscoeff]], sep="" ), quote=F )
            thiscoeff <- thiscoeff + 1
        }
        print( "", quote=F )                
    }   
}
