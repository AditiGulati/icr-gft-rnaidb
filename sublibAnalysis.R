normalizePlatesSublib <- function(object, scale="additive", log = FALSE, method="median", varianceAdjust="none",
                            posControls, negControls,...)
{

    if(!is(object, "cellHTS"))
        stop("'object' should be of class 'cellHTS'.")
    ## Check the status of the 'cellHTS' object
    if(!state(object)[["configured"]])
        stop("Please configure 'object' (using the function 'configure') before normalization.")

    ## Check the conformity between the scale of the data and the chosen preprocessing
    if(scale=="additive" & log)
        stop("For data on the 'additive' scale, please do not set 'log=TRUE'. ",
             "Please have a look at the documentation of the 'scale' and 'log' options ",
             "of the 'normalizePlates' function.") 

    if(!(varianceAdjust %in% c("none", "byPlate", "byBatch", "byExperiment"))) 
        stop(sprintf("Undefined value %s for 'varianceAdjust'.", varianceAdjust))

    ## 1. Log transformation: 
    oldRawData <- Data(object)
    if(log)
    {
        Data(object) <- suppressWarnings(log2(oldRawData))
        if(any(oldRawData[!is.na(oldRawData)]==0))
            warning("Data contains 0 values.\n",
                    "Log transformation for those values resulted in -Inf",
                    call.=FALSE)
        if(min(oldRawData, na.rm=TRUE)<0)
            warning("Data contains negative values.\n",
                    "Log transformation for those values resulted in NA",
                    call.=FALSE) 
        scale <- "additive"
    }
    
    ## 2. Plate-by-plate adjustment:
    allowedFunctions <- c("mean", "median", "shorth", "negatives", "POC", "NPI", "Bscore", "locfit")
    ## overwrite assayData with the new data 
    object <- switch(method,
                     "mean" = perPlateScaling(object, scale, method),
                     "median" = perPlateSublibScaling(object, scale, method),
                     "shorth" = perPlateScaling(object, scale, method),
                     "negatives" = perPlateScaling(object, scale, method, negControls),
                     "POC" = controlsBasedNormalization(object, method, posControls, negControls),
                     "NPI" = controlsBasedNormalization(object, method, posControls, negControls),
                     "Bscore" = Bscore(object, ...),
                     "locfit" = spatialNormalization(object, ...), 
                     stop(sprintf("Invalid value '%s' for argument 'method'.\n Allowed values are: %s.", 
                                  method, paste(allowedFunctions, collapse=", ")))
                     )

    ## 3. Variance adjustment (optional):
    if(varianceAdjust!="none")
        object <- adjustVariance(object, method=varianceAdjust)

    object@state[["normalized"]] <- TRUE
    object@processingInfo[["normalized"]] <- method
    validObject(object)
    return(object)
}

## ===========================================================
## 		----	perPlateSublibScaling ------
##
## perPlateSublibScaling centres each sublibrary on its plate median!

perPlateSublibScaling <- function(object, scale, stats="median", negControls){
	if(stats!="median") stop ("perPlateSublibScaling centres each sublibrary on its plate median!");
  	xnorm <- Data(object)
  	d <- dim(xnorm)
  	nrPlates <- max(plate(object))
  	nrSamples <- d[2]
  	nrChannels <- d[3]
  	annotation<-pData(featureData(object))
  	for(p in 1:nrPlates) 
  	{
    	indSample <- which((annotation$controlStatus == "sample") & (annotation$plate==p))
    	indControl <- which((annotation$controlStatus != "sample") & (annotation$plate==p))
    	for(ch in 1:nrChannels)
    	{
        	for(r in 1:nrSamples) 
        	{
          		if(!all(is.na(xnorm[indSample, r, ch])))
          		{
          			xnorm[indControl, r, ch] <- xnorm[indControl, r, ch]- median(xnorm[indSample, r, ch], na.rm=TRUE)
          		}
  			}# r
		}# ch	
  	}# p
  	
  	for(p in 1:nrPlates) 
  	{
    	plateInds <- which(annotation$plate==p)
    	plate_sublib<-annotation$sublib[plateInds]
    	sublib_name<-setdiff(plate_sublib, NA)
    	#sublib_name<-unique(plate_sublib)
    	for(isublib in 1:length(sublib_name))
    	{
    		inds <- which((annotation$sublib == sublib_name[isublib]) & (annotation$plate==p) & (annotation$controlStatus == "sample"))
    		for(ch in 1:nrChannels)
    		{
        		for(r in 1:nrSamples) 
        		{
          			if(!all(is.na(xnorm[inds, r, ch])))
          			{
          				xnorm[inds,r,ch] <- xnorm[inds, r, ch]- median(xnorm[inds, r, ch], na.rm=TRUE)
          			}
  				}# r
			}# ch
		}#isublib		
  	}# p
  
  	Data(object) <- xnorm
  	return(object)
}

scoreReplicatesSublib <- function(object, sign="+", method="zscore", ...)
{
    methodArgs <- list(...)
    ## 1) Score each replicate using the selected method:
    xnorm <- if(method=="none") Data(object) else do.call(paste("scoreReplicatesSublib", method, sep="By"),
                args=c(list(object), methodArgs))
    ## Store the scores in 'assayData' slot.

    ## 2) Use "sign" to make the meaning of the replicates summarization
    ## independent of the type of the assay
    sg <- switch(sign,
                 "+" = 1,
                 "-" = -1,
                 stop(sprintf("Invalid value '%s' for argument 'sign'", sign)))

    Data(object) <- sg*xnorm
    validObject(object)
    object@processingInfo[["scored"]] <- method
    return(object)
}

scoreReplicatesSublibByzscore <- function(object)
{
    xnorm <- Data(object)
    d <- dim(xnorm)
  	nrSamples <- d[2]
  	nrChannels <- d[3]
    annotation<-pData(featureData(object))
    
    indSample <- which((annotation$controlStatus == "sample"))
    indControl <- which((annotation$controlStatus != "sample"))
    for(ch in 1:nrChannels)
    {
        for(r in 1:nrSamples) 
        {
          	xnorm[indControl, r, ch] <- (xnorm[indControl, r, ch]- median(xnorm[indSample, r, ch], na.rm=TRUE))/mad(xnorm[indSample, r, ch], na.rm=TRUE)
  		}# r
	}# ch

    sublib<-annotation$sublib
    sublib_name<-setdiff(sublib, NA)
    #sublib_name<-unique(sublib)
    for(isublib in 1:length(sublib_name))
    {
    	inds <- which((annotation$sublib == sublib_name[isublib]) & (annotation$controlStatus == "sample"))
    	for(ch in 1:nrChannels)
    	{
        	for(r in 1:nrSamples) 
        	{
          		xnorm[inds,r,ch] <- (xnorm[inds, r, ch]- median(xnorm[inds, r, ch], na.rm=TRUE))/mad(xnorm[inds, r, ch], na.rm=TRUE)
  			}# r
		}# ch
	}#isublib	
        	
    return(xnorm)
}