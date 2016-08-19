separateZprime<-function(name, 
		  datapath, 
		  poscontrols="pos",
		  negcontrols="neg", 
		  descripFile="Description.txt", 
		  plateconf="plateconf_384.txt", 
		  screenlog="Screenlog.txt", 
		  platelist="platelist_384.txt", 
		  separateZprimeFile="_separate_zprime.txt"
		  )
{
  	require(cellHTS2)

  	# add a trailing slash if necessary
  	datapath=paste(datapath,"/", sep="")

  	# keep this very compact
  	# raw data
  	x<-readPlateList(platelist, name=name, path=datapath);
  	#	configure
  	cat("-----",name,"------ Z score\n")
  	# check if we have a screenlog 
  	if (file.exists( paste(datapath,screenlog,sep="") )){
		cat("Screenlog found.\n")
		x<-configure(x, descripFile=descripFile, confFile=plateconf, logFile=screenlog, path=datapath);
  	} else{
		cat("No Screenlog found. Proceeding without.\n")
		x<-configure(x, descripFile=descripFile, confFile=plateconf, path=datapath);
  	}


  	###########################################################	
  	# normalize	
  	xn<-normalizePlates(x,scale="multiplicative", log=TRUE, method="median", varianceAdjust = "none", negControls=negcontrols, posControls=poscontrols);

  	###########################################################
  	# write zprime for this screen in the same folder
  	nrSamples <- dim(Data(xn))[2]	
  	annotation<-pData(featureData(x))
	controlStatus<-as.character(unique(annotation$controlStatus))
	poscontrols <-"(?i)^siPLK1$|plk1"
	idxPos<-grep(poscontrols, controlStatus)
	if(length(idxPos) < 1) stop("There is no PLK1 positive control!")
	
	negcontrols<-"(?i)^siCON1$"
	idxNeg<-grep(negcontrols, controlStatus)
	result<-data.frame(replicates=paste("Replicate", 1: nrSamples, sep=""), zp_sicon1=rep(NA, nrSamples), zp_sicon2=rep(NA, nrSamples), zp_allstar=rep(NA, nrSamples))
	if(length(idxNeg) > 0)
	{
		zp_sicon1<-	getZfactor(xn,
           posControls=poscontrols,
           negControls=negcontrols) 
        zp_sicon1<-as.data.frame(zp_sicon1)
        result$zp_sicon1[1:nrow(zp_sicon1)]<-zp_sicon1$Channel1[1:nrow(zp_sicon1)]
	}
	
	negcontrols<-"(?i)^siCON2$"
	idxNeg<-grep(negcontrols, controlStatus)
	if(length(idxNeg) > 0)
	{		
		zp_sicon2<-	getZfactor(xn,
           posControls=poscontrols,
           negControls=negcontrols) 
       	zp_sicon2<-as.data.frame(zp_sicon2)
        result$zp_sicon2[1:nrow(zp_sicon2)]<-zp_sicon2$Channel1[1:nrow(zp_sicon2)]
	}
		
	negcontrols<-"(?i)^allstar$"
	idxNeg<-grep(negcontrols, controlStatus)
	if(length(idxNeg) > 0)
	{		
		zp_allstar<- getZfactor(xn,
           posControls=poscontrols,
           negControls=negcontrols) 
        zp_allstar<-as.data.frame(zp_allstar)
        result$zp_allstar[1:nrow(zp_allstar)]<-zp_allstar$Channel1[1:nrow(zp_allstar)]
	}
	zprimefile<-paste(datapath, separateZprimeFile, sep="")
	write.table(result, zprimefile, sep = "\t", row.names=FALSE)	
}





