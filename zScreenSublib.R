source("/Rnaidb_git/icr-gft-rnaidb/sublibAnalysis.R")
zScreenSublib<-function(
		  name, 
		  datapath, 
		  poscontrols="pos",
		  negcontrols="neg", 
		  annotate=TRUE, 
		  annotationfile="KS_TS_384_Template.txt", 
		  descripFile="Description.txt", 
		  plateconf="plateconf_384.txt", 
		  screenlog="Screenlog.txt", 
		  platelist="platelist_384.txt", 
		  reportHTML=TRUE, 
		  replicate_summary="median", 
		  summaryName="summary.txt", 
		  zscoreName="zscore.txt", 
		  zprimeName="zprime.txt", 
		  reportdirName="zscore"
		  ){

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
  	# annotate
  	if(annotate){
  	# annotate the genes (necessary for sublib analysis)
  	#cat("Annotating from", annotationfile, "\n")
  		x<-cellHTS2::annotate(x, geneIDFile=annotationfile, path=datapath)
		compounds<-geneAnno(x)
  	}else{
		stop ("Annotation is necessary for sublibrary analysis!")
  	}
  	annotation<-pData(featureData(x))
  	if(length(grep("sublib", colnames(annotation)))<1){stop ("Annotation should have a sublib column. Please replace sub_lib or Sub_lib with sublib!")}
  	indNonSublib <-which(is.na(annotation$sublib) | (annotation$controlStatus != "sample"))
  	
  	###########################################################	
 	# normalize
 	xn_combined<-normalizePlates(x,scale="multiplicative", log=TRUE, method="median", varianceAdjust = "none", negControls=negcontrols, posControls=poscontrols);
  	xn<-normalizePlatesSublib(x,scale="multiplicative", log=TRUE, method="median", varianceAdjust = "none", negControls=negcontrols, posControls=poscontrols);
	Data(xn)[indNonSublib,,]<-Data(xn_combined)[indNonSublib,,]
	
  	###########################################################	
  	# calculate zscores
  	xsc_combined<-scoreReplicates(xn_combined, method="zscore", sign="+")
  	xsc_combined<-summarizeReplicates(xsc_combined, summary=replicate_summary)
  	
  	xsc <- scoreReplicatesSublib(xn, method="zscore", sign="+")
  	xsc<-summarizeReplicates(xsc, summary=replicate_summary)
	Data(xsc)[indNonSublib,,]<-Data(xsc_combined)[indNonSublib,,]
	
  	# top table
  	summaryFile <- paste(datapath,summaryName, sep="")
  		
  	summary_info<-getTopTable(list("raw"=x, "normalized"=xn, "scored"=xsc), file=summaryFile)
  	setSettings(list(
        plateList=list(reproducibility=list(include=TRUE, map=TRUE),
      	intensities=list(include=TRUE, map=TRUE)),
      	screenSummary=list(scores=list(range=c(-20, 10), map=TRUE)))) 

  	###########################################################	
  	# write a QC report for zscores
  	if(reportHTML){
  		reportdir<-paste(datapath, reportdirName, sep="")
  		writeReport(raw=x, normalized=xn, scored=xsc, outdir=reportdir, force=TRUE, posControls=poscontrols, negControls=negcontrols, mainScriptFile="/Rnaidb_git/icr-gft-rnaidb/zScreenSublib.R")
  	}

  	###########################################################
  	# write z scores for this screen in the same folder
  	scorefile<-paste(datapath, zscoreName, sep="")

  	plates<-plate(xsc)
  	wells<-well(xsc)
  	scores<-Data(xsc)
  	# prepare a simple text report
  	combinedz<-data.frame(compound=compounds, plate=plates, well=wells, zscore=scores)
  	names(combinedz)<-c("Compound", "Plate", "Well", "Zscore")
  	write.table(combinedz, scorefile, sep="\t", quote=FALSE, row.names=FALSE)

  	###########################################################
 	z_summary <- read.table(
    summaryFile,
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE
  )
  
  sublibs <- na.omit(unique(z_summary$sublib))

  summary_with_z <- NULL
  summary_with_z <- cbind(
    z_summary, 
    zscore_rep1=rep(NA, times=nrow(z_summary)),
    zscore_rep2=rep(NA, times=nrow(z_summary)),
    zscore_rep3=rep(NA, times=nrow(z_summary))
  )

  for (sublib in sublibs) {
    lib_rows <- which(z_summary$sublib == sublib)
    sublib_mad_rep1 <- mad(z_summary$normalized_r1_ch1[lib_rows], na.rm=TRUE)
    sublib_mad_rep2 <- mad(z_summary$normalized_r2_ch1[lib_rows], na.rm=TRUE)
    sublib_mad_rep3 <- mad(z_summary$normalized_r3_ch1[lib_rows], na.rm=TRUE)
    summary_with_z$zscore_rep1[lib_rows] <- z_summary$normalized_r1_ch1[lib_rows]/sublib_mad_rep1
    summary_with_z$zscore_rep1[lib_rows] <- round(summary_with_z$zscore_rep1[lib_rows], digits=3)
    summary_with_z$zscore_rep2[lib_rows] <- z_summary$normalized_r2_ch1[lib_rows]/sublib_mad_rep2
    summary_with_z$zscore_rep2[lib_rows] <- round(summary_with_z$zscore_rep2[lib_rows], digits=3)
    summary_with_z$zscore_rep3[lib_rows] <- z_summary$normalized_r3_ch1[lib_rows]/sublib_mad_rep3
    summary_with_z$zscore_rep3[lib_rows] <- round(summary_with_z$zscore_rep3[lib_rows], digits=3)
  }
  z_summaryName <- gsub(".txt","_with_rep_zscores.txt", summaryName, fixed=TRUE)
  write.table(
    summary_with_z,
    file=z_summaryName,
    col.names=TRUE,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )
  	###########################################################
  	# write zprime for this screen in the same folder	
  	zp <- getZfactor(xn,
           posControls=poscontrols,
           negControls=negcontrols)     

  	zprimefile<-paste(datapath, zprimeName, sep="")
  	write.table(zp, zprimefile, sep = "\t", row.names=FALSE)
}

