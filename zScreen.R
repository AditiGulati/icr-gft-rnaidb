# genericScreen2
# updated and improved version of genericScreen. Keep each screen separate and contained within the screen folder

#setwd( "/home/agulati/data/New_screens_from_July_2014_onwards/Breast_BT549_KS_TS_384_template_2014-07-01" )

#setwd(datapath)

#zScreen<-function(name, datapath, posControls, 
#negControls, annotate, annotationfile, descripFile, plateconf, screenlog, platelist, reportHTML, replicate_summary, summaryName, zscoreName, zprimeName, reportdirName){

zScreen<-function(name, 
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
  # normalize	
  xn<-normalizePlates(x,scale="multiplicative", log=TRUE, method="median", varianceAdjust = "none", negControls=negcontrols, posControls=poscontrols);

  ###########################################################	
  # calculate zscores
  xsc <- scoreReplicates(xn, method="zscore", sign="+")
  xsc<-summarizeReplicates(xsc, summary=replicate_summary)

  ###########################################################	
  # annotate
  if(annotate){
  # annotate the genes (not entirely necessary)
  #cat("Annotating from", annotationfile, "\n")
	xsc<-cellHTS2::annotate(xsc, geneIDFile=annotationfile, path=datapath)
	compounds<-geneAnno(xsc)
  }else{
	compounds<-rep(NA,length(wells))
  }

	# top table
	
  getTopTable(list("raw"=x, "normalized"=xn, "scored"=xsc), file=paste(datapath,summaryName, sep=""))

  ###########################################################	
  # write a QC report for zscores
  if(reportHTML){
  reportdir<-paste(datapath, reportdirName, sep="")
  writeReport(raw=x, normalized=xn, scored=xsc, outdir=reportdir, force=TRUE, posControls=poscontrols, negControls=negcontrols, mainScriptFile="/home/agulati/scripts/zScreen.R")
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
  # write zprime for this screen in the same folder	
  # calculates zprime for each replicate separately 
  zp <- getZfactor(xn,
           posControls=poscontrols,
           negControls=negcontrols)     

  zprimefile<-paste(datapath, zprimeName, sep="")
  write.table(zp, zprimefile, sep = "\t", row.names=FALSE)
}





