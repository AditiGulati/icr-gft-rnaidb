xls2platelist<-function(xls, datapath){
	require(gdata)
	setwd(datapath)
	xls<-unlist(strsplit(xls, split=" "))
	screen_plate_number <- 0
	screen_name<-basename(datapath)
	for (iFile in 1:length(xls))
	{
		mydata<-read.xls(xls[iFile],1)
		plates<-unique(mydata$Plate)	
		for (plate in plates)
		{
			temp<-mydata[mydata$Plate==plate,c(1,3,6)]
			screen_plate_number <- screen_plate_number + 1
    		filename<-paste(screen_name,"_P",screen_plate_number,".txt", sep="")
			write.table(temp, file=filename, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
		}
	}
}
