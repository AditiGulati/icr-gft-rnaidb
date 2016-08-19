xls2txt_file<-function(xls, datapath){
	require(gdata)
	mydata<-read.xls(xls,1)
    		temp<-mydata[,c(1,3,6)]
		filename<-paste(my.prefix,".txt",sep="")
		write.table(temp, file=filename, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}
my.file<-list.files(pattern="xls")
my.prefix<-substr((my.file),1,nchar(my.file) - 4)
