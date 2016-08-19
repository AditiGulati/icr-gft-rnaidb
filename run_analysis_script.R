#
## call zScreen
#

#setwd( "/home/agulati/data/New_screens_from_July_2014_onwards/Breast_BT549_KS_TS_384_template_2014-07-01" )
#getwd()

source( "/Rnaidb_git/icr-gft-rnaidb/xls2platetxtfiles.R" )
#source( "/Rnaidb_git/icr-gft-rnaidb/zScreen.R" )
source( "/Rnaidb_git/icr-gft-rnaidb/zScreenSublib.R" )
source( "/Rnaidb_git/icr-gft-rnaidb/xls2txt_file.R" )
source( "/Rnaidb_git/icr-gft-rnaidb/Zscreen_qc.R")
source( "/Rnaidb_git/icr-gft-rnaidb/separateZprime.R")

#
## Get the guide file with information on new screen 
#

guideFile<-list.files(pattern="^.*guide_file.*.txt$")

## check for the correct guide file ##

#write.table(guideFile, file="/home/agulati/scripts/res.txt", sep="\t")

#guide <- read.table("/home/agulati/data/New_screens_from_July_2014_onwards/Breast_BT549_KS_TS_384_template_2014-07-01/Breast_BT549_KS_TS_384_template_2014-07-01_guide_file.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
guide = read.table( guideFile, header=TRUE, sep="\t", stringsAsFactors=FALSE )

Xls=guide$xls_file
Datapath=guide$Datapath
Cell_line=guide$Cell_line
Negcontrols = "(?i)^neg$|^siCON1$|^siCON2$|^allstar$"
Poscontrols = "(?i)^pos$|^siPLK1$|plk1"
Annotationfile=guide$Template_library_file
DescripFile=guide$Descrip_file		
Platelist=guide$Platelist_file
Plateconf=guide$Plateconf_file
ReportHTML=guide$report_html
ZscoreName=guide$zscore_file
SummaryName=guide$summary_file
ZprimeName=guide$zprime_file
ReportdirName=guide$reportdir_file
Controls_qc=guide$qc_file
Qc_plot_1=guide$plot_1_file
Qc_plot_2=guide$plot_2_file
Qc_plot_3=guide$plot_3_file
Corr_coeff=guide$corr_file
separateZprimeFile=guide$separate_zprime_file

#
## function using gdata for xls2platelistfiles conversion
#

res_1 <- xls2platelist(
  datapath=Datapath,
  xls=Xls
)

#
## function for cellHTS2 analysis
#

res_2 <- zScreenSublib(
  name=Cell_line,
  datapath=Datapath,
  poscontrols=Poscontrols,
  negcontrols=Negcontrols,
  descripFile=DescripFile,
  annotationfile=Annotationfile,
  reportHTML=ReportHTML,
  plateconf=Plateconf,
  platelist=Platelist,
  replicate_summary="median",
  zscoreName=ZscoreName, 
  summaryName=SummaryName,
  zprimeName=ZprimeName,
  reportdirName=ReportdirName
)

#
## function using gdata for converting xls to txt file for storing in the database
#

#res_3 <- xls2txt_file(
#		 datapath=Datapath,
#		 xls=Xls)

separateZprime(
  name= Cell_line,
  datapath= Datapath,
  poscontrols= Poscontrols,
  negcontrols= Negcontrols,
  descripFile= DescripFile,
  plateconf= Plateconf,
  platelist= Platelist,
  separateZprimeFile = separateZprimeFile
)

#
## function for qc - creating boxplots and calculating correlation coefficient
#
		 
res_4 <- Zscreen_qc(
  datapath=Datapath,
  summaryName=SummaryName,
  controls_qc=Controls_qc,
  qc_plot_1=Qc_plot_1,
  qc_plot_2=Qc_plot_2,
  qc_plot_3=Qc_plot_3,
  corr_coeff=Corr_coeff
)

