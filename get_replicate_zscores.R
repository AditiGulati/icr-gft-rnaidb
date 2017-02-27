##
##  add replicate Z-scores to the summary file 
##

summaryFileList <- "summary_file_list.txt"

temp_info <- NULL

summaryFiles <- read.table(
  summaryFileList,
  header=FALSE,
  sep="\t",
  stringsAsFactors=FALSE
)

temp_info <- c(
	temp_info,
	paste("nrow in file:", nrow(summaryFiles))
	)
	
summaryFile <- NULL

for (summaryFile in 1:nrow(summaryFiles)) {
  summary <- read.table(
    summaryFiles[summaryFile,],
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE
  )

  sublibs <- na.omit(unique(summary$sublib))

  temp_info <- c(
  	temp_info,
  	sublibs
  	)

  summary_with_z <- NULL
  summary_with_z <- cbind(
    summary, 
    zscore_rep1=rep(NA, times=nrow(summary)),
    zscore_rep2=rep(NA, times=nrow(summary)),
    zscore_rep3=rep(NA, times=nrow(summary))
  )
  for (sublib in sublibs) {
    lib_rows <- which(summary$sublib == sublib)
    
    sublib_mad_rep1 <- mad(summary$normalized_r1_ch1[lib_rows], na.rm=TRUE)
    sublib_mad_rep2 <- mad(summary$normalized_r2_ch1[lib_rows], na.rm=TRUE)
    sublib_mad_rep3 <- mad(summary$normalized_r3_ch1[lib_rows], na.rm=TRUE)
    
    summary_with_z$zscore_rep1[lib_rows] <- summary$normalized_r1_ch1[lib_rows]/sublib_mad_rep1
    summary_with_z$zscore_rep1[lib_rows] <- round(summary_with_z$zscore_rep1[lib_rows], digits=3)
    
    summary_with_z$zscore_rep2[lib_rows] <- summary$normalized_r2_ch1[lib_rows]/sublib_mad_rep2
    summary_with_z$zscore_rep2[lib_rows] <- round(summary_with_z$zscore_rep2[lib_rows], digits=3)
    
    summary_with_z$zscore_rep3[lib_rows] <- summary$normalized_r3_ch1[lib_rows]/sublib_mad_rep3
    summary_with_z$zscore_rep3[lib_rows] <- round(summary_with_z$zscore_rep3[lib_rows], digits=3)
  }
   
  summaryName <- gsub("_summary.txt","_summary_with_rep_zscores.txt", summaryFiles[summaryFile,], fixed=TRUE)
  
  write.table(
    summary_with_z,
    file=summaryName,
    col.names=TRUE,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )
}

write.table(
	temp_info,
	"temp_info.txt"
	)