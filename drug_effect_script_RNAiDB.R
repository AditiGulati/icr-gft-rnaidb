##
## Drug effect analysis script
##

inputFiles = list.files()

file <- grep("DMSO", inputFiles, ignore.case=T)

for (i in 1:length(inputFiles)) {
  if (i == file) {
    dmsoFilePath = inputFiles[i]
  }else{
    drugFilePath = inputFiles[i]
  }
}

drug_effect <- NULL

# read in the dmso file
dmso <- read.table(
  dmsoFilePath,
  head=T,
  sep="\t",
  check.names=FALSE,
  stringsAsFactors=FALSE
)
dmso <- dmso[ order(dmso[,1],dmso[,4]), ]

# read in the dmso file
drug <- read.table(
  drugFilePath,
  head=T,
  sep="\t",
  check.names=FALSE,
  stringsAsFactors=FALSE
)
drug <- drug[ order(drug[,1],drug[,4]), ]

dmso$wellAnno <-  toupper(dmso$wellAnno)

# populate the output file
value <- NULL
Sample <- NULL
for (value in 1:length(dmso$wellAnno)) {
  sample <- NULL
  if (
    dmso$wellAnno[value] == "SAMPLE"
    #dmso$wellAnno[value] == "sample" | dmso$wellAnno[value] == "SAMPLE"
  ) {
    sample[value] = "SAMPLE"
  }
  if (
    dmso$wellAnno[value] == "NEG" |
    dmso$wellAnno[value] == "SICON1" |
    dmso$wellAnno[value] == "SICON2" |
    dmso$wellAnno[value] == "ALLSTAR" |
    dmso$wellAnno[value] == "SIALLSTAR" |
    dmso$wellAnno[value] == "SIPLK1" |
    dmso$wellAnno[value] == "POS" |
    dmso$wellAnno[value] == "PLK1" |
    dmso$wellAnno[value] == "SICON1 DRUG" |
    dmso$wellAnno[value] == "SICON2 DRUG" |
    dmso$wellAnno[value] == "SICON1_DRUG" |
    dmso$wellAnno[value] == "SICON2_DRUG" |
    dmso$wellAnno[value] == "ALLSTAR DRUG" |
    dmso$wellAnno[value] == "ALLSTAR_DRUG" |
    dmso$wellAnno[value] == "PLK1 DRUG" | 
    dmso$wellAnno[value] == "SIPLK1_DRUG"
  ) {    
    sample[value] = "CONTROL" 
  }
  if (
    dmso$wellAnno[value] == "EMPTY"
  ) {    
    sample[value] = "EMPTY"
  }
  Sample <- rbind(Sample, sample[value])
}

value <- NULL
for (value in 1:length(dmso$wellAnno)) {
  if (
    dmso$wellAnno[value] == "NEG" |
    dmso$wellAnno[value] == "SICON1" |
    dmso$wellAnno[value] == "SICON2" |
    dmso$wellAnno[value] == "ALLSTAR" |
    dmso$wellAnno[value] == "SIALLSTAR" |
    dmso$wellAnno[value] == "SIPLK1" |
    dmso$wellAnno[value] == "POS" |
    dmso$wellAnno[value] == "PLK1" |
    dmso$wellAnno[value] == "SICON1 DRUG" |
    dmso$wellAnno[value] == "SICON2 DRUG" |
    dmso$wellAnno[value] == "SICON1_DRUG" |
    dmso$wellAnno[value] == "SICON2_DRUG" |
    dmso$wellAnno[value] == "ALLSTAR DRUG" |
    dmso$wellAnno[value] == "ALLSTAR_DRUG" |
    dmso$wellAnno[value] == "PLK1 DRUG" | 
    dmso$wellAnno[value] == "SIPLK1_DRUG"
  ) {  
    dmso$GeneID[value] = dmso$wellAnno[value]
  }
}
value <- NULL
for (value in 1:length(dmso$wellAnno)) {  
  if (
    dmso$wellAnno[value] == "NEG" |
    dmso$wellAnno[value] == "SICON1" |
    dmso$wellAnno[value] == "SICON2" |
    dmso$wellAnno[value] == "ALLSTAR" |
    dmso$wellAnno[value] == "SIALLSTAR" |
    dmso$wellAnno[value] == "SICON1 DRUG" |
    dmso$wellAnno[value] == "SICON2 DRUG" |
    dmso$wellAnno[value] == "SICON1_DRUG" |
    dmso$wellAnno[value] == "SICON2_DRUG" |
    dmso$wellAnno[value] == "ALLSTAR DRUG" |
    dmso$wellAnno[value] == "ALLSTAR_DRUG"  
  ) {    
    dmso$Entrez_gene_ID[value] = "NEG"
  }  
  if (
    dmso$wellAnno[value] == "PLK1" |
    dmso$wellAnno[value] == "POS" |
    dmso$wellAnno[value] == "SIPLK1" |
    dmso$wellAnno[value] == "PLK1 DRUG" | 
    dmso$wellAnno[value] == "SIPLK1_DRUG"
  ) {    
    dmso$Entrez_gene_ID[value] = "POS"
  }
  if (
    dmso$wellAnno[value] == "EMPTY"
  ) {    
    dmso$Entrez_gene_ID[value] = "EMPTY"
  }
}

# add columns 1 to 12 
drug_effect <- data.frame(
  dmso[ ,"plate"], 
  dmso[ ,"well"], 
  dmso$GeneID,
  dmso$Entrez_gene_ID,
  dmso$sublib,
  Sample, 
  dmso[ ,"raw_r1_ch1"], 
  dmso[ ,"raw_r2_ch1"], 
  dmso[ ,"raw_r3_ch1"], 
  drug[ ,"raw_r1_ch1"], 
  drug[ ,"raw_r2_ch1"], 
  drug[ ,"raw_r3_ch1"],
  stringsAsFactors=FALSE
)

colnames(drug_effect) <- cbind(
  "Plate", 
  "Well",
  "Mature_Sanger_ID",
  "PreCursor",
  "sublib",
  "Sample", 
  "dmso_r1_ch1", 
  "dmso_r2_ch1", 
  "dmso_r3_ch1", 
  "drug_r1_ch1", 
  "drug_r2_ch1", 
  "drug_r3_ch1"
)

##
## log2 transform raw scores 
## for all reps in DMSO and 
## DRUG screens
##

i <- NULL 
j <- NULL
dmso_r1_log2 <- NULL
dmso_r2_log2 <- NULL
dmso_r3_log2 <- NULL
drug_r1_log2 <- NULL
drug_r2_log2 <- NULL
drug_r3_log2 <- NULL
for (j in 7:12) {
  for (i in 1:nrow(drug_effect)) {  
    dmso_r1_log2[i] <- log2(drug_effect[i,7])
    dmso_r2_log2[i] <- log2(drug_effect[i,8])
    dmso_r3_log2[i] <- log2(drug_effect[i,9])
    drug_r1_log2[i] <- log2(drug_effect[i,10])
    drug_r2_log2[i] <- log2(drug_effect[i,11])
    drug_r3_log2[i] <- log2(drug_effect[i,12])
  }
}

## add columns 12-17 to the drug_effect file ##
drug_effect <- cbind(
  drug_effect,                      
  dmso_r1_log2, 
  dmso_r2_log2, 
  dmso_r3_log2, 
  drug_r1_log2, 
  drug_r2_log2, 
  drug_r3_log2
)
drug_effect <- data.frame(drug_effect)

##
## 1. calculate per plate median for each log2 rep in
##    DMSO and DRUG screens excluding controls so that 
##    the median is not skewed by extreme effects
##
## 2. Add rows with plate centered values
##    for each of DMSO and DRUG reps
##    to the drug_effect file
##    

sublibs <- levels(as.factor(drug_effect$sublib))
sublibs <- na.omit(sublibs)

plates <- NULL
plates <- levels(as.factor(drug_effect$Plate))
plate <- NULL
plate_sample_rows <- NULL
plate_control_rows <- NULL
plate_sample_rows_for_controls <- NULL
pp_sample_median <- NULL

# calculate plate centered log2 values for samples only

plate_centered_log_intensities <- matrix(data=rep(NA, times=6*length(drug_effect$Sample)), ncol=6, nrow=length(drug_effect$Sample)) # assume we have 2*3 reps
column <- NULL
log2_columns <- NULL
log2_columns <- cbind(drug_effect[,13:18])
for(column in 1:ncol(log2_columns)) {
  for (plate in plates) {
    plate_control_rows <- which(
      drug_effect$Plate == plate & 
        drug_effect$Sample != "SAMPLE" &
        drug_effect$Sample != "EMPTY"
    )
    plate_sample_rows_for_controls <- which(
      drug_effect$Plate == plate & 
        drug_effect$Sample == "SAMPLE"
    )
    pp_control_median <- median(log2_columns[plate_sample_rows_for_controls, column],na.rm=TRUE) 
    row <- NULL
    for (row in plate_control_rows) {
      plate_centered_log_intensities[row,column] <- (log2_columns[row,column] - pp_control_median)
    }
    for (sublib in sublibs) {
      plate_sample_rows <- which(
        drug_effect$Plate == plate & 
          drug_effect$Sample == "SAMPLE" &
          drug_effect$sublib == sublib
      )
      pp_sample_median <- median(log2_columns[plate_sample_rows, column],na.rm=TRUE) 
      row <- NULL
      for (row in plate_sample_rows) {
        plate_centered_log_intensities[row,column] <- (log2_columns[row,column] - pp_sample_median)
      }
    }
  }
}
#add columns 19-24 to the drug_effect file
drug_effect <- cbind(drug_effect, plate_centered_log_intensities)

names(drug_effect)[19] <- paste("dmso_r1_log2_plate_centered") 
names(drug_effect)[20] <- paste("dmso_r2_log2_plate_centered") 
names(drug_effect)[21] <- paste("dmso_r3_log2_plate_centered") 
names(drug_effect)[22] <- paste("drug_r1_log2_plate_centered") 
names(drug_effect)[23] <- paste("drug_r2_log2_plate_centered") 
names(drug_effect)[24] <- paste("drug_r3_log2_plate_centered") 

##
## Calculate zscores  
## for each rep
##

column <- NULL
zscore_columns <- matrix(data=rep(NA, times=6*length(drug_effect$Sample)), ncol=6, nrow=length(drug_effect$Sample))
plate_centered_columns <- cbind(drug_effect[,19:24])

#control_MEDIAN <- NULL
control_MAD <- NULL

#sample_MEDIAN <- NULL
sample_MAD <- NULL

sample_rows_for_control <- NULL

for (column in 1:ncol(plate_centered_columns)) {
  control_rows <- which(
    drug_effect$Sample != "SAMPLE" & 
      drug_effect$Sample != "EMPTY"
  )
  sample_rows_for_control <- which(
    drug_effect$Sample == "SAMPLE"
  )
  #control_MEDIAN <- median(plate_centered_columns[sample_rows_for_control,column], na.rm = TRUE)
  #print(control_MEDIAN)
  control_MAD <- mad(plate_centered_columns[sample_rows_for_control,column], na.rm = TRUE)
  i <- NULL
  for (i in control_rows) {
    zscore_columns[i,column] <- (plate_centered_columns[i,column])/control_MAD
  }
  for (sublib in sublibs) {
    sample_rows <- NULL
    zscore_sublib_column <- NULL
    sample_rows <- which(
      drug_effect$Sample == "SAMPLE" & 
        drug_effect$sublib == sublib)
    #sample_MEDIAN <- median(plate_centered_columns[sample_rows,column], na.rm = TRUE)
    #print(sample_MEDIAN)
    sample_MAD <- mad(plate_centered_columns[sample_rows,column], na.rm = TRUE)
    i <- NULL
    for (i in sample_rows) {
      zscore_columns[i,column] <- (plate_centered_columns[i,column])/sample_MAD
    }
  }
}
#add columns 25_30 to the drug_effect file
drug_effect <- cbind(drug_effect, zscore_columns)

names(drug_effect)[25] <- paste("dmso_r1_zscore") 
names(drug_effect)[26] <- paste("dmso_r2_zscore")
names(drug_effect)[27] <- paste("dmso_r3_zscore")  
names(drug_effect)[28] <- paste("drug_r1_zscore")
names(drug_effect)[29] <- paste("drug_r2_zscore") 
names(drug_effect)[30] <- paste("drug_r3_zscore")  

##
## calculate median of   
## zscores in DMSO
##

dmso_z_median <- NULL
drug_z_median <- NULL
de_z_median <- NULL

columns25_27 <- cbind(drug_effect[,25:27])
columns25_27 <- as.matrix(columns25_27)
row <- NULL
for (row in 1:nrow(columns25_27)) {
  drug_z_med <- NULL
  dmso_z_med <- median(columns25_27[row,], na.rm=TRUE)
  dmso_z_median <- rbind(dmso_z_median, dmso_z_med)
  colnames(dmso_z_median) <- "dmso_median_Zscore"
}

columns28_30 <- cbind(drug_effect[,28:30])
columns28_30 <- as.matrix(columns28_30)
row <- NULL
for (row in 1:nrow(columns28_30)) {
  drug_z_med <- NULL
  drug_z_med <- median(columns28_30[row,], na.rm=TRUE)
  drug_z_median <- rbind(drug_z_median, drug_z_med)
  colnames(drug_z_median) <- "drug_median_Zscore"
}

#add column 31 & 32 to the drug_effect file
drug_effect <- cbind(drug_effect, dmso_z_median, drug_z_median)

##
## Calculate DRUG EFFECT (DE) 
## 3*(drug_rep - dmso_rep)
##

median_zscores <- NULL
median_zscores <- cbind(drug_effect[31], drug_effect[32])
row <- NULL
DE_med <- NULL
drug_value <- NULL
dmso_value <- NULL
for (row in 1:nrow(median_zscores)){
  de_med <- NULL
  drug_value <- median_zscores[row,2]
  dmso_value <- median_zscores[row,1]
  de_med <- (drug_value) - (dmso_value)
  DE_med <- rbind(DE_med, de_med)
  nrow(DE_med)
  colnames(DE_med) <- "DE_median_Zscore" 
}
## add column 25 to the drug_effect file ##
drug_effect <- cbind(drug_effect, DE_med)

##
## round off the values in drug effect file to 3 decimal places
##

drug_effect$dmso_median_zscore <- round(drug_effect$dmso_median_Zscore, digits=2)
drug_effect$drug_median_zscore <- round(drug_effect$drug_median_Zscore, digits=2)
drug_effect$DE_median_zscore <- round(drug_effect$DE_median_Zscore, digits=2)

##
## write drug_effect results to a text file
##

screen_DE_dir_name = gsub("_summary.txt","",drugFilePath)
screen_DE_dir_name = paste("drugEffect_", screen_DE_dir_name, sep="")

write.table(
  drug_effect, 
  file=paste(screen_DE_dir_name,".txt", sep=""), 
  row.names=FALSE,
  sep="\t"
)