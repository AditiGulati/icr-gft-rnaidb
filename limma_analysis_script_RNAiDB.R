require(limma)

inputFiles = list.files()

file <- grep("summaryFile1_", inputFiles, ignore.case=T)

for (i in 1:length(inputFiles)) {
  if (i == file) {
    summaryFilePath_1 = inputFiles[i]
  }else{
    summaryFilePath_2 = inputFiles[i]
  }
}

#
# limma analysis
#

isogenics<-factor(
  rep(
    c(
      "parental",
      "mutant"
    ),
    c(3,3)),
  levels=c("parental","mutant")
)

design<-model.matrix(~isogenics)

run_limma <- function(x, design, prefix="limma_analysis"){
  fit<-lmFit(x, design)
  fit<-eBayes(fit)
  result<-topTable(
    fit,
    coef="isogenicsmutant",
    adjust="BH",
    number=nrow(x)
  )
  
  write.table(
    result,
    file=paste("limma_", prefix, ".txt", sep=""),
    col.names=TRUE,
    row.names=FALSE,
    sep="\t",
    quote=FALSE
  )
}

wt_screen <- read.table(
  summaryFilePath_1,
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE
)

mutant_screen <- read.table(
  summaryFilePath_2,
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE
)

# Join up the columns we need from the par and
# mut data frames

# join_the compound (eg. FLUOROURACIL_50000) and
# function (e.g. anti-metabolite) to make row ids

screen_rownames <- paste(
  wt_screen[,"GeneID"],
  wt_screen[,"Entrez_gene_ID"],
  wt_screen[,"sublib"],
  sep=":"
)

#
#
#

screen_wt_mutant <- as.matrix(
  cbind(
    wt_screen[,c(ncol(wt_screen)-2,ncol(wt_screen)-1,ncol(wt_screen))],
    mutant_screen[,c(ncol(wt_screen)-2,ncol(wt_screen)-1,ncol(wt_screen))]
  )
)

colnames(screen_wt_mutant) <- c(
  "wt1",
  "wt2",
  "wt3",
  "mutant1",
  "mutant2",
  "mutant3"
)
rownames(screen_wt_mutant) <- screen_rownames

short_summary1 = gsub("_summary_with_rep_zscores.txt","",summaryFilePath_1)
shorter_summary1 = gsub("summaryFile1_","",short_summary1)
short_summary2 = gsub("_summary_with_rep_zscores.txt","",summaryFilePath_2)
shorter_summary2 = gsub("summaryFile2_","",short_summary2)

# 
run_limma(
  na.omit(screen_wt_mutant), # table of data
  design=design,
  prefix=paste(shorter_summary1,"_and_",shorter_summary2,sep="")
)

pdf(
  paste("limma_",shorter_summary1,"_and_",shorter_summary2,".pdf",sep=""),
  width=9,
  height=9
)
plot(
  as.data.frame(screen_wt_mutant),
  col=rgb(0,0,0,0.5),
  pch=19,
  cex=0.5,
  main="wildtype vs. mutant"
)

dev.off()
