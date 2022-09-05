#!/usr/bin/env Rscript
#Júlia Perera
#7/1/2021
#Script to generate QC table from STAR output and Picard RNASeqMetrics

# Run like this:
#Rscript table4QCpresentation.R  QCdir num.lanes paired color fastv output
# if not installed, install it
if(!("QualityGraphs" %in% installed.packages()[,"Package"])){
    #library(devtools)
    #install_github("margenomics/QualityGraphs")
    source("https://raw.githubusercontent.com/margenomics/QualityGraphs/master/R/RNASeq.metrics.R")
}else{library(QualityGraphs)}

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("Please provide required arguments: 
       1) path to ReadMapping folder 
       2) number of lanes 
       3) STAR directory (directory where STAR *Log.final.out files are stored)
       4) RNA Metrics directory (directory where CollectRnaSeqMetrics output files are stored)
       5) name of ouput files (optional)", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  args[5] = "ReadMappingQC"
  
}

# Read arguments
QCDir <- args[1]
lanes=as.numeric(args[2])
R=1 # after mapping R1 and R2 are together
t=R*lanes
starDir=args[3]
rnametricsDir=args[4]

# Set working dir depending on user and machine
setwd(QCDir)

# Load libraries
library(openxlsx)
library(gridExtra)
library(grid)


## create and add a style to the column headers
headerStyle <- createStyle(
  fontSize = 11, fontColour = "#FFFFFF", halign = "center",
  fgFill = "#4F81BD", border = "TopBottom", borderColour = "#4F81BD",
  valign = "center",wrapText = TRUE
)

star=read.STAR.Logs(starDir)
rnametrics=read.RNA_Metrics(rnametricsDir)
rnametrics=rnametrics[match(rownames(star),rownames(rnametrics)),]

RNASeq.metrics.plots(starDir=starDir,rnametricsDir=rnametricsDir,counts=NULL,
                                 resDir=QCDir,picname="RNASeq.metrics",toPNG=F)
  
samples=unique(gsub("_cDNA.*","",rownames(star)))

# Prepare table
df=data.frame(matrix(NA,ncol = ncol(star)+ncol(rnametrics),nrow=length(samples)))
colnames(df)=c(colnames(star),colnames(rnametrics))
rownames(df)=samples

# Fill table with STAR
#=====================
for (i in 1:length(samples)){
  df[i,1:ncol(star)]=colSums(star[(i*t-(t-1)):(i*t),])
  df[i,(ncol(star)+1):ncol(df)]=colMeans(rnametrics[(i*t-(t-1)):(i*t),])
  
}

df$`% UNIQUELY MAPPED`=df$`UNIQUELY MAPPED`/df$`TOTAL READS`
df$`% MULTIMAPPED`=df$`MULTIMAPPED`/df$`TOTAL READS`
df$`% UNMAPPED (too short)`=df$`UNMAPPED (too short)`/df$`TOTAL READS`



# reorder columns
df=df[,c("TOTAL READS", "UNIQUELY MAPPED","% UNIQUELY MAPPED",
         "MULTIMAPPED", "% MULTIMAPPED", "MULTIMAPPED (too many)", 
         "UNMAPPED (too many mm)", "UNMAPPED (too short)", "% UNMAPPED (too short)", 
         "UNMAPPED (other)", "CHIMERIC",  
         "RIBOSOMAL_BASES", "CODING_BASES", "UTR_BASES", "INTRONIC_BASES", 
         "INTERGENIC_BASES", "MRNA_BASES")]


# Save data to file
wb <- createWorkbook()
## Add worksheets
addWorksheet(wb, "ReadMapping QC")
writeData(wb, "ReadMapping QC", df,  rowNames = T)
addStyle(wb, sheet = "ReadMapping QC", headerStyle, rows = 1, cols = 2:(ncol(df)+1), gridExpand = TRUE)
setRowHeights(wb, sheet = 1, rows = 1, heights =50)
saveWorkbook(wb, file.path(QCDir,paste0(args[5],".xlsx")), overwrite = TRUE)
