params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
peakFiles <- dir("data/peaks/",pattern="*.peaks",
                 full.names = TRUE)
peakFiles



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
macsPeaks_DF <- list()
for(i in 1:length(peakFiles)){
  macsPeaks_DF[[i]] <- read.delim(peakFiles[i],
                                  comment.char="#")
}
length(macsPeaks_DF)


## ---- include=FALSE-----------------------------------------------------------
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(tracktables)
library(limma)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
library(GenomicRanges)
macsPeaks_GR <- list()
for(i in 1:length(macsPeaks_DF)){
     peakDFtemp <- macsPeaks_DF[[i]]
     macsPeaks_GR[[i]] <- GRanges(
     seqnames=peakDFtemp[,"chr"],
     IRanges(peakDFtemp[,"start"],
             peakDFtemp[,"end"]
     )
  )
}
macsPeaks_GR[[1]]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
fileNames <- basename(peakFiles)
fileNames
sampleNames <- gsub("_peaks.xls","",fileNames)
sampleNames


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
macsPeaks_GRL <- GRangesList(macsPeaks_GR)
names(macsPeaks_GRL) <- sampleNames
class(macsPeaks_GRL)
names(macsPeaks_GRL)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
lengths(macsPeaks_GRL)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
library(rtracklayer)
macsPeaks_GRLCentred <- resize(macsPeaks_GRL,10,fix="center")
width(macsPeaks_GRLCentred)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
Mel_1_Peaks <- macsPeaks_GRL$Mel_1
Mel_2_Peaks <- macsPeaks_GRL$Mel_2
length(Mel_1_Peaks)
length(Mel_2_Peaks)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
Mel_1_Unique <- Mel_1_Peaks[!Mel_1_Peaks %over% Mel_2_Peaks]
Mel_2_Unique <- Mel_2_Peaks[!Mel_2_Peaks %over% Mel_1_Peaks]
length(Mel_1_Unique)
length(Mel_2_Unique)
export.bed(Mel_1_Unique,"Mel_1_Unique.bed")
export.bed(Mel_2_Unique,"Mel_2_Unique.bed")


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
Mel_1_Common <- Mel_1_Peaks[Mel_1_Peaks %over% Mel_2_Peaks]
Mel_2_Common <- Mel_2_Peaks[Mel_2_Peaks %over% Mel_1_Peaks]
length(Mel_1_Common)
length(Mel_2_Common)
export.bed(Mel_1_Common,"Mel_1_Common.bed")
export.bed(Mel_2_Common,"Mel_2_Common.bed")


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
allPeaksSet_Overlapping <- unlist(macsPeaks_GRL)
allPeaksSet_Overlapping


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
allPeaksSet_nR <- reduce(allPeaksSet_Overlapping)
allPeaksSet_nR
export.bed(allPeaksSet_nR,"allPeaksSet_nR.bed")


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
commonPeaks <- allPeaksSet_nR[allPeaksSet_nR %over% Mel_1_Peaks &
                               allPeaksSet_nR %over% Mel_2_Peaks]
commonPeaks
export.bed(commonPeaks,"commonPeaks.bed")



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
mel1_Only <- allPeaksSet_nR[allPeaksSet_nR %over% Mel_1_Peaks &
                             !allPeaksSet_nR %over% Mel_2_Peaks]
mel2_Only <- allPeaksSet_nR[!allPeaksSet_nR %over% Mel_1_Peaks &
                             allPeaksSet_nR %over% Mel_2_Peaks]
length(mel1_Only)
length(mel2_Only)
export.bed(mel1_Only,"mel1_Only.bed")
export.bed(mel2_Only,"mel2_Only.bed")



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
overlap <- list()
for(i in 1:length(macsPeaks_GRL)){
  overlap[[i]] <- allPeaksSet_nR %over% macsPeaks_GRL[[i]]
}
overlap[[1]][1:2]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
overlapMatrix <- do.call(cbind,overlap)
colnames(overlapMatrix) <- names(macsPeaks_GRL)
overlapMatrix[1:2,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
mcols(allPeaksSet_nR) <- overlapMatrix
allPeaksSet_nR[1:2,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,fig.height=5,fig.width=5----
library(limma)
vennDiagram(mcols(allPeaksSet_nR))


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
vennCounts(mcols(allPeaksSet_nR))


## ----eval=T,echo=T,warning=FALSE----------------------------------------------
ch12_HC_Peaks <- allPeaksSet_nR[rowSums(as.data.frame(mcols(allPeaksSet_nR)[,c("ch12_1","ch12_2")])) >= 2]

export.bed(ch12_HC_Peaks,"ch12_HC_Peaks.bed")

ch12_HC_Peaks[1:2,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
ch12_HC_UniquePeaks <- allPeaksSet_nR[
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("ch12_1","ch12_2")])) >= 2 &
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("Mel_1","Mel_2")])) == 0  
  ]
export.bed(ch12_HC_UniquePeaks,"ch12_HC_UniquePeaks.bed")
ch12_HC_UniquePeaks[1,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
HC_Peaks <- allPeaksSet_nR[
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("ch12_1","ch12_2")])) >= 2 |
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("Mel_1","Mel_2")])) >= 2  
  ]
HC_Peaks
export.bed(HC_Peaks,"HC_Peaks.bed")


## ----eval=F, echo=T,warning=FALSE---------------------------------------------
## 
## library(Rsamtools)
## 
## bams <- c("~/Projects/Results/chipseq/testRun/BAMs/Sorted_Myc_Ch12_1.bam",
##           "~/Projects/Results/chipseq/testRun/BAMs/Sorted_Myc_Ch12_2.bam",
##           "~/Projects/Results/chipseq/testRun/BAMs/Sorted_Myc_Mel_1.bam",
##           "~/Projects/Results/chipseq/testRun/BAMs/Sorted_Myc_Mel_2.bam")
## bamFL <- BamFileList(bams,yieldSize = 5000000)
## bamFL


## ----eval=F, echo=T, warning=FALSE--------------------------------------------
## library(GenomicAlignments)
## myMycCounts <- summarizeOverlaps(HC_Peaks,
##                               reads = bamFL,
##                               ignore.strand = TRUE)
## class(myMycCounts)
## save(myMycCounts,file="data/MycCounts.RData")


## ---- eval=T, echo=F, warning=FALSE-------------------------------------------
suppressPackageStartupMessages(library(GenomicAlignments))
load("data/MycCounts.RData")
class(myMycCounts)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
metaDataFrame <- data.frame(CellLine=c("Ch12","Ch12","Mel","Mel"))
rownames(metaDataFrame) <- colnames(myMycCounts)
metaDataFrame


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
library(DESeq2)
deseqMyc <- DESeqDataSetFromMatrix(countData = assay(myMycCounts),
                              colData = metaDataFrame,
                              design = ~ CellLine,
                              rowRanges=HC_Peaks)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
deseqMyc <- DESeq(deseqMyc)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
deseqMyc


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
MelMinusCh12 <- results(deseqMyc,
                        contrast = c("CellLine","Mel","Ch12"),
                        format="GRanges")
MelMinusCh12 <- MelMinusCh12[order(MelMinusCh12$pvalue),]
class(MelMinusCh12)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
MelMinusCh12[1,]


## ----eval=T,echo=T, eval=F, echo=T, warning=FALSE-----------------------------
## MelMinusCh12Filt <- MelMinusCh12[!is.na(MelMinusCh12$pvalue) | !is.na(MelMinusCh12$padj)]
## UpinMel <-  MelMinusCh12[MelMinusCh12$padj < 0.05 & MelMinusCh12$log2FoldChange > 0]
## DowninMel <-  MelMinusCh12[MelMinusCh12$padj < 0.05 & MelMinusCh12$log2FoldChange < 0]
## export.bed(UpinMel,"UpinMel.bed")
## export.bed(DowninMel,"DowninMel.bed")


## ----eval=T,echo=T, eval=F, echo=T, warning=FALSE-----------------------------
## library(tracktables)
## myReport <- makebedtable(MelMinusCh12Filt,"MelMinusCh12.html",
##                          basedirectory = getwd())
## 
## browseURL(myReport)

