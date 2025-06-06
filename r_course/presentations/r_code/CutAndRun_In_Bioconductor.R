params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides != "yes"){
  cat("# Cut-And-Run

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Set Up

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Set Up

---
"    
  )
  
}



## ----setwd_introtoR,eval=F----------------------------------------------------------------
## setwd("/PathToMyDownload/RU_Course_template/r_course")
## # e.g. setwd("~/Downloads/Intro_To_R_1Day/r_course")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# ChIP Approaches

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# ChIP Approaches

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Fastq QC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Fastq QC

---
"    
  )
  
}



## ----shortreada,include=FALSE-------------------------------------------------------------
library(ShortRead)



## ----shortread, warning=F, message=F------------------------------------------------------
library(ShortRead)


## ----echo=F,eval=F------------------------------------------------------------------------
## fqSample <- FastqSampler("~/Downloads/SOX9CNR_D0_rep1_R1.fastq.gz",n=10^6)
## fastq <- yield(fqSample)
## 
## writeFastq(fastq,file = "../data/SOX9CNR_D0_rep1_R1_subsample.fastq.gz",mode = "w")
## 


## ----eval=F, echo=F-----------------------------------------------------------------------
## fastq <- readFastq(dirPath = "data/SOX9CNR_D0_rep1_R1_subsample.fastq.gz")


## ----mycRep1Reads,echo=T,eval=F-----------------------------------------------------------
## fqSample <- FastqSampler("~/Downloads/SOX9CNR_D0_rep1_R1.fastq.gz",n=10^6)
## fastq <- yield(fqSample)


## ----eval=T, echo=T-----------------------------------------------------------------------
fastq <- readFastq(dirPath = "data/SOX9CNR_D0_rep1_R1_subsample.fastq.gz")


## ----mycRep1ReadsShortReadQ---------------------------------------------------------------
fastq


## ----mycRep1ReadsAccessor-----------------------------------------------------------------
readSequences <- sread(fastq)
readQuality <- quality(fastq)
readIDs <- id(fastq)
readSequences


## ----mycRep1ReadsQScores------------------------------------------------------------------
readQuality <- quality(fastq)
readQualities <- alphabetScore(readQuality)
readQualities[1:10]


## ----mycRep1ReadsQScoresPlot--------------------------------------------------------------
library(ggplot2)
toPlot <- data.frame(ReadQ=readQualities)
ggplot(toPlot,aes(x=ReadQ))+geom_histogram()+theme_minimal()


## ----mycRep1ReadsAlpFreq------------------------------------------------------------------
readSequences <- sread(fastq)
readSequences_AlpFreq <- alphabetFrequency(readSequences)
readSequences_AlpFreq[1:3,]


## ----mycRep1ReadsAlpFreqSum---------------------------------------------------------------
summed__AlpFreq  <- colSums(readSequences_AlpFreq)
summed__AlpFreq[c("A","C","G","T","N")]


## ----mycRep1ReadsAlpByCycle---------------------------------------------------------------
readSequences_AlpbyCycle <- alphabetByCycle(readSequences)
readSequences_AlpbyCycle[1:4,1:10]


## ----mycRep1ReadsAlpByCyclePlot-----------------------------------------------------------
AFreq <- readSequences_AlpbyCycle["A",]
CFreq <- readSequences_AlpbyCycle["C",]
GFreq <- readSequences_AlpbyCycle["G",]
TFreq <- readSequences_AlpbyCycle["T",]
toPlot <- data.frame(Count=c(AFreq,CFreq,GFreq,TFreq),
                     Cycle=rep(1:40,4),
                     Base=rep(c("A","C","G","T"),each=40))



## ----mycRep1ReadsAlpByCyclePlot2,cache=TRUE,eval=FALSE,dependson="mycRep1ReadsAlpByCyclePlot",fig.height=4,fig.width=8----
## 
## ggplot(toPlot, aes(y=Count,x=Cycle,colour=Base)) + geom_line() +
##   theme_bw()


## ----mycRep1ReadsAlpByCyclePlot3----------------------------------------------------------

ggplot(toPlot,aes(y=Count,x=Cycle,colour=Base)) + geom_line() + ylim(150000,400000) +
  theme_bw()


## ----mycRep1ReadsQByCycle,cache=TRUE,dependson="mycRep1ReadsAlpFreq"----------------------
qualAsMatrix <- as(readQuality,"matrix")
qualAsMatrix[1:2,]


## ----mycRep1ReadsQByCyclePlot,cache=TRUE,dependson="mycRep1ReadsQByCycle",fig.width=8,fig.height=4----
boxplot(qualAsMatrix[1:1000,])


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides != "yes"){
  cat("# CUT&RUN/ATAC (part 2) - Alignment and Peak Calling

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Aligning CUT&RUN and ATACseq reads

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Aligning CUT&RUN and ATACseq reads

---
"    
  )
  
}



## ----fa1, echo=TRUE, message = F, warning=F-----------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10


## ----fa2,cache=FALSE,echo=TRUE------------------------------------------------------------
mainChromosomes <- paste0("chr",c(1:19,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Mmusculus.UCSC.mm10[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
mainChrSeqSet


## ----fa3, echo=TRUE,eval=FALSE------------------------------------------------------------
## writeXStringSet(mainChrSeqSet,
##                 "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")


## ----echo=TRUE,eval=FALSE-----------------------------------------------------------------
## library(Rsubread)
## buildindex("mm10_mainchrs","BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa",
##            memory=8000,
##            indexSplit=TRUE)
## 


## ----echo=F,eval=FALSE, include=F---------------------------------------------------------
## 
## chr18Seq <- BSgenome.Mmusculus.UCSC.mm10[["chr18"]]
## chr18SeqSet <- DNAStringSet(chr18Seq)
## writeXStringSet(chr18SeqSet,
##                 "BSgenome.Mmusculus.UCSC.mm10.chr18.fa")
## buildindex("mm10_chr18","BSgenome.Mmusculus.UCSC.mm10.chr18.fa",
##            memory=8000,
##            indexSplit=TRUE)


## ----eval=FALSE,include=FALSE, echo=FALSE, message = F------------------------------------
## require(ShortRead)
## read1 <- readFastq("~/Desktop/BRC/training/ATAC.Cut-Run.ChIP/pipeline_files/CnR/SOX9CNR_W6_rep1_allFiles/FQ_QC/SOX9CNR_W6_rep1_QC_R1.fastq.gz")
## read2 <- readFastq("~/Desktop/BRC/training/ATAC.Cut-Run.ChIP/pipeline_files/CnR/SOX9CNR_W6_rep1_allFiles/FQ_QC/SOX9CNR_W6_rep1_QC_R2.fastq.gz")
## 
## writeFastq(read1[1:1000,],"data/SOX9CNR_W6_rep1_1K_R1.fastq.gz")
## writeFastq(read2[1:1000,],"data/SOX9CNR_W6_rep1_1K_R2.fastq.gz")
## # id(read2[1:1000,])
## # myRes <- bamQC("~/Downloads/Sorted_ATAC_50K_2.bam")


## ----eval=TRUE, message = F---------------------------------------------------------------
require(ShortRead)

# first 1000 reads in each file
read1 <- readFastq("data/SOX9CNR_W6_rep1_1K_R1.fastq.gz")
read2 <- readFastq("data/SOX9CNR_W6_rep1_1K_R2.fastq.gz")
id(read1)[1:2]
id(read2)[1:2]


## -----------------------------------------------------------------------------------------
read1_toAlign <- "~/Downloads/SOX9CNR_W6_rep1_QC_R1.fastq.gz"
read2_toAlign <- "~/Downloads/SOX9CNR_W6_rep1_QC_R2.fastq.gz"


## ----echo=F,eval=TRUE, warning=F----------------------------------------------------------
library(Rsubread)


## ----align, echo=TRUE,eval=F, warning=F---------------------------------------------------
## myMapped <- align("mm10_mainchrs",
##                   readfile1 = read1_toAlign,
##                   readfile2 = read2_toAlign,
##                   type = "dna",
##                   output_file = "SOX9CNR_W6_rep1.bam",
##                   nthreads = 4,
##                   minFragLength = 0, maxFragLength = 2000)
## 


## ----sortindex, echo=TRUE,eval=FALSE, message=FALSE---------------------------------------
## library(Rsamtools)
## 
## sortBam("SOX9CNR_W6_rep1.bam", "SOX9CNR_W6_rep1_sorted")
## indexBam("SOX9CNR_W6_rep1_sorted.bam")


## ----coverage, echo=TRUE,eval=FALSE-------------------------------------------------------
## forBigWig <- coverage("SOX9CNR_W6_rep1_sorted.bam")
## forBigWig


## ----bw, echo=TRUE,eval=FALSE, message=FALSE----------------------------------------------
## library(rtracklayer)
## export.bw(forBigWig,con="SOX9CNR_W6_rep1.bw")


## ----weightedCover1, echo=TRUE,eval=FALSE-------------------------------------------------
## mappedReads <- idxstatsBam("SOX9CNR_W6_rep1_sorted.bam")
## TotalMapped <- sum(mappedReads[,"mapped"])
## 
## TotalMapped


## ----weightedCover2, echo=TRUE,eval=FALSE-------------------------------------------------
## 
## forBigWig <- coverage("SOX9CNR_W6_rep1_sorted.bam",
##                       weight = (10^6)/TotalMapped)
## forBigWig
## export.bw(forBigWig,con="SOX9CNR_W6_rep1_weighted.bw")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Post alignment processing - CUT&RUN

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Post alignment processing - CUT&RUN

---
"    
  )
  
}



## ----eval=F-------------------------------------------------------------------------------
## BiocManager::install("Herper")
## library(Herper)
## 


## ----makeCondaEnv, echo=T, eval=F---------------------------------------------------------
## dir.create("miniconda")
## macs_paths <- install_CondaTools(tools= c("macs3", "samtools", "bedtools"), env="CnR_analysis", pathToMiniConda = "miniconda")


## ----eval=F, echo=F-----------------------------------------------------------------------
## tempdir2 <- function() {
##     tempDir <- tempdir()
##     if(dir.exists(tempDir)){
##       tempDir <- file.path(tempDir,"rr")
##     }
##     tempDir <- gsub("\\", "/", tempDir, fixed = TRUE)
##     tempDir
## }
## 
## myMiniconda <- file.path(tempdir2(), "Test")
## install_CondaTools(tools=c("macs3", "samtools", "bedtools"), env="CnR_analysis", pathToMiniConda = myMiniconda)
## 


## ----makeCondaEnv2, echo=T, eval=F--------------------------------------------------------
## macs_paths


## ----testCondaEnv, echo=T, eval=F---------------------------------------------------------
## Herper::with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
##                       system("samtools sort -h"))


## ----processData_readingInDatad, echo=TRUE,eval=TRUE,cache=FALSE--------------------------
library(GenomicAlignments)
flags=scanBamFlag(isProperPair = TRUE)



## ----processData_readingInDatas, echo=TRUE,eval=TRUE,cache=FALSE--------------------------
myParam=ScanBamParam(flag=flags,
                   what=c("qname","mapq","isize", "flag"))
myParam



## ----processData_readingInDataa, echo=TRUE,eval=T-----------------------------------------
sortedBAM <- "data/SOX9CNR_W6_rep1_chr18_sorted.bam"
cnrReads <- readGAlignmentPairs(sortedBAM,
                                 param=myParam)
cnrReads[1:2,]



## ----processData_readingInData2, echo=TRUE,eval=TRUE,cache=FALSE--------------------------
read1 <- GenomicAlignments::first(cnrReads)
read2 <- GenomicAlignments::second(cnrReads)
read2[1:2,]


## ----processData_readingInData3, echo=TRUE,eval=TRUE,cache=FALSE--------------------------
read1MapQ <- mcols(read1)$mapq
read2MapQ <- mcols(read2)$mapq
read1MapQ[1:5]


## ----processData_readingInData4, echo=TRUE,eval=TRUE,cache=FALSE--------------------------
read1MapQFreqs <- table(read1MapQ)
read2MapQFreqs <- table(read2MapQ)
read1MapQFreqs
read2MapQFreqs


## ----processData_readingInData5,fig.width=9,fig.height=4,  echo=TRUE,eval=TRUE,cache=FALSE----
library(ggplot2)
toPlot <- data.frame(MapQ=c(names(read1MapQFreqs),names(read2MapQFreqs)),
           Frequency=c(read1MapQFreqs,read2MapQFreqs),
           Read=c(rep("Read1",length(read1MapQFreqs)),rep("Read2",length(read2MapQFreqs))))
toPlot$MapQ <- factor(toPlot$MapQ,levels = unique(sort(as.numeric(toPlot$MapQ))))
ggplot(toPlot,aes(x=MapQ,y=Frequency,fill=MapQ))+
  geom_bar(stat="identity")+
  facet_grid(~Read)


## ----processData_extractingRead1, echo=T,eval=TRUE,cache=FALSE----------------------------
insertSizes <- abs(mcols(read1)$isize)
head(insertSizes)


## ----processData_plottingFrffagmentLengths, echo=TRUE,eval=TRUE,cache=FALSE---------------

fragLenSizes <- table(insertSizes)
fragLenSizes[1:5]



## ----processData_plottingFrdagmentLengths, echo=TRUE,eval=TRUE,cache=FALSE, fig.height=4, fig.width=8----
library(ggplot2)
toPlot <- data.frame(InsertSize=as.numeric(names(fragLenSizes)),
                            Count=as.numeric(fragLenSizes))
fragLenPlot <- ggplot(toPlot,aes(x=InsertSize,y=Count))+geom_line()
fragLenPlot+theme_bw()


## ----processData_plottingFragmentLengths24, echo=TRUE,eval=TRUE,cache=FALSE, fig.height=4.5, fig.width=8----
fragLenPlot+ 
  geom_vline(xintercept = c(120),colour="darkgreen")+theme_bw()



## ----processData_createOpenRegionBAM, echo=TRUE,eval=TRUE,cache=FALSE---------------------
cnrReads_filter <- cnrReads[insertSizes < 1000 & (!mcols(read1)$mapq == 0 | !mcols(read2)$mapq == 0)]



## ----readinBL, echo=T, eval=T-------------------------------------------------------------
library(rtracklayer)
blacklist <- "data/mm10-blacklist.v2.bed"
bl_regions <- rtracklayer::import(blacklist)
bl_regions


## ----removeBl, echo=T, eval=T-------------------------------------------------------------

fragment_spans <- granges(cnrReads_filter)
bl_remove <- overlapsAny(fragment_spans, bl_regions)
table(bl_remove)



## ----writeFilteredBAM, echo=T, eval=F-----------------------------------------------------
## cnrReads_filter_noBL <- cnrReads_filter[!bl_remove]
## cnrReads_unlist <- unlist(cnrReads_filter_noBL)
## names(cnrReads_unlist) <- mcols(cnrReads_unlist)$qname
## 
## filter_bam <- gsub("sorted.bam","filter.bam", basename(sortedBAM))
## rtracklayer::export(cnrReads_unlist, filter_bam,format = "bam")
## 


## ----fixmate, echo=T, eval=F--------------------------------------------------------------
## 
## forPeak_bam <- gsub("_filter.bam", "_forPeak.bam", filter_bam)
## Herper::with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
##                       {
##                         tempBam <- paste0(tempfile(), ".bam")
##                         system(paste("samtools sort", "-n", "-o", tempBam, filter_bam, sep = " "))
##                         system(paste("samtools fixmate", "-m", tempBam, forPeak_bam, sep = " "))
##                       })
## 


## ----bamPrepMacs, echo=T, eval=F----------------------------------------------------------
## 
## forMacs_bam <- gsub("_forPeak.bam", "_macs.bam", forPeak_bam)
## sortBam(forPeak_bam, gsub(".bam", "", forMacs_bam))
## indexBam(forMacs_bam)
## 


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Peak calling - CUT&RUN

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Peak calling - CUT&RUN

---
"    
  )
  
}



## ----callMacs, echo=TRUE,eval=F, warning=F------------------------------------------------
## peaks_name <- gsub("_macs.bam", "", basename(forMacs_bam))
## with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
##                       system2(command="macs3",args =c("callpeak",
##                       "-t", forMacs_bam,
##                       "-f", "BAMPE",
##                       "--outdir", ".",
##                       "-n", peaks_name)))
## 


## ----makeBedpe, echo=TRUE,eval=F, warning=F-----------------------------------------------
## # use the BAM with blacklist reads removed (from MACS section)
## bedpe <- gsub("\\.bam", "\\.bed", forPeak_bam)
## with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
##                 system(paste("bedtools bamtobed -bedpe -i", forPeak_bam, ">", bedpe, sep = " "))
## )
## 


## ----echo=T,eval=F, include = T, warning=F------------------------------------------------
## library(dplyr)
## library(tibble)
## 
## indexFa("BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")
## seqlengths(Rsamtools::FaFile("BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")) %>%
##   as.data.frame() %>%
##   rownames_to_column() %>%
##   write.table(file = "chrom.lengths.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


## ----makeBedgraph, echo=T,eval=F, warning=F-----------------------------------------------
## 
## bedgraph <- gsub("\\.bed", "\\.bedgraph", bedpe)
## with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
##                 system(paste("bedtools genomecov -bg -i", bedpe, "-g", "data/chrom.lengths.txt", ">", bedgraph, sep = " "))
## )


## ----getSEACR, echo=T,eval=F, warning=F---------------------------------------------------
## download.file("https://github.com/FredHutch/SEACR/archive/refs/tags/v1.4-beta.2.zip", destfile = "~/Downloads/SEACR_v1.4.zip")
## unzip("~/Downloads/SEACR_v1.4.zip", exdir = "~/Downloads/SEACR_v1.4" )
## seacr_path <- "~/Downloads/SEACR_v1.4/SEACR-1.4-beta.2/SEACR_1.4.sh"
## system(paste(seacr_path, "-h"))


## ----SEACRpermissions, echo=T,eval=F, warning=F-------------------------------------------
## system(paste("chmod 777", seacr_path))
## system(paste(seacr_path, "-h"))


## ----runSEACR, echo=T,eval=F, warning=F---------------------------------------------------
## seacr_prefix <- gsub("\\.bedgraph", "_top01", bedgraph)
## 
## system(paste(seacr_path, "-b", bedgraph, "-c 0.01", "-n non", "-m  stringent", "-o", seacr_prefix))


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Post alignment processing and peak calling - ATACseq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Post alignment processing and peak calling - ATACseq

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides != "yes"){
  cat("# CUT&RUN/ATAC (part 3) - Quality Control

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Quality control - CUT&RUN

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Quality control - CUT&RUN

---
"    
  )
  
}


## ----eval=F, echo=F-----------------------------------------------------------------------
## mappedReads <- idxstatsBam("../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_W6_rep1_sorted.bam")
## mappedReads <- mappedReads[mappedReads$seqnames %in% paste0("chr", c(1:19, "X", "Y", "M")), ]
## saveRDS(mappedReads, file="data/idxstatsBam_sox9_CnR_W6R1.rds")


## ----mapped1, echo=TRUE,eval=FALSE--------------------------------------------------------
## mappedReads <- idxstatsBam("~/Downloads/SOX9CNR_W6_rep1_sorted.bam")
## mappedReads <- mappedReads[mappedReads$seqnames %in% paste0("chr", c(1:19, "X", "Y", "M")), ]
## TotalMapped <- sum(mappedReads[,"mapped"])
## ggplot(mappedReads,aes(x=seqnames,y=mapped))+
##   geom_bar(stat="identity")+coord_flip()


## ----mapped, echo=FALSE,eval=TRUE,fig.width=4,fig.height=4--------------------------------
mappedReads <- readRDS("data/idxstatsBam_sox9_CnR_W6R1.rds")
TotalMapped <- sum(mappedReads[,"mapped"])
suppressPackageStartupMessages(library(ggplot2))
ggplot(mappedReads,aes(x=seqnames,y=mapped))+geom_bar(stat="identity")+coord_flip()


## ----mycQCdwdwshowL,include=FALSE---------------------------------------------------------
library(ChIPQC)


## ----eval=F-------------------------------------------------------------------------------
## QCresult <- ChIPQCsample(reads="/pathTo/myCnRreads.bam",
##                          genome="mm10",
##                          peaks = "/pathTo/myCnRpeaks.bed",
##                          blacklist = "/pathTo/mm10_Blacklist.bed")


## ----mycQC,cache=F,eval=TRUE, message=F---------------------------------------------------
library(ChIPQC)
blklist <- rtracklayer::import.bed("data/mm10-blacklist.v2.bed")
qc_sox9_rep1 <- ChIPQCsample("data/SOX9CNR_W6_rep1_chr18_sorted.bam",
                             annotation = "mm10",
                             peaks = "data/SOX9CNR_W6_rep1_chr18_peaks.narrowPeak",
                             blacklist = blklist,
                             chromosomes = "chr18")


## ----mycQC2,cache=F,eval=TRUE-------------------------------------------------------------
class(qc_sox9_rep1)


## ----mycQCshow,eval=TRUE------------------------------------------------------------------
qc_sox9_rep1


## ----mycQCshow2,cache=F,eval=FALSE, echo=T------------------------------------------------
## bamsToQC <- c("~/Downloads/SOX9CNR_D0_rep1_sorted.bam",
##               "~/Downloads/SOX9CNR_D0_rep2_sorted.bam",
##               "~/Downloads/SOX9CNR_W6_rep1_sorted.bam",
##               "~/Downloads/SOX9CNR_W6_rep2_sorted.bam")
## 
## peaksToQC <- c("data/peaks/SOX9CNR_D0_rep1_macs_peaks.narrowPeak",
##               "data/peaks/SOX9CNR_D0_rep2_macs_peaks.narrowPeak",
##               "data/peaks/SOX9CNR_W6_rep1_macs_peaks.narrowPeak",
##               "data/peaks/SOX9CNR_W6_rep2_macs_peaks.narrowPeak")
## 


## ----mycQCshow4,cache=F,eval=FALSE, echo=T------------------------------------------------
## 
## myQC <- lapply(seq_along(bamsToQC),function(x){
##   ChIPQCsample(
##     bamsToQC[x],
##     annotation = "mm10",
##     peaks = peaksToQC[x],
##     blacklist = blklist,
##     chromosomes = "chr18"
##   )
## })
## names(myQC) <- basename(bamsToQC)
## 


## ----mycQCshow3,cache=F,eval=FALSE, echo=F------------------------------------------------
## bamsToQC <- c("../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_D0_rep1_sorted.bam",
##               "../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_D0_rep2_sorted.bam",
##               "../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_W6_rep1_sorted.bam",
##               "../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_W6_rep2_sorted.bam")
## 
## peaksToQC <- c("data/peaks/SOX9CNR_D0_rep1_macs_peaks.narrowPeak",
##               "data/peaks/SOX9CNR_D0_rep2_macs_peaks.narrowPeak",
##               "data/peaks/SOX9CNR_W6_rep1_macs_peaks.narrowPeak",
##               "data/peaks/SOX9CNR_W6_rep2_macs_peaks.narrowPeak")
## 
## 
## myQC <- lapply(seq_along(bamsToQC),function(x){
##   ChIPQCsample(
##     bamsToQC[x],
##     annotation = "mm10",
##     peaks = peaksToQC[x],
##     blacklist = blklist,
##     chromosomes = "chr18"
##   )
## })
## names(myQC) <- basename(bamsToQC)
## 
## saveRDS(myQC, "sox9_QC_withPeaks.rds")


## ----qcmetricsA,include=FALSE, echo=F, eval=T---------------------------------------------
myQC <- readRDS("data/sox9_QC_withPeaks.rds")


## ----qcmetrics,cache=FALSE,eval=TRUE------------------------------------------------------
QCmetrics(myQC)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Blacklists and SSD

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Blacklists and SSD

---
"    
  )
  
}



## ----fig.width=6,fig.height=2,warning=FALSE,message=FALSE---------------------------------
plotSSD(myQC)+xlim(0,7)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Library complexity and enrichment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Library complexity and enrichment

---
"    
  )
  
}



## ----fig.width=6,fig.height=3,warning=FALSE,message=FALSE---------------------------------
myFlags <- flagtagcounts(myQC)
myFlags["DuplicateByChIPQC",]/myFlags["Mapped",]


## ----warning=FALSE,message=FALSE,fig.width=8,fig.height=4---------------------------------
p <- plotRegi(myQC)


## ----warning=FALSE,fig.width=12,fig.height=6----------------------------------------------
p

