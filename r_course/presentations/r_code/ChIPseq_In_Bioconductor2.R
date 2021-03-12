params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)




## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# ChIPseq (part 2)

---
"    
  )
  
}



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Quality Control

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Quality Control

---
"    
  )
  
}



## ----mycQCdwdwshowL,include=FALSE---------------------------------------------
library(ChIPQC)


## ----eval=F-------------------------------------------------------------------
## QCresult <- ChIPQCsample(reads="/pathTo/myChIPreads.bam",
##                          genome="mm10",
##                          blacklist = "/pathTo/mm10_Blacklist.bed")


## ----mycQC,cache=TRUE,eval=FALSE----------------------------------------------
## library(ChIPQC)
## toBlkList<-"~/Downloads/ENCFF547MET.bed.gz"
## chipqc_MycMel_rep1 <- ChIPQCsample("SR_Myc_Mel_rep1.bam",
##                          annotation = "mm10",
##                          blacklist = toBlkList,
##                          chromosomes = paste0("chr",1:10))
## class(chipqc_MycMel_rep1)
## 


## ----mycQCsecret,eval=FALSE,echo=F--------------------------------------------
## library(ChIPQC)
## toBlkList<-"~/Documents/Box Sync/RU/Teaching/Compilation/Genomes_And_Datasets/mm10/ENCFF547MET.bed.gz"
## chipqc_MycMel_rep1 <- ChIPQCsample("SR_Myc_Mel_rep1.bam",
##                          annotation = "mm10",
##                          blacklist = toBlkList,
##                          chromosomes = paste0("chr",1:10))
## save(chipqc_MycMel_rep1,file='~/Documents/Box Sync/RU/Teaching/RU_side/RU_ChIPseq/chipseq/inst/extdata/data/rep1.RData')


## ----mycQCshowLa,echo=FALSE,eval=TRUE-----------------------------------------
toBlkList<-"data/ENCFF547MET.bed.gz"
library(ChIPQC)
load(file='data/rep1.RData')
class(chipqc_MycMel_rep1)


## ----mycQCshow,eval=TRUE------------------------------------------------------
chipqc_MycMel_rep1


## ---- echo=F, eval=F----------------------------------------------------------
## FQ_FILES<-paste0("~/Documents/Box Sync/RU/Teaching/Compilation/Genomes_And_Datasets/ChIPseq_course/",c("ENCFF001NQP.fastq.gz","ENCFF001NQP.fastq.gz","ENCFF001NGC.fastq.gz","ENCFF001NGO.fastq.gz","ENCFF001NCH.fastq.gz","ENCFF001NCF.fastq.gz","ENCFF001NIM.fastq.gz"))
## 
## FQ_NAMES<-c("Myc_Mel_1.bam","Myc_Mel_2.bam","Myc_Ch12_1.bam","Myc_Ch12_2.bam","input_Mel_1.bam","input_Mel_2.bam","input_Ch12_1.bam")
## 
## myMapped <- align("~/Documents/Box Sync/RU/Teaching/Compilation/Genomes_And_Datasets/mm10/mm10_mainchrs",
##                     FQ_FILES,
##                     output_format = "BAM",
##                     output_file = FQ_NAMES,
##                     nthreads = 4)
## 
## library(Rsamtools)
## library(stringr)
## 
## SR_FQ_NAMES<-paste0("SR_",FQ_NAMES)
## SR_FQ_NAMES_1<-paste0("SR_",str_split(FQ_NAMES,".bam", simplify = T)[,1])
## 
## bplapply(1:length(SR_FQ_NAMES), function(x){
## sortBam(FQ_NAMES[x], SR_FQ_NAMES_1[x])
## indexBam(SR_FQ_NAMES[x])
## })
## 


## ----mycQCshowd2,cache=TRUE,eval=FALSE,include=FALSE, echo=F------------------
## FQ_NAMES<-c("Myc_Mel_1.bam","Myc_Mel_2.bam","Myc_Ch12_1.bam","Myc_Ch12_2.bam","input_Mel_1.bam","input_Mel_2.bam","input_Ch12_1.bam")
## SR_FQ_NAMES<-paste0("SR_",FQ_NAMES)
## bamsToQC <- SR_FQ_NAMES
## myQC <- bplapply(bamsToQC,ChIPQCsample,
##         annotation = "mm10",
##         blacklist = toBlkList,
##         chromosomes = paste0("chr",1:10))
## names(myQC)<-bamsToQC
## save(myQC, file="data/myQCnoPeaks.RData")
## # tried to update, but ChIPQC is upset. so leave it for now and owrk with old chipqc object


## ----mycQCshow2,cache=TRUE,eval=FALSE-----------------------------------------
## bamsToQC <- c("Sorted_Myc_Ch12_1.bam","Sorted_Myc_Ch12_2.bam",
##              "Sorted_Myc_MEL_1.bam","Sorted_Myc_MEL_2.bam",
##              "Sorted_Input_MEL.bam","Sorted_Input_Ch12.bam")
## myQC <- bplapply(bamsToQC,ChIPQCsample,
##         annotation = "mm10",
##         blacklist = toBlkList,
##         chromosomes = paste0("chr",1:10))
## names(myQC) <- bamsToQC


## ----qcmetricsA,include=FALSE-------------------------------------------------
load(file="data/myQCnoPeaks.RData")


## ----qcmetrics,cache=FALSE,eval=TRUE------------------------------------------
QCmetrics(myQC)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Assessing fragment length

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Assessing fragment length

---
"    
  )
  
}



## ----qcmetridedecs,cache=FALSE,eval=TRUE,fig.width=6,fig.height=4-------------
plotCC(myQC,facetBy = "Sample")


## ----qcmetridecs,cache=FALSE,eval=TRUE----------------------------------------
myMeta <- data.frame(Sample= names(myQC),
                     Tissue=c("Ch12","Ch12","MEL","MEL","MEL","Ch12"),
                     Antibody=c(rep("Myc",4),rep("Input",2)))
myMeta


## ----qcmetricsede,cache=FALSE,eval=TRUE,fig.width=6,fig.height=3--------------
plotCC(myQC,facetBy = "Tissue",addMetaData = myMeta,
       colourBy="Antibody")


## ----qcmetricsrf,cache=FALSE,eval=TRUE,fig.width=6,fig.height=3---------------
plotCC(myQC,facetBy = "Tissue",addMetaData = myMeta,
       colourBy="Antibody")+theme_bw()+
  ggtitle("ChIPQC results")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
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



## ----fig.width=6,fig.height=2,warning=FALSE,message=FALSE---------------------
plotSSD(myQC)+xlim(0,5)


## ----fig.width=6,fig.height=3,warning=FALSE,message=FALSE---------------------
plotSSD(myQC)+xlim(0.2,0.8)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
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



## ----fig.width=6,fig.height=3,warning=FALSE,message=FALSE---------------------
myFlags <- flagtagcounts(myQC)
myFlags["DuplicateByChIPQC",]/myFlags["Mapped",]


## ----warning=FALSE,message=FALSE,fig.width=8,fig.height=4---------------------
p <- plotRegi(myQC)


## ----warning=FALSE,fig.width=12,fig.height=6----------------------------------
p


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Peak Calling

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Peak Calling

---
"    
  )
  
}



## macs2 callpeak -t Sorted_Myc_MEL_1.bam

##                –name Mel_Rep1

##                –-outdir PeakDirectory

##                -c Sorted_Input_MEL.bam

## 

## ----fig.height=5, fig.width=15,eval=FALSE------------------------------------
## myChIP <- "Sorted_Myc_MEL_1.bam"
## myControl <- "Sorted_Input_MEL.bam"
## 
## macsCommand <- paste0("macs2 callpeak -t ", myChIP,
##                       " -n ", "Mel_Rep1",
##                       " –-outdir ","PeakDirectory",
##                       " -c ", myControl)
## system(macsCommand)


## ----eval=T,echo=T,  warning=FALSE,collapse=T---------------------------------
macsPeaks <- "data/Mel1_peaks.xls"

macsPeaks_DF <- read.delim(macsPeaks)
macsPeaks_DF[1:8,]


## ----eval=T,echo=T,  warning=FALSE,collapse=T---------------------------------
macsPeaks <- "data/Mel1_peaks.xls"

macsPeaks_DF <- read.delim(macsPeaks, comment.char = "#")
macsPeaks_DF[1:2, ]


## ----eval=T,echo=T,  warning=FALSE,collapse=T---------------------------------
library(GenomicRanges)
macsPeaks_GR <- GRanges(seqnames = macsPeaks_DF[, "chr"], IRanges(macsPeaks_DF[, "start"], macsPeaks_DF[, "end"]))
macsPeaks_GR


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
seqnames(macsPeaks_GR)
ranges(macsPeaks_GR)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
mcols(macsPeaks_GR) <- macsPeaks_DF[, c("abs_summit", "fold_enrichment")]
macsPeaks_GR


## ----eval=F,echo=T,  warning=FALSE,collapse=T---------------------------------
## library(rtracklayer)
## macsPeaks_GR <- import("data/Mel1_peaks.narrowPeak", format = "narrowPeak")


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
library(rtracklayer)
blkList <- import.bed(toBlkList)
macsPeaks_GR <- macsPeaks_GR[!macsPeaks_GR %over% blkList] 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Peak Annotation

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Peak Annotation

---
"    
  )
  
}



## ----include=FALSE------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)



## ----eval=F,echo=T, eval=T, echo=T, warning=FALSE,tidy=T,message=FALSE--------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)



## ----eval=T,echo=T, message=FALSE,messages=FALSE, eval=T, echo=T, warning=FALSE----
peakAnno <- annotatePeak(macsPeaks_GR, tssRegion=c(-500, 500), 
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb="org.Mm.eg.db")
class(peakAnno)


## ----eval=T,echo=T, message=F,messages=F, eval=T, echo=T, warning=FALSE,tidy=T----
peakAnno


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
peakAnno_GR <- as.GRanges(peakAnno)
peakAnno_DF <- as.data.frame(peakAnno)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
peakAnno_GR[2,]


## ---- eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T-------
plotAnnoBar(peakAnno)


## ----eval=T,echo=T, eval=F, echo=T, warning=FALSE,fig.height=5, fig.width=15,tidy=T----
## plotDistToTSS(peakAnno)


## ---- eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T-------
upsetplot(peakAnno, vennpie=F)

