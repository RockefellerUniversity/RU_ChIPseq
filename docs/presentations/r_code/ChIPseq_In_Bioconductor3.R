params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
library(TFBSTools)
library(GSEABase)
knitr::opts_chunk$set(echo = TRUE, tidy = T) 



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides != "yes"){
  cat("# ChIPseq (part 3)

---
"    
  )
  
}



## ----eval=T,echo=T, message=FALSE,messages=FALSE, eval=T, echo=T, warning=FALSE,tidy=FALSE----
library(GenomicRanges)
macsPeaks <- "data/peaks/Mel_1_peaks.xls"
macsPeaks_DF <- read.delim(macsPeaks,comment.char="#")
macsPeaks_GR <- GRanges(seqnames=macsPeaks_DF[,"chr"],
                        IRanges(macsPeaks_DF[,"start"],macsPeaks_DF[,"end"]))
mcols(macsPeaks_GR) <- macsPeaks_DF[,c("abs_summit", "fold_enrichment")]
macsPeaks_GR[1:5,]


## ----eval=T,echo=T, message=FALSE,messages=FALSE, eval=T, echo=T, warning=FALSE,tidy=FALSE----
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
peakAnno <- annotatePeak(macsPeaks_GR, tssRegion=c(-1000, 1000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb="org.Mm.eg.db")


## ----eval=T,echo=T, message=FALSE,messages=FALSE, eval=T, echo=T, warning=FALSE-----------
annotatedPeaksGR <- as.GRanges(peakAnno)
annotatedPeaksDF <- as.data.frame(peakAnno)
annotatedPeaksDF[1:2,]


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Gene Set Enrichment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Gene Set Enrichment

---
"    
  )
  
}



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------------------
annotatedPeaksGR[1,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------------------
annotatedPeaksGR_TSS <- annotatedPeaksGR[
  annotatedPeaksGR$annotation == "Promoter",]
genesWithPeakInTSS <- unique(annotatedPeaksGR_TSS$geneId)
genesWithPeakInTSS[1:2]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T, message = F---------------------
allGeneGR <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
allGeneGR[1:2,]
allGeneIDs <- allGeneGR$gene_id


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
GO_result <- enrichGO(gene = genesWithPeakInTSS, 
                      universe = allGeneIDs,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
GO_result_df <- data.frame(GO_result)
GO_result_df[1:5, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T, fig.height=4, fig.width=8-------
library(enrichplot)
GO_result_plot <- pairwise_termsim(GO_result)
emapplot(GO_result_plot, showCategory = 20)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
library(msigdbr)
msigdbr_collections()


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
library(msigdbr)
msig_t2g <- msigdbr(species = "Mus musculus", 
                    category = "H", 
                    subcategory = NULL)
msig_t2g <- msig_t2g[ , colnames(msig_t2g) %in% c("gs_name", "entrez_gene")]
msig_t2g[1:3, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
hallmark <- enricher(gene = genesWithPeakInTSS, 
                     universe = allGeneIDs,
                     TERM2GENE = msig_t2g)
hallmark_df <- data.frame(hallmark)
hallmark_df[1:3, ]


## ----eval=T,echo=T, warning=FALSE,tidy=T--------------------------------------------------
allGenesForGOseq <- as.integer(allGeneIDs %in% genesWithPeakInTSS)
names(allGenesForGOseq) <- allGeneIDs
allGenesForGOseq[1:3]


## ----include=FALSE------------------------------------------------------------------------
library(goseq)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
library(goseq)
pwf=nullp(allGenesForGOseq,"mm10","knownGene",plot.fit=FALSE)



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE, message = F, tidy=T--------------------
Myc_hallMarks <- goseq(pwf,"mm10","knownGene",
                       gene2cat = data.frame(msig_t2g))

Myc_hallMarks[1:3, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
library(rGREAT)


## ----eval=T,echo=T, eval=T, echo=T,messages=F,message=F,warning=FALSE,tidy=T--------------
great_Job <- submitGreatJob(macsPeaks_GR,species="mm10",version = "3.0.0",request_interval = 1)
availableCategories(great_Job)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE, message =F-----------------------------
great_ResultTable = getEnrichmentTables(great_Job,
                                        category="Regulatory Motifs")
names(great_ResultTable)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------------------
msigProMotifs <- great_ResultTable[["MSigDB Predicted Promoter Motifs"]]
msigProMotifs[1:4,]


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Motifs

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Motifs

---
"    
  )
  
}



## ----echo=TRUE,include=FALSE--------------------------------------------------------------

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10


## ----echo=TRUE,collapse=F-----------------------------------------------------------------

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10


## ----echo=TRUE,collapse=F-----------------------------------------------------------------
macsSummits_GR <- GRanges(seqnames(macsPeaks_GR),
                          IRanges(macsPeaks_GR$abs_summit,
                                  macsPeaks_GR$abs_summit),
                          score=macsPeaks_GR$fold_enrichment)
macsSummits_GR <- resize(macsSummits_GR,100,fix="center")



## ----echo=TRUE,collapse=F-----------------------------------------------------------------
macsSummits_GR


## ----echo=TRUE,collapse=F-----------------------------------------------------------------
peaksSequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10,
                         macsSummits_GR)
names(peaksSequences) <- paste0(seqnames(macsSummits_GR),":",
                                         start(macsSummits_GR),
                                         "-",
                                         end(macsSummits_GR))

peaksSequences[1:2,]


## ----echo=TRUE,collapse=F-----------------------------------------------------------------
writeXStringSet(peaksSequences,file="mycMel_rep1.fa")



## ----echo=TRUE,collapse=F,eval=FALSE------------------------------------------------------
## library(rtracklayer)
## motifGFF <- import("~/Downloads/fimo.gff")


## ----echo=TRUE,collapse=F,eval=FALSE------------------------------------------------------
## motifGFF$Name <- paste0(seqnames(motifGFF),":",
##                         start(motifGFF),"-",end(motifGFF))
## motifGFF$ID <- paste0(seqnames(motifGFF),":",
##                       start(motifGFF),"-",end(motifGFF))
## export.gff3(motifGFF,con="~/Downloads/fimoUpdated.gff")


## ----echo=TRUE,collapse=F,eval=TRUE-------------------------------------------------------
library(JASPAR2020)
JASPAR2020


## ----echo=TRUE,collapse=F,eval=TRUE-------------------------------------------------------
library(TFBSTools)
pfm <- getMatrixByName(JASPAR2020, 
                       name="MYC")
pfm


## ----echo=TRUE,collapse=F,eval=TRUE-------------------------------------------------------
library(motifmatchr)
MycMotifs <- matchMotifs(pfm,
                         macsSummits_GR,BSgenome.Mmusculus.UCSC.mm10, 
                         out = "positions")
MycMotifs


## -----------------------------------------------------------------------------------------
export.bed(MycMotifs[[1]],con = "MycMotifs.bed")

