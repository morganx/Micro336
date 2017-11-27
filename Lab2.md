---
title: "Micro336 Lab 2"
author: "XCM"
date: "7/6/2017"
output:
  pdf_document: default
  html_document: default
  word_document: default
---

This week's questions and your answers:

1. Describe the data you are analyzing today. What was sequenced to generate these samples, and how are these four samples different from one another?



```r
ptm <- proc.time()
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 3.5 (BiocInstaller 1.26.0), ?biocLite for help
```

```r
#My code for Dada2 should go in this box. End this chunk with your "plotQualityProfile" plots
if(!require(ggplot2)){
    install.packages(c("ggplot2"))
    
}
library(ggplot2)
if(!require("dada2")){
biocLite("dada2")
}
```

```
## Loading required package: dada2
```

```
## Loading required package: Rcpp
```

```
## Warning: package 'Rcpp' was built under R version 3.4.1
```

```r
if(!require("Rqc")){
biocLite("Rqc")
}
```

```
## Loading required package: Rqc
```

```
## Warning: package 'Rqc' was built under R version 3.4.1
```

```
## Loading required package: BiocParallel
```

```
## Loading required package: ShortRead
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: Biostrings
```

```
## Warning: package 'Biostrings' was built under R version 3.4.1
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:plyr':
## 
##     rename
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:plyr':
## 
##     desc
```

```
## The following object is masked from 'package:phyloseq':
## 
##     distance
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'XVector'
```

```
## The following object is masked from 'package:plyr':
## 
##     compact
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: Rsamtools
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicRanges
```

```
## Warning: package 'GenomicRanges' was built under R version 3.4.1
```

```
## Loading required package: GenomicAlignments
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:phyloseq':
## 
##     sampleNames
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:plyr':
## 
##     count
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following object is masked from 'package:Biostrings':
## 
##     type
```

```
## The following object is masked from 'package:base':
## 
##     apply
```

```
## 
## Attaching package: 'ShortRead'
```

```
## The following object is masked from 'package:plyr':
## 
##     id
```

```r
if(!require("phyloseq")){
biocLite("phyloseq")
}

library(phyloseq)
library(dada2)
library(Rqc)
setwd("~/micro336/")
```




```r
# First, change our path to the directory containing the fastq files after unzipping so R can find these files

path <- "data" 
qa<-rqc(path=path, pattern=".fq.gz")
```

```
## Warning: closing unused connection 12 (data/RS-R2.fq.gz)
```

```
## Warning: closing unused connection 11 (data/RS-R1.fq.gz)
```

```
## Warning: closing unused connection 10 (data/RE-R2.fq.gz)
```

```
## Warning: closing unused connection 9 (data/RE-R1.fq.gz)
```

```
## Warning: closing unused connection 8 (data/PS-R2.fq.gz)
```

```
## Warning: closing unused connection 7 (data/PS-R1.fq.gz)
```

```
## Warning: closing unused connection 6 (data/PE-R2.fq.gz)
```

```
## Warning: closing unused connection 5 (data/PE-R1.fq.gz)
```

```
## '/private/var/folders/mp/pqjfht1j5vl23kmxqzxmwtr8ktty5t/T/RtmpasAHUn/rqc_report.html' has been created.
```


```r
#If I did this correctly, running the next command should show me the list of files in the data folder that I downloaded at the beginning of class 

list.files(path)
```

```
##  [1] "databases"    "HMPv35.RData" "PE-R1.fq.gz"  "PE-R2.fq.gz" 
##  [5] "PS-R1.fq.gz"  "PS-R2.fq.gz"  "RE-R1.fq.gz"  "RE-R2.fq.gz" 
##  [9] "README"       "RS-R1.fq.gz"  "RS-R2.fq.gz"  "truth.txt"
```

```r
#This piece of code assigns files that end with –R1.fq.gz to the list of forward fastq files (fnFs), and files that end with –R2.fq.gz to the list of reverse fastq  files (FnRs)

fnFs <- sort(list.files(path, pattern="-R1.fq.gz"))
fnRs <- sort(list.files(path, pattern="-R2.fq.gz"))

# What are the actual names of our samples? Let’s parse them out of the files.

sample.names <- sapply(strsplit(fnFs, "-"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
plotQualityProfile(fnFs)
```

![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-1.png)

```r
plotQualityProfile(fnRs)
```

![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-2.png)

```r
#folder <- system.file(package="ShortRead", "path")
#rqc(path = folder, pattern = ".fq.gz")
rqcReadWidthPlot(qa)
```

![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-3.png)

```r
pairs <- unique(perFileInformation(qa)$pair)
for(pair in pairs) 
{  
qa.sub <- subsetByPair(qa, pair)
f<-rqcReadWidthPlot(qa.sub) + theme(text=element_text(size=6), axis.text.x=element_text(angle=90, hjust=1)) + scale_x_discrete( breaks=seq(0, 250, 5))
print(f)
}
```

![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-4.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-5.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-6.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-7.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-8.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-9.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-10.png)![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-11.png)


2. 	Most Illumina sequencing runs for 16S rRNA are either 150 bp x 2, 250 bp x 2, or 300 bp x 2 paired ends. What is the length of paired end reads used for this Illumina run?


3.	Approximately what is the median read quality and quality range at the end of read1 in the raw data?(Just describe the general tendency of all four forward plots, and estimate).


4.	Approximately what is the median read quality and quality range at the end of read2 in the raw data? (Just describe the geeral tendency of all four reverse plots, and estimate )


5.	What is different about the quality of read1 and read2 in the raw data?


6.	Given that higher read quality indicates a lower sequencing error rate, does read1 or read2 in the raw data have a higher error rate?


7.	In the processed data (top plots), what is the median read quality at the ends of read 1 and 2?




```r
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#I'm saving the name "filterandtrim" but we're just filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
```

```
## Creating output directory:data/filtered
```

```r
head(out)
```

```
##             reads.in reads.out
## PE-R1.fq.gz    38536       481
## PS-R1.fq.gz    19113       180
## RE-R1.fq.gz    38537     24234
## RS-R1.fq.gz    19113     13304
```

```r
errF <- learnErrors(filtFs, multithread=TRUE)
```

```
## Initializing error rates to maximum possible estimate.
## Sample 1 - 481 reads in 65 unique sequences.
## Sample 2 - 180 reads in 30 unique sequences.
## Sample 3 - 24234 reads in 2933 unique sequences.
## Sample 4 - 13304 reads in 1199 unique sequences.
##    selfConsist step 2 
##    selfConsist step 3 
##    selfConsist step 4 
## 
## 
## Convergence after  4  rounds.
## Total reads used:  38199
```

```r
errR <- learnErrors(filtRs, multithread=TRUE)
```

```
## Initializing error rates to maximum possible estimate.
## Sample 1 - 481 reads in 61 unique sequences.
## Sample 2 - 180 reads in 20 unique sequences.
## Sample 3 - 24234 reads in 5075 unique sequences.
## Sample 4 - 13304 reads in 2624 unique sequences.
##    selfConsist step 2 
##    selfConsist step 3 
##    selfConsist step 4 
## 
## 
## Convergence after  4  rounds.
## Total reads used:  38199
```

```r
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
```

Infer sequence variants in each sample

```r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/PE_F_filt.fastq.gz
```

```
## Encountered 65 unique sequences from 481 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/PS_F_filt.fastq.gz
```

```
## Encountered 30 unique sequences from 180 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/RE_F_filt.fastq.gz
```

```
## Encountered 2933 unique sequences from 24234 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/RS_F_filt.fastq.gz
```

```
## Encountered 1199 unique sequences from 13304 total sequences read.
```

```r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/PE_R_filt.fastq.gz
```

```
## Encountered 61 unique sequences from 481 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/PS_R_filt.fastq.gz
```

```
## Encountered 20 unique sequences from 180 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/RE_R_filt.fastq.gz
```

```
## Encountered 5075 unique sequences from 24234 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: data/filtered/RS_R_filt.fastq.gz
```

```
## Encountered 2624 unique sequences from 13304 total sequences read.
```

```r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

```
## Sample 1 - 481 reads in 65 unique sequences.
## Sample 2 - 180 reads in 30 unique sequences.
## Sample 3 - 24234 reads in 2933 unique sequences.
## Sample 4 - 13304 reads in 1199 unique sequences.
```

```r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

```
## Sample 1 - 481 reads in 61 unique sequences.
## Sample 2 - 180 reads in 20 unique sequences.
## Sample 3 - 24234 reads in 5075 unique sequences.
## Sample 4 - 13304 reads in 2624 unique sequences.
```

```r
print(dadaFs[[1]])
```

```
## dada-class: object describing DADA2 denoising results
## 12 sample sequences were inferred from 65 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, BAND_SIZE = 16, USE_QUALS = TRUE
```

```r
print(dadaRs[[1]])
```

```
## dada-class: object describing DADA2 denoising results
## 11 sample sequences were inferred from 61 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, BAND_SIZE = 16, USE_QUALS = TRUE
```

```r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```

```
## 470 paired-reads (in 12 unique pairings) successfully merged out of 481 (in 22 pairings) input.
```

```
## 176 paired-reads (in 7 unique pairings) successfully merged out of 180 (in 10 pairings) input.
```

```
## 23877 paired-reads (in 61 unique pairings) successfully merged out of 24234 (in 234 pairings) input.
```

```
## 13175 paired-reads (in 30 unique pairings) successfully merged out of 13304 (in 84 pairings) input.
```

```r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

```
##                                                                                                                                                                                                                                                                                                                                                                                                                                        sequence
## 1              TTAGGAATCTTCCACAATGGGCGCAAGCCTGATGGAGCGACGCCGCGTGAGGGATGAAGGTTTTCGGATCGTAAACCTCTGAATCTGGGACGAAAGAGCCTTAGGGCAGATGACGGTACCAGAGTAATAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTACCCGGAATCACTGGGCGTAAAGGGCGTGTAGGCGGAAATTTAAGTCTGGTTTTAAAGACCGGGGCTCAACCTCGGGGATGGACTGGATACTGGATTTCTTGACCTCTGGAGAGGTAACTGGAATTCCTGGTGTAGCGGTGGAATGCGTAGATACCAGGAGGAACACCAATGGCGAAGGCAAGTTACTGGACAGAAGGTGACGCTGAGGCGCGAAAGTGTGGGGAGCAAACCGG
## 2 TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG
## 3 TGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTCTGAAGAAGGCCTTCGGGTTGTAAAGGACTTTTGTCAGGGAAGAAAAGGCTGTTGCTAATATCAGCGGCTGATGACGGTACCTGAAGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGGGCGCAGACGGTTACTTAAGCAGGATGTGAAATCCCCGGGCTCAACCCGGGAACTGCGTTCTGAACTGGGTGACTCGAGTGTGTCAGAGGGAGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCTCCTGGGACAACACTGACGTTCATGCCCGAAAGCGTGGGTAGCAAACAGG
## 4 TAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACATATGTGTAAGTAACTGTGCACATCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGG
## 5 TAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGTATGAAGAAGGTTTTCGGATCGTAAAGTACTGTTGTTAGAGAAGAACAAGGATAAGAGTAACTGCTTGTCCCTTGACGGTATCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG
## 6      TGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTCGGGTATGGATACCCGTTTGCATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATATCTTGAGTGCAGTTGAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCCTGCTAAGCTGCAACTGACATTGAGGCTCGAAAGTGTGGGTATCAAACAGG
##   abundance forward reverse nmatch nmismatch nindel prefer accept
## 1       112       1       1     84         0      0      2   TRUE
## 2        68       2       2     71         0      0      1   TRUE
## 3        62       3       4     71         0      0      2   TRUE
## 4        39       5       3     71         0      0      2   TRUE
## 5        37       4       5     71         0      0      1   TRUE
## 6        34       8       8     76         0      0      2   TRUE
```
Make merged data into table & remove chimeras

```r
seqtab <- makeSequenceTable(mergers)
```

```
## The sequences being tabled vary in length.
```

```r
dim(seqtab)
```

```
## [1]  4 74
```

```r
table(nchar(getSequences(seqtab)))
```

```
## 
## 403 404 409 411 416 417 422 424 428 429 430 
##   2   7   1   1   4   1   3   6   3  37   9
```

```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

```
## Identified 44 bimeras out of 74 input sequences.
```

```r
dim(seqtab.nochim)
```

```
## [1]  4 30
```

```r
sum(seqtab.nochim)/sum(seqtab)
```

```
## [1] 0.9917768
```
Bookkeeping of how many reads were lost @ pipeline stages:

```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
```

```
##    input filtered denoised merged tabled nonchim
## PE 38536      481      481    470    470     470
## PS 19113      180      180    176    176     176
## RE 38537    24234    24234  23877  23877   23715
## RS 19113    13304    13304  13175  13175   13027
```

```r
silva<-("data/databases/silva_v128/silva_nr_v128_train_set.fa.gz")
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=TRUE)
```

```
## Error: Input/Output
##   no input files found
##   dirPath: data/databases/silva_v128/silva_nr_v128_train_set.fa.gz
##   pattern: character(0)
```

```r
unname(head(taxa))
```

```
## Error in head(taxa): object 'taxa' not found
```

Send to Phyloseq

```r
samples.out <- rownames(seqtab.nochim)
subject <- samples.out
samdf <- data.frame(Subject=subject)
rownames(samdf) <- samples.out

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
```

```
## Error in tax_table(taxa): object 'taxa' not found
```

```r
ps
```

```
## Error in eval(expr, envir, enclos): object 'ps' not found
```

```r
# What do our reads per sample & taxonomic distribution look like?
 plot_bar(ps, x="Subject", fill="Phylum") 
```

```
## Error in psmelt(physeq): object 'ps' not found
```

```r
 ps.tx <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
```

```
## Error in taxa_are_rows(physeq): object 'ps' not found
```

```r
 plot_bar(ps.tx, x="Subject", fill="Phylum") 
```

```
## Error in psmelt(physeq): object 'ps.tx' not found
```

```r
 proc.time() - ptm
```

```
##    user  system elapsed 
## 138.145  12.458 134.246
```

```r
  object.size(x=lapply(ls(), get))  
```

```
## 35704872 bytes
```

```r
 print(object.size(x=lapply(ls(), get)), units="Mb")
```

```
## 34.1 Mb
```
