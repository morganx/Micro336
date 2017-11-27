---
output: pdf_document
---

========================================================

## Load [the phyloseq package](http://joey711.github.com/phyloseq/)

```r
library("phyloseq")
setwd("~/micro336/")
ptm <- proc.time()
```

### Load the HMPv35 dataset.
Create a temporary local file where there large file will be downloaded and read-from.
#This file downloadable from: https://raw.github.com/joey711/phyloseq-demo/gh-pages/HMPv35.RData but is included in your class data directory

```r
load("data/HMPv35.RData")
HMPv35
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45336 taxa and 4743 samples ]
## sample_data() Sample Data:       [ 4743 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 45336 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 45336 tips and 45099 internal nodes ]
## refseq()      DNAStringSet:      [ 45336 reference sequences ]
```




## Investigate Some Basic Features of HMPv35
For more advanced tools for exploring the data, see [additional phyloseq documentation](http://joey711.github.com/phyloseq/)


```r
ntaxa(HMPv35)
```

```
## [1] 45336
```

```r
nsamples(HMPv35)
```

```
## [1] 4743
```

```r
rank_names(HMPv35)
```

```
## [1] "Rank1"  "Phylum" "Class"  "Order"  "Family" "Genus"
```

```r
head(tax_table(HMPv35))
```

```
## Taxonomy Table:     [6 taxa by 6 taxonomic ranks]:
##              Rank1  Phylum       Class     Order            
## OTU_97.15099 "Root" "Firmicutes" "Bacilli" "Lactobacillales"
## OTU_97.13686 "Root" "Firmicutes" "Bacilli" "Lactobacillales"
## OTU_97.30326 "Root" "Firmicutes" "Bacilli" "Lactobacillales"
## OTU_97.26112 "Root" "Firmicutes" "Bacilli" "Lactobacillales"
## OTU_97.34719 "Root" "Firmicutes" "Bacilli" "Lactobacillales"
## OTU_97.12776 "Root" "Firmicutes" "Bacilli" "Lactobacillales"
##              Family             Genus          
## OTU_97.15099 "Streptococcaceae" "Streptococcus"
## OTU_97.13686 "Streptococcaceae" "Streptococcus"
## OTU_97.30326 "Streptococcaceae" "Streptococcus"
## OTU_97.26112 "Streptococcaceae" "Streptococcus"
## OTU_97.34719 "Streptococcaceae" "Streptococcus"
## OTU_97.12776 "Streptococcaceae" "Streptococcus"
```

```r
sample_variables(HMPv35)
```

```
## [1] "X.SampleID"     "RSID"           "visitno"        "sex"           
## [5] "RUNCENTER"      "HMPbodysubsite" "Mislabeled"     "Contaminated"  
## [9] "Description"
```

```r
levels(sample_data(HMPv35)$HMPbodysubsite)
```

```
##  [1] "Anterior_nares"               "Attached_Keratinized_gingiva"
##  [3] "Buccal_mucosa"                "Hard_palate"                 
##  [5] "Left_Antecubital_fossa"       "Left_Retroauricular_crease"  
##  [7] "Mid_vagina"                   "Palatine_Tonsils"            
##  [9] "Posterior_fornix"             "Right_Antecubital_fossa"     
## [11] "Right_Retroauricular_crease"  "Saliva"                      
## [13] "Stool"                        "Subgingival_plaque"          
## [15] "Supragingival_plaque"         "Throat"                      
## [17] "Tongue_dorsum"                "Vaginal_introitus"
```

```r
HMPv35saliva = subset_samples(HMPv35, HMPbodysubsite == "Saliva")
HMPv35saliva
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45336 taxa and 290 samples ]
## sample_data() Sample Data:       [ 290 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 45336 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 45336 tips and 45099 internal nodes ]
## refseq()      DNAStringSet:      [ 45336 reference sequences ]
```

```r
nsamples(HMPv35saliva)
```

```
## [1] 290
```

```r
sample_names(HMPv35saliva)[1:10]
```

```
##  [1] "700114707" "700111749" "700108534" "700108165" "700114330"
##  [6] "700105883" "700114048" "700114372" "700113014" "700111300"
```

##Lets explore our data

```r
taxa_sums(HMPv35)[1:10]
```

```
## OTU_97.15099 OTU_97.13686 OTU_97.30326 OTU_97.26112 OTU_97.34719 
##          111          616         1469          114          161 
## OTU_97.12776 OTU_97.31377 OTU_97.25805 OTU_97.21215  OTU_97.2726 
##           12           57           60           33          134
```

```r
sample_sums(HMPv35)[1:10]
```

```
## 700114607 700114424 700114712 700114713 700114554 700114606 700114335 
##         2        24      5137      5096      4571      5247      6295 
## 700114707 700114442 700114610 
##      3333      4628      1642
```


###Compare sequences per sample or OTU

```r
library(ggplot2)
readsumsdf = data.frame(nreads = sort(taxa_sums(HMPv35), TRUE), sorted = 1:ntaxa(HMPv35), 
                        type = "OTU")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(HMPv35), 
                                                        TRUE), sorted = 1:nsamples(HMPv35), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p<-p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
print(p)
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Removed 14 rows containing missing values (geom_bar).
```

![plot of chunk Compare sequences per sample or OTU](figure/Compare sequences per sample or OTU-1.png)

###Compare sequences per sample or OTU saliva only

```r
readsumsdf = data.frame(nreads = sort(taxa_sums(HMPv35saliva), TRUE), sorted = 1:ntaxa(HMPv35saliva), 
                        type = "OTU")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(HMPv35saliva), 
                                                        TRUE), sorted = 1:nsamples(HMPv35saliva), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p<-p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
print(p)
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Removed 22380 rows containing missing values (geom_bar).
```

![plot of chunk Compare sequences per sample or OTU saliva](figure/Compare sequences per sample or OTU saliva-1.png)


```r
alpha_summary_H35_sal<- estimate_richness(HMPv35saliva, measures = c("Observed", "Shannon"))
alpha_meta_H35_sal <- data.frame(alpha_summary_H35_sal, sample_data(HMPv35saliva))
p = plot_richness(HMPv35saliva, x = "sex", measures=c("Observed", "Shannon"))
p<-p + geom_boxplot()
print(p)
```

![plot of chunk Calculate alpha diversity [Richness, Shannon] for untrimmed](figure/Calculate alpha diversity [Richness, Shannon] for untrimmed-1.png)

Check that all OTUs have representative counts  
For here taxa = OTU  
__Commands interpretation:__  
_Total number of taxa in dataset:_ sum(taxa_sums(HMPv35saliva) > 0)   

_Any taxa with no hits:_ any(taxa_sums(HMPv35saliva)== 0)


```r
sum(taxa_sums(HMPv35saliva) > 0)
```

```
## [1] 22956
```

```r
any(taxa_sums(HMPv35saliva)== 0)
```

```
## [1] TRUE
```

```r
sum(taxa_sums(HMPv35saliva) == 0)
```

```
## [1] 22380
```

```r
any(taxa_sums(HMPv35saliva) > 1)
```

```
## [1] TRUE
```

```r
sum(taxa_sums(HMPv35saliva) > 1)
```

```
## [1] 18897
```

```r
any(taxa_sums(HMPv35saliva) < 1)
```

```
## [1] TRUE
```

```r
sum(taxa_sums(HMPv35saliva) < 1)
```

```
## [1] 22380
```

#####Prune taxa with less than 1 count and check taxa numbers again

```r
#Create new file with only present (no zeroes) taxa

H35_sal_trim = prune_taxa(taxa_sums(HMPv35saliva) > 1, HMPv35saliva)
any(sample_sums(H35_sal_trim) == 0)
```

```
## [1] FALSE
```

```r
any(sample_sums(H35_sal_trim) > 0)
```

```
## [1] TRUE
```

```r
sum(taxa_sums(H35_sal_trim) > 0)
```

```
## [1] 18897
```

```r
any(sample_sums(H35_sal_trim) < 1)
```

```
## [1] FALSE
```

```r
sum(taxa_sums(H35_sal_trim) < 1)
```

```
## [1] 0
```

###Compare sequences per sample or OTU for trimmed

```r
readsumsdf = data.frame(nreads = sort(taxa_sums(H35_sal_trim), TRUE), sorted = 1:ntaxa(H35_sal_trim), 
                        type = "OTU")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(H35_sal_trim), 
                                                        TRUE), sorted = 1:nsamples(H35_sal_trim), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p<-p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
```
Needs bugfix - these dfs do not join together 
#```{r Calculate alpha diversity [Richness, Shannon] for trimmed }
alpha_summary_H35_sal_trim<- estimate_richness(H35_sal_trim, measures = c("Observed", "Shannon"))
alpha_meta_H35_sal_trim <- data.frame(alpha_summary_H35_sal_trim, sample_data(H35_sal_trim))
p = plot_richness(H35_sal_trim, x = "sex", measures=c("Observed", "Shannon"))
p<-p + geom_boxplot()
print(p)
summary(alpha_summary_H35_sal_trim)
#```



```r
# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(HMPv35saliva))

# Histogram of sample read counts
p<-ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + scale_x_continuous(breaks = seq(0, 30000, by = 5000)) # Ticks from 0-30000, every 5000
print(p)
```

![plot of chunk Plot sequences using a histogram for ease of viewing](figure/Plot sequences using a histogram for ease of viewing-1.png)

```r
summary(sample_sum_df)
```

```
##       sum       
##  Min.   :   57  
##  1st Qu.: 2944  
##  Median : 4024  
##  Mean   : 4944  
##  3rd Qu.: 5856  
##  Max.   :25108
```

#disabled d/t bug
#```{r rarefy and Calculate alpha diversity [Richness, Shannon] for rarified dataset }

H35_sal_rare = rarefy_even_depth(H35_sal_trim, 1000)

alpha_summary_H35_sal_rare<- estimate_richness(H35_sal_rare, measures = c("Observed", "Shannon"))
alpha_meta_H35_sal_rare <- data.frame(alpha_summary_H35_sal_rare, sample_data(H35_sal_trim))
summary(alpha_summary_H35_sal_rare)

p = plot_richness(H35_sal_rare, x = "sex", measures=c("Observed", "Shannon"))
p<-p + geom_boxplot()
print(p)
#calculate t-test to compare trimmed vs rarified results
ttest <- t(sapply(alpha_summary_H35_sal_rare, function(x) unlist(t.test(x~sample_data(H35_sal_rare)$sex)[c("estimate","p.value","statistic","conf.int")])))
ttest

ttest <- t(sapply(alpha_summary_H35_sal_trim, function(x) unlist(t.test(x~sample_data(H35_sal_trim)$sex)[c("estimate","p.value","statistic","conf.int")])))
ttest
#```



```r
#library(vegan)

#trimmed
NMDS.ord <- ordinate(H35_sal_trim, "NMDS", "bray")
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.2012579 
## Run 1 stress 0.2012107 
## ... New best solution
## ... Procrustes: rmse 0.04571728  max resid 0.474442 
## Run 2 stress 0.2035867 
## Run 3 stress 0.2026903 
## Run 4 stress 0.197504 
## ... New best solution
## ... Procrustes: rmse 0.02998204  max resid 0.2781044 
## Run 5 stress 0.2011684 
## Run 6 stress 0.2017195 
## Run 7 stress 0.2009814 
## Run 8 stress 0.2015434 
## Run 9 stress 0.2021713 
## Run 10 stress 0.1994501 
## Run 11 stress 0.1986269 
## Run 12 stress 0.1982034 
## Run 13 stress 0.2009033 
## Run 14 stress 0.1976077 
## ... Procrustes: rmse 0.02677383  max resid 0.2875298 
## Run 15 stress 0.2011641 
## Run 16 stress 0.1968044 
## ... New best solution
## ... Procrustes: rmse 0.01852908  max resid 0.2667802 
## Run 17 stress 0.2058797 
## Run 18 stress 0.2032119 
## Run 19 stress 0.1976448 
## Run 20 stress 0.2056447 
## *** No convergence -- monoMDS stopping criteria:
##      7: no. of iterations >= maxit
##     13: stress ratio > sratmax
```

```r
sampleplot<-plot_ordination(H35_sal_trim, NMDS.ord, type = "samples", color = "sex") + scale_colour_manual(values=c("black", "blue"))
#stressplot(NMDS.ord)
sampleplot <-sampleplot + geom_point()
#print(sampleplot)
                             
#rarified
NMDS.ord <- ordinate(H35_sal_rare, "NMDS", "bray")
```

```
## Error in ordinate(H35_sal_rare, "NMDS", "bray"): object 'H35_sal_rare' not found
```

```r
sampleplot = plot_ordination(H35_sal_rare, NMDS.ord, type = "samples", color = "sex") + scale_colour_manual(values=c("black", "blue"))
```

```
## Error in plot_ordination(H35_sal_rare, NMDS.ord, type = "samples", color = "sex"): object 'H35_sal_rare' not found
```

```r
print(sampleplot)
```

![plot of chunk compare beta div using NMDS](figure/compare beta div using NMDS-1.png)

```r
#print(stressplot(NMDS.ord))
```


Exploring community composition


```r
library(plyr)
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.4.1
```

```
## Warning in .registerS3method(fin[i, 1], fin[i, 2], fin[i, 3], fin[i, 4], :
## restarting interrupted promise evaluation
```

```
## Warning in get(method, envir = home): restarting interrupted promise
## evaluation
```

```
## Warning in get(method, envir = home): internal error -3 in R_decompress1
```

```
## Error: package or namespace load failed for 'dplyr' in get(method, envir = home):
##  lazy-load database '/Library/Frameworks/R.framework/Versions/3.4/Resources/library/dplyr/R/dplyr.rdb' is corrupt
```

```r
#library(scales)


# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample
H35_sal_rare_over1 <- H35_sal_rare %>%
  tax_glom(taxrank = "Phylum")  %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%      
  filter(Abundance > 0.01)   %>%                      
  arrange(Phylum)
```

```
## Error in H35_sal_rare %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {: could not find function "%>%"
```

```r
d<- ggplot(H35_sal_rare_over1,aes(x=sex,y=Abundance,fill=Phylum))+ geom_bar(stat="identity",position="fill") +  scale_y_continuous()
```

```
## Error in ggplot(H35_sal_rare_over1, aes(x = sex, y = Abundance, fill = Phylum)): object 'H35_sal_rare_over1' not found
```

```r
b<- d+ ylab("Relative Abundance")
```

```
## Error in d + ylab("Relative Abundance"): non-numeric argument to binary operator
```

```r
a<- b+ scale_fill_manual(name="Phylum",values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                                 "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                                                 "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z<- a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y<- z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
print(y)
```

```
## Error in print(y): object 'y' not found
```



```r
HPR.fuso = subset_taxa(H35_sal_rare , Phylum == "Fusobacteria")
```

```
## Error in tax_table(physeq): object 'H35_sal_rare' not found
```

```r
# melt to long format (for ggploting) 
# prune out genera  with no hits in each sample
HPR.fuso.RA <- HPR.fuso %>%
  tax_glom(taxrank = "Genus")  %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.00)   %>%
  arrange(Genus)
```

```
## Error in HPR.fuso %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {: could not find function "%>%"
```

```r
d= ggplot(HPR.fuso.RA,aes(x=visitno,y=Abundance,fill=Genus))+ facet_grid(sex~.) + geom_bar(stat="identity",position="fill") +  scale_y_continuous()
```

```
## Error in ggplot(HPR.fuso.RA, aes(x = visitno, y = Abundance, fill = Genus)): object 'HPR.fuso.RA' not found
```

```r
b= d+ ylab("Relative Abundance")
```

```
## Error in d + ylab("Relative Abundance"): non-numeric argument to binary operator
```

```r
a= b+ scale_fill_manual(name="Genus",values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                                 "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                                                 "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z= a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y= z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
print(y)
```

```
## Error in print(y): object 'y' not found
```


```r
HPR.fuso.AA = subset_taxa(H35_sal_rare , Phylum == "Fusobacteria")
```

```
## Error in tax_table(physeq): object 'H35_sal_rare' not found
```

```r
f<-plot_bar(HPR.fuso, x="visitno", facet_grid = sex~Genus)
```

```
## Error in psmelt(physeq): object 'HPR.fuso' not found
```

```r
print(f)
```

```
## Error in print(f): object 'f' not found
```

Focusing on differences in community

```r
#subset data by saliva and stool
HMPSvS = subset_samples(HMPv35, HMPbodysubsite == "Saliva" | HMPbodysubsite == "Stool")
HMPSvS = prune_taxa(taxa_sums(HMPSvS) > 1, HMPSvS)
HMPSvS
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 28297 taxa and 609 samples ]
## sample_data() Sample Data:       [ 609 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 28297 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 28297 tips and 28112 internal nodes ]
## refseq()      DNAStringSet:      [ 28297 reference sequences ]
```

```r
HMPSvS_rare = rarefy_even_depth(HMPSvS, 1000)
```

```
## You set `rngseed` to FALSE. Make sure you've set & recorded
##  the random seed of your session for reproducibility.
## See `?set.seed`
```

```
## ...
```

```
## 31 samples removedbecause they contained fewer reads than `sample.size`.
```

```
## Up to first five removed samples are:
```

```
## 700114048700106946700105973700105621700111811	
```

```
## ...
```

```
## 5166OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
HMPSvS_rare
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 23131 taxa and 578 samples ]
## sample_data() Sample Data:       [ 578 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 23131 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 23131 tips and 22970 internal nodes ]
## refseq()      DNAStringSet:      [ 23131 reference sequences ]
```

```r
alpha_summary_HMPSvS_rare<- estimate_richness(HMPSvS_rare, measures = c("Observed", "Shannon"))
alpha_meta_HMPSvS_rare <- data.frame(alpha_summary_HMPSvS_rare, sample_data(HMPSvS_rare))
summary(alpha_summary_HMPSvS_rare)
```

```
##     Observed        Shannon     
##  Min.   :128.0   Min.   :3.358  
##  1st Qu.:335.0   1st Qu.:5.081  
##  Median :390.0   Median :5.362  
##  Mean   :393.7   Mean   :5.321  
##  3rd Qu.:452.8   3rd Qu.:5.591  
##  Max.   :659.0   Max.   :6.268
```

```r
p = plot_richness(HMPSvS_rare, x = "HMPbodysubsite", measures=c("Observed", "Shannon"))
p + geom_boxplot()
```

![plot of chunk Resubset based on saliva vs stool and redo analysis](figure/Resubset based on saliva vs stool and redo analysis-1.png)

```r
print(p)
```

![plot of chunk Resubset based on saliva vs stool and redo analysis](figure/Resubset based on saliva vs stool and redo analysis-2.png)

```r
# p = plot_richness(HMPSvS_rare, x = "sex", color = "HMPbodysubsite", measures=c("Observed", "Shannon"))
# p + geom_boxplot() 

NMDS.ord <- ordinate(HMPSvS_rare, "NMDS", "bray")
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.06912995 
## Run 1 stress 0.07023644 
## Run 2 stress 0.07007064 
## Run 3 stress 0.07027506 
## Run 4 stress 0.06980516 
## Run 5 stress 0.06955668 
## ... Procrustes: rmse 0.004917026  max resid 0.02307949 
## Run 6 stress 0.06997057 
## Run 7 stress 0.06983178 
## Run 8 stress 0.06973685 
## Run 9 stress 0.06985235 
## Run 10 stress 0.07111422 
## Run 11 stress 0.07045788 
## Run 12 stress 0.06965872 
## Run 13 stress 0.07032548 
## Run 14 stress 0.06964159 
## Run 15 stress 0.06934456 
## ... Procrustes: rmse 0.001465774  max resid 0.0230952 
## Run 16 stress 0.07035092 
## Run 17 stress 0.0700846 
## Run 18 stress 0.07023031 
## Run 19 stress 0.07059629 
## Run 20 stress 0.06961725 
## ... Procrustes: rmse 0.00462997  max resid 0.02070488 
## *** No convergence -- monoMDS stopping criteria:
##     19: stress ratio > sratmax
##      1: scale factor of the gradient < sfgrmin
```

```r
sampleplot = plot_ordination(HMPSvS_rare, NMDS.ord, type = "samples", color = "HMPbodysubsite") + scale_colour_manual(values=c("black", "blue"))
print(sampleplot)
```

![plot of chunk Resubset based on saliva vs stool and redo analysis](figure/Resubset based on saliva vs stool and redo analysis-3.png)

```r
#print(stressplot(NMDS.ord))
```


```r
#calculate t-test to compare saliva vs stool results
ttest <- t(sapply(alpha_summary_HMPSvS_rare, function(x) unlist(t.test(x~sample_data(HMPSvS_rare)$HMPbodysubsite)[c("estimate","p.value","statistic","conf.int")])))
ttest
```

```
##          estimate.mean in group Saliva estimate.mean in group Stool
## Observed                    442.956044                     349.5803
## Shannon                       5.508005                       5.1529
##               p.value statistic.t  conf.int1   conf.int2
## Observed 4.558908e-50    16.46955 82.2396668 104.5117654
## Shannon  5.212384e-29    11.83584  0.2961742   0.4140361
```


```r
#confirm using ANOSIM

HMPbodysubsite_group = get_variable(HMPSvS_rare, "HMPbodysubsite")
HMPbodysubsite_ano = anosim(distance(HMPSvS_rare, "bray"), HMPbodysubsite_group)
```

```
## Error in anosim(distance(HMPSvS_rare, "bray"), HMPbodysubsite_group): could not find function "anosim"
```

```r
HMPbodysubsite_ano$signif
```

```
## Error in eval(expr, envir, enclos): object 'HMPbodysubsite_ano' not found
```

```r
HMPbodysubsite_ano$statistic
```

```
## Error in eval(expr, envir, enclos): object 'HMPbodysubsite_ano' not found
```


```r
#confirm using ADONIS

df = as(sample_data(HMPSvS_rare), "data.frame")
d = distance(HMPSvS_rare, "bray")
bodysiteadonis = adonis(d ~ HMPbodysubsite, df)
```

```
## Error in adonis(d ~ HMPbodysubsite, df): could not find function "adonis"
```

```r
bodysiteadonis
```

```
## Error in eval(expr, envir, enclos): object 'bodysiteadonis' not found
```


```r
#SIMPER
HMPSvS_rare_df <- data.frame(otu_table(HMPSvS_rare))
HMPSvS_rare_df_t <- t(HMPSvS_rare_df)
HMPSvS_rare_sd <- data.frame(sample_data(HMPSvS_rare))
(sim <- with(HMPSvS_rare_sd, simper(HMPSvS_rare_df_t, HMPbodysubsite),permutations = 999))
```

```
## Error in simper(HMPSvS_rare_df_t, HMPbodysubsite): could not find function "simper"
```

```r
result <- summary(sim) 
```

```
## Error in summary(sim): object 'sim' not found
```

```r
#It is also useful to know the average dissimilarity between the two treatment communities being compared
#The results showed that the community composition in treatment at 16S rDNA level are 68% different from each other.
lapply(sim, FUN = function(x){x$overall})
```

```
## Error in lapply(sim, FUN = function(x) {: object 'sim' not found
```

```r
#export data
simper <- result$"Saliva_Stool"
```

```
## Error in eval(expr, envir, enclos): object 'result' not found
```

```r
write.table(simper, "/Users/sergiomorales/Downloads/simper_SalivavStool.txt", sep="\t")
```

```
## Error in is.data.frame(x): object 'simper' not found
```


#####To allow quantitiative (visual comparison) I subsetted the full data by significantly different OTUs by SIMPER [only 25% of variance] 


```r
#Rename OTU ID and add to tax_table.

tax_table(HMPSvS_rare) <- cbind(tax_table(HMPSvS_rare), OTU=taxa_names(HMPSvS_rare))

simper_25 <- subset(simper, cumsum < 0.25) 
```

```
## Error in subset(simper, cumsum < 0.25): object 'simper' not found
```

```r
simper_25 <- data.frame(OTU = rownames(simper_25), simper_25)
```

```
## Error in rownames(simper_25): object 'simper_25' not found
```

```r
keepTaxa_25 <- simper_25$OTU
```

```
## Error in eval(expr, envir, enclos): object 'simper_25' not found
```

```r
print(keepTaxa_25)
```

```
## Error in print(keepTaxa_25): object 'keepTaxa_25' not found
```

```r
SIMPERtop25 <- subset_taxa(HMPSvS_rare, OTU %in% keepTaxa_25)
```

```
## Error in OTU %in% keepTaxa_25: object 'keepTaxa_25' not found
```

```r
#manual version of subsetting by name
# SIMPERtop25 <- subset_taxa(HMPSvS_rare, OTU=="OTU_97.42864"| OTU=="OTU_97.45365"| OTU=="OTU_97.43343"| OTU=="OTU_97.40451"| OTU=="OTU_97.44941"| OTU=="OTU_97.40560"| OTU=="OTU_97.43147"| OTU=="OTU_97.44851"| OTU=="OTU_97.43346"| OTU=="OTU_97.40551"| OTU=="OTU_97.45246"| OTU=="OTU_97.44594"| OTU=="OTU_97.42356"| OTU=="OTU_97.37770"| OTU=="OTU_97.39795"| OTU=="OTU_97.40359"| OTU=="OTU_97.45429"| OTU=="OTU_97.40593"| OTU=="OTU_97.43133"| OTU=="OTU_97.43685"| OTU=="OTU_97.29789"| OTU=="OTU_97.40681"| OTU=="OTU_97.55"| OTU=="OTU_97.39526"| OTU=="OTU_97.43012"| OTU=="OTU_97.44958"| OTU=="OTU_97.40474"| OTU=="OTU_97.19587"| OTU=="OTU_97.39258"| OTU=="OTU_97.42626"| OTU=="OTU_97.30909"| OTU=="OTU_97.40494"| OTU=="OTU_97.39182"| OTU=="OTU_97.27044"| OTU=="OTU_97.30314"| OTU=="OTU_97.66"| OTU=="OTU_97.51"| OTU=="OTU_97.42854"| OTU=="OTU_97.43188"| OTU=="OTU_97.45310"| OTU=="OTU_97.42663"| OTU=="OTU_97.44836"| OTU=="OTU_97.43028"| OTU=="OTU_97.158"| OTU=="OTU_97.40329"| OTU=="OTU_97.33295"| OTU=="OTU_97.42838"| OTU=="OTU_97.156"| OTU=="OTU_97.39904"| OTU=="OTU_97.108"| OTU=="OTU_97.41849"| OTU=="OTU_97.40656"| OTU=="OTU_97.159"| OTU=="OTU_97.33358"| OTU=="OTU_97.29875"| OTU=="OTU_97.106"| OTU=="OTU_97.22457"| OTU=="OTU_97.62"| OTU=="OTU_97.39352"| OTU=="OTU_97.15266"| OTU=="OTU_97.20"| OTU=="OTU_97.37990"| OTU=="OTU_97.31213"| OTU=="OTU_97.9131"| OTU=="OTU_97.43881"| OTU=="OTU_97.40838"| OTU=="OTU_97.38187"| OTU=="OTU_97.44427"| OTU=="OTU_97.34259"| OTU=="OTU_97.44693"| OTU=="OTU_97.9764"| OTU=="OTU_97.100"| OTU=="OTU_97.45366"| OTU=="OTU_97.221"| OTU=="OTU_97.43049"| OTU=="OTU_97.39302"| OTU=="OTU_97.38639"| OTU=="OTU_97.40621"| OTU=="OTU_97.38446"| OTU=="OTU_97.30614"| OTU=="OTU_97.274"| OTU=="OTU_97.24387"| OTU=="OTU_97.39937"| OTU=="OTU_97.40251"| OTU=="OTU_97.39509"| OTU=="OTU_97.27432"| OTU=="OTU_97.41593"| OTU=="OTU_97.70"| OTU=="OTU_97.40331"| OTU=="OTU_97.42422"| OTU=="OTU_97.154"| OTU=="OTU_97.2041"| OTU=="OTU_97.5145"| OTU=="OTU_97.43072"| OTU=="OTU_97.40055"| OTU=="OTU_97.43017"| OTU=="OTU_97.40"| OTU=="OTU_97.44564"| OTU=="OTU_97.18196"| OTU=="OTU_97.4632"| OTU=="OTU_97.45328"| OTU=="OTU_97.192"| OTU=="OTU_97.164"| OTU=="OTU_97.40308"| OTU=="OTU_97.3221"| OTU=="OTU_97.39413"| OTU=="OTU_97.35552")

dat_SIMPERtop25 <- psmelt(SIMPERtop25)
```

```
## Error in psmelt(SIMPERtop25): object 'SIMPERtop25' not found
```

```r
dat_SIMPERtop25 <- dat_SIMPERtop25[order(dat_SIMPERtop25$Phylum),]
```

```
## Error in eval(expr, envir, enclos): object 'dat_SIMPERtop25' not found
```


###Plot significantly different OTUs at Phylum level (top 25%)

```r
d= ggplot(dat_SIMPERtop25, aes(x = HMPbodysubsite, y = Abundance, fill = Phylum))+ geom_bar(stat="identity",position="fill") +  scale_y_continuous(labels = percent_format())
```

```
## Error in ggplot(dat_SIMPERtop25, aes(x = HMPbodysubsite, y = Abundance, : object 'dat_SIMPERtop25' not found
```

```r
c= d+ xlab("HMPbodysubsite")
```

```
## Error in d + xlab("HMPbodysubsite"): non-numeric argument to binary operator
```

```r
b= c+ ylab("Abundance")
```

```
## Error in c + ylab("Abundance"): non-numeric argument to binary operator
```

```r
a= b+ scale_fill_manual(name="Phylum",values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#330000", "#CCCCFF", "#3399FF", "#66FFFF", "#B5EAAA","#FFE87C"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z= a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y= z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
y
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```


##Too many OTUs to plot Genera so reduced to only top 16 OTUs representing 10% of variance (cumulative sum)


```r
simper_10 <- subset(simper, cumsum < 0.10) 
```

```
## Error in subset(simper, cumsum < 0.1): object 'simper' not found
```

```r
simper_10 <- data.frame(OTU = rownames(simper_10), simper_10)
```

```
## Error in rownames(simper_10): object 'simper_10' not found
```

```r
keepTaxa_10 <- simper_10$OTU
```

```
## Error in eval(expr, envir, enclos): object 'simper_10' not found
```

```r
print(keepTaxa_10)
```

```
## Error in print(keepTaxa_10): object 'keepTaxa_10' not found
```

```r
SIMPERtop10 <- subset_taxa(HMPSvS_rare, OTU %in% keepTaxa_10)
```

```
## Error in OTU %in% keepTaxa_10: object 'keepTaxa_10' not found
```

```r
#manual subsetting
# SIMPERtop10 <- subset_taxa(HMPSvS_rare, OTU=="OTU_97.42864"| OTU=="OTU_97.45365"| OTU=="OTU_97.43343"| OTU=="OTU_97.40451"| OTU=="OTU_97.44941"| OTU=="OTU_97.40560"| OTU=="OTU_97.43147"| OTU=="OTU_97.44851"| OTU=="OTU_97.43346"| OTU=="OTU_97.40551"| OTU=="OTU_97.45246"| OTU=="OTU_97.44594"| OTU=="OTU_97.42356"| OTU=="OTU_97.37770"| OTU=="OTU_97.39795")

dat_SIMPERtop10 <- psmelt(SIMPERtop10)
```

```
## Error in psmelt(SIMPERtop10): object 'SIMPERtop10' not found
```

```r
dat_SIMPERtop10 <- dat_SIMPERtop10[order(dat_SIMPERtop10$Phylum),]
```

```
## Error in eval(expr, envir, enclos): object 'dat_SIMPERtop10' not found
```


###Plot significantly different OTUs at Genus level (top 10%)

```r
d= ggplot(dat_SIMPERtop10, aes(x = HMPbodysubsite, y = Abundance, fill = Genus))+ geom_bar(stat="identity",position="fill") +  scale_y_continuous(labels = percent_format())
```

```
## Error in ggplot(dat_SIMPERtop10, aes(x = HMPbodysubsite, y = Abundance, : object 'dat_SIMPERtop10' not found
```

```r
c= d+ xlab("HMPbodysubsite")
```

```
## Error in d + xlab("HMPbodysubsite"): non-numeric argument to binary operator
```

```r
b= c+ ylab("Abundance")
```

```
## Error in c + ylab("Abundance"): non-numeric argument to binary operator
```

```r
a= b+ scale_fill_manual(name="Genus",values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#330000", "#CCCCFF", "#3399FF", "#66FFFF", "#B5EAAA","#FFE87C"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z= a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y= z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
y
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

```r
SIMPERtop10_gen_glom <- tax_glom(SIMPERtop10, taxrank = 'Genus')
```

```
## Error in inherits(x, get.component.classes()): object 'SIMPERtop10' not found
```

```r
dat_SIMPERtop10_gen_glom <- psmelt(SIMPERtop10_gen_glom)
```

```
## Error in psmelt(SIMPERtop10_gen_glom): object 'SIMPERtop10_gen_glom' not found
```

```r
dat_SIMPERtop10_gen_glom <- dat_SIMPERtop10_gen_glom[order(dat_SIMPERtop10_gen_glom$Genus),]
```

```
## Error in eval(expr, envir, enclos): object 'dat_SIMPERtop10_gen_glom' not found
```

```r
d= ggplot(dat_SIMPERtop10_gen_glom, aes(x = HMPbodysubsite, y = Abundance, fill = Genus))+ geom_bar(stat="identity",position="fill") +  scale_y_continuous(labels = percent_format())
```

```
## Error in ggplot(dat_SIMPERtop10_gen_glom, aes(x = HMPbodysubsite, y = Abundance, : object 'dat_SIMPERtop10_gen_glom' not found
```

```r
c= d+ xlab("HMPbodysubsite")
```

```
## Error in d + xlab("HMPbodysubsite"): non-numeric argument to binary operator
```

```r
b= c+ ylab("Abundance")
```

```
## Error in c + ylab("Abundance"): non-numeric argument to binary operator
```

```r
a= b+ scale_fill_manual(name="Genus",values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#330000", "#CCCCFF", "#3399FF", "#66FFFF", "#B5EAAA","#FFE87C"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z= a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y= z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
y
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```




```r
library(Rmisc)
```

```
## Error in library(Rmisc): there is no package called 'Rmisc'
```

```r
sse <- summarySE(dat_SIMPERtop10_gen_glom, measurevar="Abundance", groupvars=c("HMPbodysubsite", "Phylum", "Class", "Family", "Genus"))
```

```
## Error in summarySE(dat_SIMPERtop10_gen_glom, measurevar = "Abundance", : could not find function "summarySE"
```

```r
sse
```

```
## Error in eval(expr, envir, enclos): object 'sse' not found
```

```r
d = ggplot(sse, aes(x = HMPbodysubsite, y = Abundance, fill = Genus)) +  
  geom_bar(stat = "identity", width = 0.85) + theme_bw() + theme(legend.position = "right", strip.background = element_blank(), panel.grid.major = element_blank()
                                                                 ,panel.grid.minor = element_blank()
                                                                 ,panel.background = element_blank())
```

```
## Error in ggplot(sse, aes(x = HMPbodysubsite, y = Abundance, fill = Genus)): object 'sse' not found
```

```r
c= d+ xlab("HMPbodysubsite")
```

```
## Error in d + xlab("HMPbodysubsite"): non-numeric argument to binary operator
```

```r
b= c+ ylab("Abundance")
```

```
## Error in c + ylab("Abundance"): non-numeric argument to binary operator
```

```r
a= b+ scale_fill_manual(name="Genus",values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#330000", "#CCCCFF", "#3399FF", "#66FFFF", "#B5EAAA","#FFE87C"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z= a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y= z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
y
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

###Plot OTUs (SIMPER Top 11) at Genus level (absolute abundance)

```r
pd <- position_dodge(0.1) # move dots .01 to the left and right to avoid overlap
d = ggplot(sse, aes(x = HMPbodysubsite, y = Abundance, fill = sort(Genus))) + facet_grid(~Genus)+
  geom_bar(stat = "identity", width = 0.85) + theme_bw() + theme(legend.position = "right", strip.background = element_blank(), panel.grid.major = element_blank()
                                                                 ,panel.grid.minor = element_blank()
                                                                 ,panel.background = element_blank())
```

```
## Error in ggplot(sse, aes(x = HMPbodysubsite, y = Abundance, fill = sort(Genus))): object 'sse' not found
```

```r
c= d+ xlab("HMPbodysubsite") + geom_errorbar(aes(ymin=Abundance-se, ymax=Abundance+se), colour="black", position=pd)
```

```
## Error in d + xlab("HMPbodysubsite"): non-numeric argument to binary operator
```

```r
b= c+ ylab("Abundance")
```

```
## Error in c + ylab("Abundance"): non-numeric argument to binary operator
```

```r
a= b+ scale_fill_manual(name="Genus",values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#330000", "#CCCCFF", "#3399FF", "#66FFFF", "#B5EAAA","#FFE87C"))
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
z= a+ theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 0.5, size=18),
            axis.text.y = element_text(colour = "black", size=18),
            axis.title.y = element_text(face="bold",size=18),
            plot.title = element_text(size = 18),
            legend.title =element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.position="right",
            legend.key.size = unit(0.50, "cm"),
            strip.text.x = element_text(size=18, face="bold"),
            strip.text.y = element_text(size=18, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black"),
            strip.background = element_rect(colour="black"))
```

```
## Error in eval(expr, envir, enclos): object 'a' not found
```

```r
y= z+ guides(fill = guide_legend(reverse = TRUE, ncol=1))
```

```
## Error in eval(expr, envir, enclos): object 'z' not found
```

```r
y
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```



Confirm SIMPER results with Exact test

#####In order to identify important organisms in each Body Site we can identify OTUs with significant differences across sites

```r
#check open packages
(.packages())
```

```
##  [1] "knitr"         "BiocInstaller" "grid"          "RColorBrewer" 
##  [5] "reshape2"      "scales"        "plyr"          "ggplot2"      
##  [9] "edgeR"         "limma"         "phyloseq"      "stats"        
## [13] "graphics"      "grDevices"     "utils"         "datasets"     
## [17] "methods"       "base"
```

```r
##Close all
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()


library("phyloseq")
library("edgeR")
```

```
## Loading required package: limma
```

```
## Warning: package 'limma' was built under R version 3.4.1
```

```
## Error in value[[3L]](cond): Package 'limma' version 3.32.5 cannot be unloaded:
##  Error in unloadNamespace(package) : namespace 'limma' is imported by 'edgeR' so cannot be unloaded
```

```r
library(ggplot2)
```

```
## 
## Attaching package: 'ggplot2'
```

```
## The following object is masked _by_ '.GlobalEnv':
## 
##     scale_fill_discrete
```

```r
library(plyr)
```

```
## -------------------------------------------------------------------------
```

```
## You have loaded plyr after dplyr - this is likely to cause problems.
## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
## library(plyr); library(dplyr)
```

```
## -------------------------------------------------------------------------
```

```r
library(scales)
```

```
## Warning: package 'scales' was built under R version 3.4.1
```

```r
library(reshape2)
library(RColorBrewer)
library(grid)
library(empiricalFDR.DESeq2)
```

```
## Error in library(empiricalFDR.DESeq2): there is no package called 'empiricalFDR.DESeq2'
```

```r
library("DESeq2")
```

```
## Loading required package: S4Vectors
```

```
## Warning: package 'S4Vectors' was built under R version 3.4.1
```

```
## Loading required package: stats4
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
## Error in value[[3L]](cond): Package 'S4Vectors' version 0.14.3 cannot be unloaded:
##  Error in unloadNamespace(package) : namespace 'S4Vectors' is imported by 'XVector', 'IRanges', 'Biostrings' so cannot be unloaded
```

```r
#Convert phyloseq OTU count data into DGEList for edgeR package
##http://joey711.github.io/phyloseq-extensions/edgeR.html

phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, 
              ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method = method)
  # Check for division by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data,\n         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}




dge = phyloseq_to_edgeR(HMPSvS_rare, group = "HMPbodysubsite")
```

```
## Loading required package: edgeR
```

```
## Loading required package: limma
```

```
## Warning: package 'limma' was built under R version 3.4.1
```

```
## Error in DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, : could not find function "DGEList"
```

```r
# Perform binary test
et = exactTest(dge)
```

```
## Error in exactTest(dge): could not find function "exactTest"
```

```r
# Extract values from test results
tt = topTags(et, n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
```

```
## Error in topTags(et, n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue"): could not find function "topTags"
```

```r
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$FDR <= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(HMPSvS_rare)[rownames(sigtab), ], "matrix"))
 
#subset by most significant only (p value <0.05 and at least a 2-fold change in abudnance)                                                                 
sigtab_2fold <- subset(sigtab, FDR < 0.05 & logFC >= 2 | logFC <= -2 & FDR < 0.05)
head(sigtab_2fold)
```

```
##              Rank1         Phylum               Class           Order
## OTU_97.42864  Root     Firmicutes          Clostridia   Clostridiales
## OTU_97.45365  Root     Firmicutes             Bacilli Lactobacillales
## OTU_97.43343  Root  Bacteroidetes         Bacteroidia   Bacteroidales
## OTU_97.40451  Root  Bacteroidetes         Bacteroidia   Bacteroidales
## OTU_97.44941  Root Proteobacteria Gammaproteobacteria  Pasteurellales
## OTU_97.40560  Root  Bacteroidetes         Bacteroidia   Bacteroidales
##                        Family         Genus          OTU     logFC
## OTU_97.42864  Veillonellaceae   Veillonella OTU_97.42864 -4.666820
## OTU_97.45365 Streptococcaceae Streptococcus OTU_97.45365 -4.288329
## OTU_97.43343   Bacteroidaceae   Bacteroides OTU_97.43343  4.042124
## OTU_97.40451   Bacteroidaceae   Bacteroides OTU_97.40451  3.941965
## OTU_97.44941  Pasteurellaceae   Haemophilus OTU_97.44941 -3.869445
## OTU_97.40560   Bacteroidaceae   Bacteroides OTU_97.40560  3.708270
##                logCPM PValue FDR Rank1.1       Phylum.1
## OTU_97.42864 9.441791      0   0    Root     Firmicutes
## OTU_97.45365 9.059983      0   0    Root     Firmicutes
## OTU_97.43343 9.093976      0   0    Root  Bacteroidetes
## OTU_97.40451 8.953752      0   0    Root  Bacteroidetes
## OTU_97.44941 8.784919      0   0    Root Proteobacteria
## OTU_97.40560 8.806963      0   0    Root  Bacteroidetes
##                          Class.1         Order.1         Family.1
## OTU_97.42864          Clostridia   Clostridiales  Veillonellaceae
## OTU_97.45365             Bacilli Lactobacillales Streptococcaceae
## OTU_97.43343         Bacteroidia   Bacteroidales   Bacteroidaceae
## OTU_97.40451         Bacteroidia   Bacteroidales   Bacteroidaceae
## OTU_97.44941 Gammaproteobacteria  Pasteurellales  Pasteurellaceae
## OTU_97.40560         Bacteroidia   Bacteroidales   Bacteroidaceae
##                    Genus.1        OTU.1
## OTU_97.42864   Veillonella OTU_97.42864
## OTU_97.45365 Streptococcus OTU_97.45365
## OTU_97.43343   Bacteroides OTU_97.43343
## OTU_97.40451   Bacteroides OTU_97.40451
## OTU_97.44941   Haemophilus OTU_97.44941
## OTU_97.40560   Bacteroides OTU_97.40560
```

###Plot significantly different Phyla between Stool and Saliva (p less than 0.05, more than FC2)

```r
brew = brewer.pal(6, "Set1")

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_2fold, !is.na(Phylum))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)

sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
ggplot(sigtabgen, aes(x = Phylum, y = logFC)) + geom_point(size = 6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + scale_fill_brewer() + coord_flip()
```

![plot of chunk Plot significantly different Phyla between Stool and Saliva](figure/Plot significantly different Phyla between Stool and Saliva-1.png)

###Plot significantly different Genera between Saliva and Stool between Stool and Saliva (p less than 0.05, more than FC2)

```r
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_2fold, !is.na(Genus))
# Genera order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)

sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size = 6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + scale_fill_brewer() + coord_flip()
```

![plot of chunk Plot significantly different Genera between Saliva and Stool](figure/Plot significantly different Genera between Saliva and Stool-1.png)


#####To allow quantitiative (visual comparison) I subsetted the full data by significantly different OTUs  


```r
#Rename OTU ID and add to tax_table.

# tax_table(HMPSvS_rare) <- cbind(tax_table(HMPSvS_rare), OTU=taxa_names(HMPSvS_rare))
# 
# onlysigOTU <- subset_taxa(HMPSvS_rare, OTU=="OTU_97.42864"| OTU=="OTU_97.45365"| OTU=="OTU_97.43343"| OTU=="OTU_97.40451"| OTU=="OTU_97.44941"| OTU=="OTU_97.43147"| OTU=="OTU_97.40560"| OTU=="OTU_97.43346"| OTU=="OTU_97.44851"| OTU=="OTU_97.40551"| OTU=="OTU_97.45246"| OTU=="OTU_97.44594"| OTU=="OTU_97.42356"| OTU=="OTU_97.37770"| OTU=="OTU_97.40359"| OTU=="OTU_97.39795"| OTU=="OTU_97.45429"| OTU=="OTU_97.40593"| OTU=="OTU_97.43685"| OTU=="OTU_97.43133"| OTU=="OTU_97.29789"| OTU=="OTU_97.40681"| OTU=="OTU_97.55"| OTU=="OTU_97.39526"| OTU=="OTU_97.43012"| OTU=="OTU_97.44958"| OTU=="OTU_97.40474"| OTU=="OTU_97.19587"| OTU=="OTU_97.42626"| OTU=="OTU_97.39258"| OTU=="OTU_97.39182"| OTU=="OTU_97.40494"| OTU=="OTU_97.30909"| OTU=="OTU_97.27044"| OTU=="OTU_97.66"| OTU=="OTU_97.30314"| OTU=="OTU_97.51"| OTU=="OTU_97.42854"| OTU=="OTU_97.42663"| OTU=="OTU_97.43028"| OTU=="OTU_97.45310"| OTU=="OTU_97.43188"| OTU=="OTU_97.44836"| OTU=="OTU_97.40329"| OTU=="OTU_97.39904"| OTU=="OTU_97.42838"| OTU=="OTU_97.156"| OTU=="OTU_97.158")

sigtabgen_OTUs <- data.frame(OTU = rownames(sigtabgen), sigtabgen)
keepTaxa_ET <- sigtabgen_OTUs$OTU
print(keepTaxa_ET)
```

```
##  [1] OTU_97.42864 OTU_97.45365 OTU_97.43343 OTU_97.40451 OTU_97.44941
##  [6] OTU_97.40560 OTU_97.43147 OTU_97.43346 OTU_97.44851 OTU_97.40551
## [11] OTU_97.45246 OTU_97.44594 OTU_97.42356 OTU_97.37770 OTU_97.39795
## [16] OTU_97.40359 OTU_97.40593 OTU_97.45429 OTU_97.43685 OTU_97.29789
## [21] OTU_97.55    OTU_97.43133 OTU_97.40681 OTU_97.44958 OTU_97.39526
## [26] OTU_97.43012 OTU_97.19587 OTU_97.40474 OTU_97.30314 OTU_97.39258
## [31] OTU_97.39182 OTU_97.42626 OTU_97.30909 OTU_97.51    OTU_97.40494
## [36] OTU_97.27044 OTU_97.44836 OTU_97.42838 OTU_97.42854 OTU_97.42663
## [41] OTU_97.45310 OTU_97.43188 OTU_97.41849 OTU_97.43028 OTU_97.158  
## 45 Levels: OTU_97.158 OTU_97.19587 OTU_97.27044 ... OTU_97.55
```

```r
ET_2fold <- subset_taxa(HMPSvS_rare, OTU %in% keepTaxa_ET)

dat_onlysigOTU <- psmelt(ET_2fold)
```

```
## Warning in psmelt(ET_2fold): The rank names: 
## OTU
##  have been renamed to: 
## taxa_OTU
## to avoid conflicts with special phyloseq plot attribute names.
```

```r
dat_onlysigOTU <- dat_onlysigOTU[order(dat_onlysigOTU$Phylum),]
```


###Plot significantly different OTUs at Phylum level (absolute abundance)

```r
p = ggplot(dat_onlysigOTU, aes(x = HMPbodysubsite, y = Abundance, fill = Phylum)) + facet_grid(.~HMPbodysubsite, drop = TRUE, space = "free", scales = "free") +
  geom_bar(stat = "identity", width = 0.85) + theme_bw() + theme(legend.position = "right", strip.background = element_blank())
p
```

![plot of chunk Phylum abundance for significantly different OTUs](figure/Phylum abundance for significantly different OTUs-1.png)

###Plot significantly different OTUs at Genus level (absolute abundance)

```r
p = ggplot(dat_onlysigOTU, aes(x = HMPbodysubsite, y = Abundance, fill = Genus)) + facet_grid(.~HMPbodysubsite, drop = TRUE, space = "free", scales = "free") +
  geom_bar(stat = "identity", width = 0.85) + theme_bw() + theme(legend.position = "right", strip.background = element_blank())
p + scale_fill_manual(values = c("royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "slateblue4", "tan4", "tomato4", "steelblue", "springgreen4", "snow4", "slategray2", "plum1", "yellow", "sienna", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#FFFFFF", "#0000EE"))
```

![plot of chunk Genus abundance for significantly different OTUs](figure/Genus abundance for significantly different OTUs-1.png)




###Plot significant OTUs at Phylum level by HMPbodysubsite (rel abund)

```r
p= ggplot(dat_onlysigOTU, aes(x = HMPbodysubsite, y = Abundance, fill = Phylum)) 

p +
  geom_bar(stat="identity",position="fill") +
  xlab("BODY SITE") +
  ylab("Abudance") +
  scale_fill_manual(values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "slateblue4", "tan4", "tomato4", "steelblue", "springgreen4", "snow4", "slategray2", "plum1", "yellow", "sienna", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#FFFFFF", "#0000EE")) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=12)
  )+
  theme(axis.title.x = element_text(face="bold",size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(size = 18),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 13),
        legend.position="right",
        #Manipulating the facet features
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black")) + # Black rectangle around facet title
  guides(fill = guide_legend(ncol = 2, title.hjust = 0.4))
```

![plot of chunk Pylum rel abundance for significantly different OTUs](figure/Pylum rel abundance for significantly different OTUs-1.png)

###Plot significant OTUs at Genus by HMPbodysubsite (rel abund)

```r
# Make each species its own physeq object
speciesList <- tapply(sample_names(onlysigOTU), get_variable(onlysigOTU, "HMPbodysubsite"), c)
```

```
## Error in sample_names(onlysigOTU): object 'onlysigOTU' not found
```

```r
speciesPhyseq <- lapply(speciesList, prune_samples, onlysigOTU)
```

```
## Error in lapply(speciesList, prune_samples, onlysigOTU): object 'speciesList' not found
```

```r
## For each item in speciesPhyseq find the average of the taxa rows

# Make a list of OTU tables for each species
speciesOTUtable <- lapply(speciesPhyseq,otu_table)
```

```
## Error in lapply(speciesPhyseq, otu_table): object 'speciesPhyseq' not found
```

```r
# Make a list of average OTU table for each species
speciesAvg <- lapply(speciesOTUtable,rowMeans)
```

```
## Error in lapply(speciesOTUtable, rowMeans): object 'speciesOTUtable' not found
```

```r
# Put every column of average otu counts back into a matrix where each column is a 
# species and each row is an OTU
pooledOTUtable = t(do.call(rbind,speciesAvg))
```

```
## Error in do.call(rbind, speciesAvg): object 'speciesAvg' not found
```

```r
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable)
```

```
## Error in row.names(pooledOTUtable): object 'pooledOTUtable' not found
```

```r
# Add in taxonomy info for OTUs
TT = tax_table(onlysigOTU)
```

```
## Error in tax_table(onlysigOTU): object 'onlysigOTU' not found
```

```r
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
```

```
## Error in eval(expr, envir, enclos): object 'TT' not found
```

```r
tdf = data.frame(TT, OTU = taxa_names(onlysigOTU))
```

```
## Error in data.frame(TT, OTU = taxa_names(onlysigOTU)): object 'TT' not found
```

```r
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
```

```
## Error in merge(pooledOTUtable, tdf, by.x = "OTU"): object 'pooledOTUtable' not found
```

```r
str(pOTUtax)
```

```
## Error in str(pOTUtax): object 'pOTUtax' not found
```

```r
# You will need to change the 6 to the number of columns of metadata you have
# Add in a column for total sequences
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:3]))
```

```
## Error in data.frame(pOTUtax, SeqTotal = rowSums(pOTUtax[, 2:3])): object 'pOTUtax' not found
```

```r
# Take only top X OTUs (user defined)
# Change 500 to whatever number of OTUs you want to look at
pOTU = pOTU[order(-pOTU$SeqTotal),]
```

```
## Error in eval(expr, envir, enclos): object 'pOTU' not found
```

```r
pOTUtop = pOTU[1:48,]
```

```
## Error in eval(expr, envir, enclos): object 'pOTU' not found
```

```r
# This calculaton will tell you what percentage of the data you are representing in the plot
sum(pOTUtop$SeqTotal)/sum(pOTU$SeqTotal)
```

```
## Error in eval(expr, envir, enclos): object 'pOTUtop' not found
```

```r
# Plot bar chart of phylum level differences
# Change the 6 to the number of metadata you have AND
# 8 determines which level of OTUs you look at 
pOTU.phylum =pOTUtop[,c(2:3,9)]
```

```
## Error in eval(expr, envir, enclos): object 'pOTUtop' not found
```

```r
melt.phylum = melt(pOTU.phylum,id.vars="Genus")
```

```
## Error in melt(pOTU.phylum, id.vars = "Genus"): object 'pOTU.phylum' not found
```

```r
colnames(melt.phylum)[2]="OTU"
```

```
## Error in colnames(melt.phylum)[2] = "OTU": object 'melt.phylum' not found
```

```r
agg.phylum=aggregate(.~Genus+OTU,melt.phylum,sum)
```

```
## Error in eval(m$data, parent.frame()): object 'melt.phylum' not found
```

```r
##extra colors (, "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066")

ggplot(agg.phylum,aes(x=OTU,y=value,fill=Genus)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(labels = percent_format())+ 
  xlab("Depth") +
  ylab("Relative Abudance") +
  scale_fill_manual(values = c("grey26", "royalblue", "chartreuse3",  "red", "darkorange","cyan2", "darkgreen", "deepskyblue", "mediumorchid3","#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "slateblue4", "tan4", "tomato4", "steelblue", "springgreen4", "snow4", "slategray2", "plum1", "yellow", "sienna", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "blue4", "yellow1", "violetred", "#990000", "#99CC00", "#003300", "#00CCCC", "#9966CC", "#993366", "#990033", "#4863A0", "#000033", "#330000", "#00CC99", "#00FF33", "#00CCFF", "#FF9933", "#660066", "#FF0066", "#FFFFFF", "#0000EE")) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=12)
  )+
  theme(axis.title.x = element_text(face="bold",size=14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 12),
        legend.position="right",
        #Manipulating the facet features
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black")) + # Black rectangle around facet title
  guides(fill = guide_legend(ncol = 2, title.hjust = 0.4))
```

```
## Error in ggplot(agg.phylum, aes(x = OTU, y = value, fill = Genus)): object 'agg.phylum' not found
```

```r
 proc.time() - ptm
```

```
##    user  system elapsed 
## 317.323 113.760 432.866
```

```r
 object.size(x=lapply(ls(), get))  
```

```
## 1943585088 bytes
```

```r
 print(object.size(x=lapply(ls(), get)), units="Mb")
```

```
## 1853.5 Mb
```
